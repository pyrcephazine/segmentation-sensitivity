"""Run the toy experiments and write a concise results report."""

import argparse
import csv
import importlib.metadata
import json
import os
from dataclasses import asdict, replace
from functools import lru_cache
from pathlib import Path
from time import perf_counter

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

import numpy as np

from config import MODELS, MODEL_LABELS, Parameters
from geometry import Defect, make_toy
from models import solve


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "results"
DEFAULT = Parameters()


@lru_cache(maxsize=192)
def reference(toy, model, parameters, defect=Defect()):
    """Compute intact tumor flow for a matched comparison."""
    return solve(toy, model, parameters, defect).tumor_flow


def cases():
    """Yield the complete experiment battery."""
    for name in ("line", "y"):
        toy = make_toy(name)
        target = 0 if name == "line" else 1
        branches = range(len(toy.branches))
        for model in MODELS:
            yield "intact", toy, model, DEFAULT, Defect()
        for branch in branches:
            for position in np.linspace(0.1, 0.9, 9):
                for model in MODELS:
                    yield "gap_position", toy, model, DEFAULT, Defect("gap", branch, float(position), 0.001)
            for length in (0.25, 0.5, 1, 2, 3, 4):
                for model in MODELS:
                    yield "gap_length", toy, model, DEFAULT, Defect("gap", branch, 0.5, length * 0.001)
            for factor in (0.5, 0.65, 0.8, 0.9, 1.1, 1.2, 1.35, 1.5):
                for model in MODELS:
                    yield "radius", toy, model, DEFAULT, Defect("radius", branch, 0.5, 0.003, factor)
        for resolution in (6, 8, 12, 16, 20):
            parameters = replace(DEFAULT, tissue_resolution=resolution)
            for model in MODELS[2:]:
                for defect in (Defect(), Defect("gap", target)):
                    yield "tissue_mesh", toy, model, parameters, defect
        for resolution in (2, 4, 8, 16, 32):
            parameters = replace(DEFAULT, vessel_resolution=resolution)
            for model in MODELS:
                for defect in (Defect(), Defect("gap", target)):
                    yield "vessel_mesh", toy, model, parameters, defect
        for radius in (0.75, 1.5, 2.5, 4.0, 6.0):
            parameters = replace(DEFAULT, terminal_radius=radius * 0.001)
            for model in MODELS[2:]:
                for defect in (Defect(), Defect("gap", target)):
                    yield "terminal_support", toy, model, parameters, defect
        for factor in (0.1, 0.3, 1.0, 3.0, 10.0):
            parameters = replace(DEFAULT, wall_conductivity=DEFAULT.wall_conductivity * factor)
            for model in ("lateral", "spread", "line_source"):
                for defect in (Defect(), Defect("gap", target)):
                    yield "wall_strength", toy, model, parameters, defect
        for factor in (1.5, 2.0, 3.0, 5.0, 8.0):
            parameters = replace(DEFAULT, spread_radius_factor=factor)
            for defect in (Defect(), Defect("gap", target)):
                yield "spread_width", toy, "spread", parameters, defect
        for factor in (0.5, 0.7, 0.85, 1.0, 1.15, 1.3, 1.5):
            for model in MODELS:
                yield "repair", toy, model, DEFAULT, Defect("bridge", target, 0.5, 0.002, factor)
        for resolution in (6, 10, 14):
            parameters = replace(DEFAULT, tissue_resolution=resolution)
            for model in ("mixed_terminal", "mixed_spread"):
                for defect in (Defect(), Defect("gap", target)):
                    yield "conservative_flux", toy, model, parameters, defect
    toy = make_toy("y")
    for length in (0.25, 0.5, 1.0, 2.0, 3.0):
        for model in MODELS:
            yield "junction_gap", toy, model, DEFAULT, Defect("junction", 0, 0.5, length * 0.001)
    for factor in (0.1, 0.25, 0.5, 0.75, 1.0):
        for model in MODELS:
            yield "extra", toy, model, replace(DEFAULT, healthy_terminal_factor=0.5), Defect("extra", 1, 0.25, 0.001, factor)


def run_case(experiment, toy, model, parameters, defect):
    """Run one case and return its recorded values."""
    start = perf_counter()
    solution = solve(toy, model, parameters, defect)
    if defect.kind in {"gap", "radius", "bridge"}:
        intact_defect = replace(defect, kind="radius", radius_factor=1.0)
    elif defect.kind == "junction":
        intact_defect = replace(defect, kind="junction_intact")
    elif defect.kind == "extra":
        intact_defect = replace(defect, kind="extra_intact", radius_factor=1.0)
    else:
        intact_defect = Defect()
    intact = (solution.tumor_flow if defect.kind == "none"
              else reference(toy, model, parameters, intact_defect))
    return {
        "experiment": experiment, "toy": toy.name, "model": model, "defect": defect.kind,
        "branch": defect.branch, "position_percent": 100 * defect.position,
        "gap_mm": defect.length * 1000, "radius_factor": defect.radius_factor,
        "radius_percent": 100 * (defect.radius_factor - 1),
        "tissue_resolution": parameters.tissue_resolution, "vessel_resolution": parameters.vessel_resolution,
        "terminal_radius_mm": parameters.terminal_radius * 1000,
        "wall_factor": parameters.wall_conductivity / DEFAULT.wall_conductivity,
        "spread_factor": parameters.spread_radius_factor,
        "healthy_terminal_factor": parameters.healthy_terminal_factor,
        "tumor_flow": solution.tumor_flow, "intact_tumor_flow": intact,
        "remaining_percent": 100 * solution.tumor_flow / intact,
        "terminal_tumor_flow": solution.terminal_tumor_flow,
        "lateral_tumor_flow": solution.lateral_tumor_flow, "inlet_flow": solution.inlet_flow,
        "dofs": solution.dofs, "vessel_spacing_mm": 1000 * solution.vessel_spacing,
        "tissue_spacing_mm": 1000 * solution.tissue_spacing, "seconds": perf_counter() - start,
    }


def load_rows(path):
    """Read a completed experiment table."""
    text_fields = {"experiment", "toy", "model", "defect"}
    with path.open() as file:
        return [{k: v if k in text_fields else float(v) for k, v in row.items()}
                for row in csv.DictReader(file)]


def report(rows, figures, elapsed):
    """Write the results and links to every figure."""
    intact = [r for r in rows if r["experiment"] == "intact"]
    mid_gap = {r["model"]: r["remaining_percent"] for r in rows if r["experiment"] == "gap_position"
               and r["toy"] == "y" and r["branch"] == 1 and abs(r["position_percent"] - 50) < 1e-6}
    healthy_gap = next(r["remaining_percent"] for r in rows if r["experiment"] == "gap_position"
                       and r["toy"] == "y" and r["model"] == "terminal" and r["branch"] == 2
                       and abs(r["position_percent"] - 50) < 1e-6)
    contact_mesh = [r["remaining_percent"] for r in rows if r["experiment"] == "tissue_mesh"
                    and r["toy"] == "y" and r["model"] == "contact" and r["defect"] == "gap"]
    lines = ["# Tumor flow sensitivity in a straight vessel and a Y-shaped vessel", "",
             f"Computed **{len(rows):,} cases** using DOLFINx 0.10.0 and networks_fenicsx. "
             "These are illustrative toy problems with a prescribed tumor region, not fitted animal data.", "",
             "## Main findings", "",
             "- At the default tissue spacing, a 1 mm gap halfway along the Y tumor branch leaves essentially zero delivery in "
             f"the vessel-only and terminal-exchange models, {mid_gap['lateral']:.2f}% under wall exchange, "
             f"{mid_gap['spread']:.2f}% under spread exchange, and {mid_gap['contact']:.1f}% under equal pressures.",
             f"- The equal-pressure prediction after that gap ranges from {min(contact_mesh):.1f}% to "
             f"{max(contact_mesh):.1f}% across tissue meshes. Its large residual delivery is therefore "
             "strongly dependent on discretization, rather than evidence of robust vessel repair.",
             f"- Cutting the healthy branch instead raises terminal-model tumor flow by {healthy_gap - 100:.1f}%, "
             "showing that off-target vessels can control the target flow.",
             "- Restoring the missing interval with its original radius recovers the matched intact "
             "solve. The radius sweeps show how quickly this changes when the bridge is too narrow.", "",
             "## Setup", "",
             "Blood enters the upstream end at 80 mmHg. Venous pressure is 10 mmHg. "
             "There is no independent tissue source: vessel exchange supplies the blood bed, "
             "and distributed venous drainage removes blood. The outside of the box has zero normal flow. "
             "Only anatomical outlets exchange through terminal spheres; newly cut vessel ends are closed.", "",
             "The straight vessel is 24 mm long. The Y has an 11 mm parent vessel and two "
             "13.9 mm daughter vessels. A sphere of radius 4.5 mm around the straight vessel tip "
             "or the upper Y tip designates tumor. In the plots, the tumor region is shaded red.", "",
             "$$J_T=\\sum_a w_aQ_a+\\int_\\Lambda w_T(s)j(s)\\,ds,\\qquad "
             "Q_a=G_a(P_a-\\overline p_a).$$", "",
             "Here the weights are the fraction of each source deposited inside the tumor sphere. "
             "The 3D field is unresolved vascular-bed pressure. All percentages compare with the "
             "same model and parameters on the intact geometry, retaining the same vessel partitions "
             "around the defect so that splitting a vessel does not masquerade as damage. Signed exchange is retained.", "",
             "[Toy geometry](figures/toy_geometry.pdf) · [Model equations](figures/model_equations.pdf)", "",
             "| Parameter | Value |", "|---|---:|",
             "| Blood viscosity | 0.0035 Pa s |", "| Straight vessel radius | 0.30 mm |",
             "| Parent and daughter radii | 0.35 mm and 0.25 mm |",
             "| Tissue mobility | $2\\times10^{-9}$ m² Pa⁻¹ s⁻¹ |",
             "| Venous drainage coefficient | $5\\times10^{-6}$ Pa⁻¹ s⁻¹ |",
             "| Wall conductivity | $5\\times10^{-7}$ m Pa⁻¹ s⁻¹ |",
             "| Terminal conductance | $2\\times10^{-11}$ m³ Pa⁻¹ s⁻¹ |",
             "| Terminal sphere radius | 2.5 mm |", "| Default tissue spacing | 2.5 mm |", "",
             "The fixed-permeability lumen model is matched to Poiseuille in the intact vessel; "
             "its permeability is held fixed when radius changes. Spread exchange uses a transverse "
             "uniform disk and the straight-cylinder pressure correction "
             "$\\beta_\\epsilon^{-1}=\\beta^{-1}+[\\log(\\epsilon/r)-1/4]/(2\\pi M)$. "
             "The default disk radius is five vessel radii. Equal pressures are enforced with "
             "one multiplier per fixed 4 mm wall interval, separated at a gap. The line-source model deposits on the centerline "
             "and samples pressure around the circumference; it is a resolution-dependent finite-element "
             "approximation of the singular source, not a smooth-source replacement.", "",
             "## Intact flow", "", "| Model | Straight vessel, mL min⁻¹ | Y-shaped vessel, mL min⁻¹ |",
             "|---|---:|---:|"]
    for model in MODELS:
        values = [next(r["tumor_flow"] for r in intact if r["model"] == model and r["toy"] == toy)
                  for toy in ("line", "y")]
        lines.append(f"| {MODEL_LABELS[model]} | {values[0]:.3f} | {values[1]:.3f} |")
    lines += ["", "[Intact tumor flow](figures/intact_tumor_flow.pdf)", "",
              "The direct vessel model discharges to a prescribed venous pressure; the coupled "
              "models also include resistance in the unresolved blood bed. The lower intact flow "
              "in the terminal model therefore reflects an added downstream load. Distributed "
              "exchange also changes where blood leaves the resolved vessels, so its tumor flow "
              "need not follow the change in total inlet flow.", "",
              "## Gap location and length", "",
              "Gaps are physical intervals removed before vessel meshing. Position is measured "
              "along the affected anatomical vessel, from its upstream end. All detached vessel pieces "
              "and original outlets remain in the solve. Location sweeps use a 1 mm gap; length "
              "sweeps use a gap centered halfway along the affected vessel.", "",
              "| Model | Straight vessel, flow remaining | Y parent gap | Y tumor-branch gap | Y healthy-branch gap |",
              "|---|---:|---:|---:|---:|"]
    for model in MODELS:
        vals = []
        for toy, branch in (("line", 0), ("y", 0), ("y", 1), ("y", 2)):
            row = next(r for r in rows if r["experiment"] == "gap_position" and r["toy"] == toy
                       and r["model"] == model and r["branch"] == branch and abs(r["position_percent"] - 50) < 1e-6)
            vals.append(row["remaining_percent"])
        lines.append(f"| {MODEL_LABELS[model]} | " + " | ".join(f"{0.0 if abs(v) < 1e-8 else v:.1f}%" for v in vals) + " |")
    lines += ["", "[Gap position, straight vessel](figures/gap_position_line.pdf) · "
              "[Gap position, Y-shaped vessel](figures/gap_position_y.pdf)", "",
              "[Gap length, straight vessel](figures/gap_length_line.pdf) · "
              "[Gap length, Y-shaped vessel](figures/gap_length_y.pdf)", "",
              "[Spatial map of gap effects](figures/gap_effect_map.pdf) · "
              "[Gap at the bifurcation](figures/bifurcation_gap.pdf)", "",
              "A bifurcation gap removes the nearby half-gap length from all three incident vessels. "
              "It is a different topology change from removing an interval in only one branch.", "",
              "Tissue-mediated flow after a gap is a different route from an intact arterial lumen. "
              "A nonzero tumor source therefore does not establish successful vessel reconnection. "
              "A healthy-branch defect can increase tumor delivery by changing competition for inflow.", "",
              "## Radius errors and explicit repair", "",
              "Radius errors affect a 3 mm interval. In the coupled models the changed radius also "
              "changes wall circumference and pressure-sampling locations. A bridge fills a 2 mm "
              "missing interval using an estimated radius. It restores the corresponding exchange "
              "terms as well as axial flow; the multiplier and line-source models retain their own "
              "exchange laws on the restored interval.", "",
              "For a concrete comparison, halving radius in the Y tumor branch leaves 46.1% "
              "of intact tumor flow under Poiseuille resistance and 81.0% under fixed lumen permeability. "
              "Both laws were matched on the intact vessel, so this difference isolates radius dependence.", "",
              "[Radius errors, straight vessel](figures/radius_change_line.pdf) · "
              "[Radius errors, Y-shaped vessel](figures/radius_change_y.pdf)", "",
              "[Bridge radius, straight vessel](figures/bridge_radius_line.pdf) · "
              "[Bridge radius, Y-shaped vessel](figures/bridge_radius_y.pdf)", "",
              "[Vessel pressure before and after repair](figures/pressure_along_straight_vessel.pdf)", "",
              "## Mesh resolution", "",
              "Tissue and vessel resolution are varied separately. Physical gap boundaries, radii, "
              "tumor location, source widths, and coefficients stay fixed. The charts show both "
              "absolute intact flow and the flow remaining after a gap; a stable percentage can "
              "otherwise hide a changing intact prediction.", "",
              "[Tissue resolution](figures/tissue_resolution.pdf) · "
              "[Vessel resolution](figures/vessel_resolution.pdf)", "",
              "| Model | Y intact flow at 5 mm tissue spacing | Y intact flow at 1.5 mm tissue spacing |",
              "|---|---:|---:|"]
    for model in MODELS[2:]:
        vals = [next(r["tumor_flow"] for r in rows if r["experiment"] == "tissue_mesh" and r["toy"] == "y"
                     and r["model"] == model and r["defect"] == "none" and r["tissue_resolution"] == n)
                for n in (6, 20)]
        lines.append(f"| {MODEL_LABELS[model]} | {vals[0]:.3f} | {vals[1]:.3f} |")
    lines += ["", "The conservative volume-exchange calculations use shared face fluxes and "
              "cell pressures in the 3D bed. They are included for terminal and spread exchange "
              "rather than applying ordinary volume flux spaces to an exact surface source. "
              "The equal-pressure constraints use fixed physical wall intervals rather than "
              "one multiplier per vessel mesh element.", "",
              f"Despite the relatively small intact-flow change, equal-pressure flow after the gap "
              f"varies from {min(contact_mesh):.1f}% to {max(contact_mesh):.1f}%. This is the clearest "
              "mesh sensitivity in the study. The finite-resistance exchange models consistently "
              "predict much less remaining tumor delivery for this gap.", "",
              "[Conservative flux and pressure formulations](figures/conservative_flux.pdf)", "",
              "## Exchange assumptions and false connections", "",
              "Each source-width or conductivity comparison uses its own intact reference. "
              "Changing terminal source radius also changes overlap with the tumor when the "
              "source becomes wider than the tumor. The false connection joins the daughter "
              "vessels one quarter of the way beyond the bifurcation. For that experiment the healthy "
              "outlet conductance is half the tumor-outlet conductance, so a false connection "
              "can redistribute blood between branches with different loads.", "",
              "Increasing wall conductivity from one tenth to ten times its default value raises "
              "remaining Y tumor flow after a gap from 0.007% to 13.7%. More exchange provides "
              "a stronger route through the blood bed, but does not restore the missing lumen. "
              "With the cylindrical pressure correction, increasing spread width from 1.5 to eight "
              "vessel radii changes remaining flow only from 0.70% to 0.62%.", "",
              "For terminal exchange, increasing source radius from 4 to 6 mm reduces intact "
              "tumor-assigned flow from 3.506 to 1.670 mL min⁻¹. The wider source exceeds the "
              "4.5 mm tumor radius and assigns more blood outside the target. A false connection "
              "as wide as a daughter vessel raises tumor flow by 1.3% in the direct vessel model "
              "but lowers it by 2.3% with wall exchange, which also adds exchange along the false vessel.", "",
              "[Terminal source radius](figures/terminal_source_radius.pdf) · "
              "[Wall conductivity](figures/wall_conductivity.pdf) · "
              "[Spread radius](figures/spread_radius.pdf) · "
              "[Extra connection](figures/extra_connection.pdf)", "",
              "## Data and model definitions", "",
              "- [All numerical results](data/experiments.csv)",
              "- [Parameters and package versions](data/run.json)",
              "- [Model derivations](../notes/tumor_flow_models_review.md)",
              "- [Implementation and reproduction](../README.md)", "",
              "All plots are vector PDFs with LaTeX-rendered text and equations. Numerical "
              "flow values are in mL min⁻¹.", "",
              "Network meshes use [networks_fenicsx](https://github.com/scientificcomputing/networks_fenicsx), "
              "by Jørgen Schartum Dokken and contributors. "
              "DOLFINx assembles the vessel and tissue finite elements; SciPy solves the coupled "
              "sparse systems.", ""]
    (OUTPUT / "sensitivity_results.md").write_text("\n".join(lines))


def update_review(rows):
    """Add the scientific findings and links to the model review."""
    def value(experiment, toy, model, field, **selection):
        return next(r[field] for r in rows if r["experiment"] == experiment and r["toy"] == toy
                    and r["model"] == model and all(r[k] == v if isinstance(v, str) else abs(r[k] - v) < 1e-7
                                                   for k, v in selection.items()))

    wall = value("gap_position", "y", "lateral", "remaining_percent", branch=1, position_percent=50)
    spread = value("gap_position", "y", "spread", "remaining_percent", branch=1, position_percent=50)
    contact = value("gap_position", "y", "contact", "remaining_percent", branch=1, position_percent=50)
    contact_mesh = [r["remaining_percent"] for r in rows if r["experiment"] == "tissue_mesh"
                    and r["toy"] == "y" and r["model"] == "contact" and r["defect"] == "gap"]
    healthy = value("gap_position", "y", "terminal", "remaining_percent", branch=2, position_percent=50)
    hp = value("radius", "y", "poiseuille", "remaining_percent", branch=1, radius_factor=0.5)
    darcy = value("radius", "y", "darcy", "remaining_percent", branch=1, radius_factor=0.5)
    coarse = value("tissue_mesh", "y", "lateral", "tumor_flow", tissue_resolution=6, defect="none")
    fine = value("tissue_mesh", "y", "lateral", "tumor_flow", tissue_resolution=20, defect="none")
    section = f"""## 4. Numerical results

The models were run on a 24 mm straight vessel and a Y-shaped vessel with an 11 mm parent and two 13.9 mm daughters. A sphere of radius 4.5 mm around the straight vessel tip or the upper daughter tip defines tumor. Blood enters at 80 mmHg and drains toward 10 mmHg. There is no independent tissue source: arterial exchange supplies the unresolved blood bed. The outer box has zero normal flow. [Geometry](../results/figures/toy_geometry.pdf). The [complete numerical study](../results/sensitivity_results.md) includes all figures and parameters.

### Gaps and tumor delivery

Gaps were removed at nine positions in each vessel, with lengths from 0.25 to 4 mm. A separate experiment removed the center of the bifurcation. Original outlets were retained; newly cut ends were closed. Each comparison used an intact vessel with the same mesh partitions around the affected interval.

A 1 mm gap halfway along the Y tumor branch left essentially zero tumor delivery under vessel-only or terminal-only flow, **{wall:.2f}%** under wall exchange, **{spread:.2f}%** under spread exchange, and **{contact:.1f}%** under equal mean pressures at the default mesh spacing.

[Gap position and tumor flow](../results/figures/gap_position_y.pdf) · [Gap length](../results/figures/gap_length_y.pdf) · [Spatial map](../results/figures/gap_effect_map.pdf) · [Bifurcation gap](../results/figures/bifurcation_gap.pdf)

Remaining flow can come from the surrounding blood bed rather than a restored lumen. Cutting the healthy branch instead increased terminal-model tumor delivery by **{healthy - 100:.1f}%**: vessels outside the target matter through competition for inflow. The equal-pressure prediction requires the mesh qualification below.

### Radius errors and repair

Radius changes were applied to a 3 mm interval. Halving the radius halfway along the Y tumor branch left **{hp:.1f}%** of the original tumor flow under Poiseuille resistance, compared with **{darcy:.1f}%** when lumen permeability was held fixed. This isolates the different radius dependence; both laws were matched on the intact vessel.

A bridge restored a 2 mm missing interval, and its radius was varied from half to one and a half times the original radius. Using the original radius recovered the matched intact solution; narrower bridges restored connectivity but could still restrict delivery. Coupled bridges restored their exchange terms as well as axial conductance.

[Radius perturbations](../results/figures/radius_change_y.pdf) · [Bridge radius and recovered flow](../results/figures/bridge_radius_y.pdf) · [Pressure before and after repair](../results/figures/pressure_along_straight_vessel.pdf)

### Mesh and exchange assumptions

Tissue spacing was varied from 5 to 1.5 mm, and vessel resolution was varied separately while physical gap boundaries and source widths stayed fixed. Intact Y tumor flow under wall exchange changed from **{coarse:.3f} to {fine:.3f} mL min⁻¹**. In contrast, the equal-pressure model's remaining flow after the tumor-branch gap ranged from **{min(contact_mesh):.1f}% to {max(contact_mesh):.1f}%**. Its apparent gap tolerance is therefore strongly dependent on discretization. Equal-pressure constraints used fixed 4 mm wall intervals, separated across a gap; substantial mesh sensitivity remained. The conservative volume formulation used shared face fluxes for terminal and spread exchange.

Terminal source radius, wall conductivity, and spread radius were also varied. Increasing wall conductivity from one tenth to ten times default raised remaining tumor flow after the gap from **0.007% to 13.7%**. Increasing terminal source radius from 4 to 6 mm lowered intact tumor-assigned flow from **3.506 to 1.670 mL min⁻¹**, as the source extended beyond the 4.5 mm tumor region. These change different assumptions: how readily blood exchanges, how pressure is averaged, and where blood is assigned. A false connection between daughters with unequal outlet loads also redistributed tumor flow; in exchange models it additionally introduced exchange along the false vessel.

[Tissue resolution](../results/figures/tissue_resolution.pdf) · [Vessel resolution](../results/figures/vessel_resolution.pdf) · [Conservative flux formulation](../results/figures/conservative_flux.pdf) · [Source radius](../results/figures/terminal_source_radius.pdf) · [Wall conductivity](../results/figures/wall_conductivity.pdf) · [Spread radius](../results/figures/spread_radius.pdf) · [False connection](../results/figures/extra_connection.pdf)

These results support separate measures of small radius errors, finite connectivity loss, and recovery after repair. They quantify sensitivity of the toy models; they do not establish physiological accuracy or reactive bolus delivery.
"""
    path = ROOT / "notes" / "tumor_flow_models_review.md"
    before, references = path.read_text().split("\n## References", 1)
    if "\n## 4. Numerical results" in before:
        before = before.split("\n## 4. Numerical results", 1)[0]
    destinations = {
        "resistive-graph-derivation": "radius-errors-and-repair",
        "lateral-average-derivation": "gaps-and-tumor-delivery",
        "terminal-exchange-derivation": "gaps-and-tumor-delivery",
        "perfect-contact-derivation": "gaps-and-tumor-delivery",
        "distributed-source-derivation": "mesh-and-exchange-assumptions",
        "lumen-bridge-derivation": "radius-errors-and-repair",
    }
    for derivation, destination in destinations.items():
        marker = f"[Derivation and discretization](#{derivation})."
        link = f" [Numerical results](#{destination})."
        if marker + link not in before:
            before = before.replace(marker, marker + link)
    path.write_text(before.rstrip() + "\n\n" + section + "\n## References" + references)


def main():
    """Run the experiments or regenerate their figures."""
    parser = argparse.ArgumentParser()
    parser.add_argument("--plots-only", action="store_true")
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()
    start = perf_counter()
    (OUTPUT / "data").mkdir(parents=True, exist_ok=True)
    path = OUTPUT / "data" / "experiments.csv"
    if args.plots_only:
        rows = load_rows(path)
    else:
        existing = load_rows(path) if args.resume and path.exists() else []
        rows = existing.copy()
        keys = {(r["experiment"], r["toy"], r["model"], r["defect"], r["branch"], r["position_percent"],
                 r["gap_mm"], r["radius_factor"], r["tissue_resolution"], r["vessel_resolution"],
                 r["terminal_radius_mm"], r["wall_factor"], r["spread_factor"]) for r in existing}
        sequence = list(cases())
        print(f"Running {len(sequence)} cases", flush=True)
        with path.open("a" if existing else "w", newline="") as file:
            writer = None
            for i, case in enumerate(sequence):
                experiment, toy, model, p, d = case
                key = (experiment, toy.name, model, d.kind, d.branch, 100 * d.position, d.length * 1000,
                       d.radius_factor, p.tissue_resolution, p.vessel_resolution, p.terminal_radius * 1000,
                       p.wall_conductivity / DEFAULT.wall_conductivity, p.spread_radius_factor)
                if key in keys:
                    continue
                row = run_case(*case)
                if writer is None:
                    writer = csv.DictWriter(file, fieldnames=list(row))
                    if not existing:
                        writer.writeheader()
                writer.writerow(row)
                rows.append(row)
                if (i + 1) % 20 == 0 or i == len(sequence) - 1:
                    file.flush()
                    print(f"{i + 1}/{len(sequence)}  {toy.name}  {experiment}  {model}  "
                          f"{perf_counter() - start:.0f} seconds", flush=True)
        versions = {name: importlib.metadata.version(name) for name in
                    ("fenics-dolfinx", "networks_fenicsx", "numpy", "scipy", "matplotlib")}
        manifest = {"parameters": asdict(DEFAULT), "versions": versions, "cases": len(rows),
                    "networks_fenicsx_revision": "4d0dd39d6a788185c589c9f4c5017701cd8d27db",
                    "elapsed_seconds": perf_counter() - start}
        (OUTPUT / "data" / "run.json").write_text(json.dumps(manifest, indent=2))
    solutions = {}
    for model in ("poiseuille", "lateral"):
        for label, defect in (("Intact", Defect()), ("Gap", Defect("gap", length=0.002)),
                              ("Restored", Defect("bridge", length=0.002))):
            solutions[(model, label)] = solve(make_toy("line"), model, DEFAULT, defect)
    import figure_manager
    names = figure_manager.make_figures(rows, OUTPUT / "figures", solutions)
    report(rows, names, perf_counter() - start)
    update_review(rows)
    print(f"Saved {len(names)} PDF figures and {OUTPUT / 'sensitivity_results.md'}", flush=True)


if __name__ == "__main__":
    main()
