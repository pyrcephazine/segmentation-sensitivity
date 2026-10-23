"""One theme and one PDF writer for all experiment figures."""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle

from config import MODEL_EQUATIONS, MODEL_LABELS, MODELS
from geometry import make_toy


COLORS = {
    "poiseuille": "#222222", "darcy": "#777777", "lateral": "#0072B2",
    "terminal": "#D55E00", "contact": "#009E73", "spread": "#CC79A7",
    "line_source": "#8C6D31", "line": "#0072B2", "y": "#D55E00",
    "tumor": "#B51F3A", "healthy": "#009E73",
}


def set_theme():
    """Use a consistent LaTeX-rendered figure theme."""
    plt.rcParams.update({
        "text.usetex": True,
        "font.family": "serif",
        "font.serif": ["Computer Modern Roman"],
        "font.size": 11,
        "axes.titlesize": 12,
        "axes.labelsize": 11,
        "legend.fontsize": 9,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": True,
        "axes.axisbelow": True,
        "grid.alpha": 0.18,
        "lines.linewidth": 1.9,
        "lines.markersize": 4,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.1,
        "pdf.fonttype": 42,
    })


def save(figure, directory, name):
    """Save a vector PDF."""
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)
    figure.savefig(directory / f"{name}.pdf", format="pdf")
    plt.close(figure)
    return name


def toy_geometry(directory):
    """Show the source, vessels, and tumor region."""
    figure, axes = plt.subplots(1, 2, figsize=(10, 3.6), layout="constrained")
    for ax, name in zip(axes, ("line", "y")):
        toy = make_toy(name)
        xy = np.asarray(toy.points)[:, :2] * 1000
        center = np.array(toy.tumor_center)[:2] * 1000
        ax.add_patch(Circle(center, toy.tumor_radius * 1000, color=COLORS["tumor"], alpha=0.16))
        for branch, (a, b) in enumerate(toy.branches):
            ax.plot(xy[[a, b], 0], xy[[a, b], 1], color="#333333", lw=3)
        ax.annotate("Blood enters", xy=xy[0], xytext=(0, 5), textcoords="offset points",
                    ha="left", va="bottom")
        ax.annotate("Tumor", xy=center, xytext=(-12, 25), textcoords="offset points",
                    color=COLORS["tumor"])
        if name == "y":
            ax.text(19, -9.7, "Healthy outlet", ha="center", fontsize=10)
        ax.scatter(*xy[0], color=COLORS["lateral"], s=35, zorder=4)
        ax.set(xlim=(0, 32), ylim=(-12, 13), aspect="equal", xlabel="Distance (mm)",
               ylabel="Distance (mm)", title="Straight vessel" if name == "line" else "Y-shaped vessel")
    return save(figure, directory, "toy_geometry")


def intact_flows(rows, directory):
    """Compare absolute tumor flow for the intact vessels."""
    data = [r for r in rows if r["experiment"] == "intact"]
    figure, ax = plt.subplots(figsize=(11, 4), layout="constrained")
    x = np.arange(len(MODELS))
    for offset, name in [(-0.19, "line"), (0.19, "y")]:
        values = [next(r["tumor_flow"] for r in data if r["toy"] == name and r["model"] == m)
                  for m in MODELS]
        ax.bar(x + offset, values, 0.36, color=COLORS[name],
               label="Straight vessel" if name == "line" else "Y-shaped vessel")
    ax.set_xticks(x, [MODEL_LABELS[m].replace(" ", "\n", 1) for m in MODELS])
    ax.set(ylabel=r"Tumor flow ($\mathrm{mL}\,\mathrm{min}^{-1}$)", title="Tumor flow in intact vessels")
    ax.legend(frameon=False)
    return save(figure, directory, "intact_tumor_flow")


def model_panels(rows, directory, experiment, toy, xkey, xlabel, name, ylabel="Tumor flow remaining (percent)"):
    """Show one simple panel for each model."""
    figure, axes = plt.subplots(2, 4, figsize=(13.3, 6.4), layout="constrained", sharey=True)
    subset = [r for r in rows if r["experiment"] == experiment and r["toy"] == toy]
    branch_labels = {0: "Straight vessel"} if toy == "line" else {
        0: "Parent vessel", 1: "Tumor branch", 2: "Healthy branch"}
    branch_colors = {0: "#333333", 1: COLORS["tumor"], 2: COLORS["healthy"]}
    for ax, model in zip(axes.ravel(), MODELS):
        for branch, label in branch_labels.items():
            data = sorted([r for r in subset if r["model"] == model and r["branch"] == branch],
                          key=lambda r: r[xkey])
            if data:
                ax.plot([r[xkey] for r in data], [r["remaining_percent"] for r in data],
                        "o-", color=branch_colors[branch], label=label)
        ax.axhline(100, color="#888888", ls=":", lw=1)
        ax.set(title=MODEL_LABELS[model], xlabel=xlabel, ylabel=ylabel)
        ax.set_ylim(bottom=min(0, ax.get_ylim()[0]))
    legend_ax = axes.ravel()[-1]
    legend_ax.axis("off")
    handles, labels = axes.ravel()[0].get_legend_handles_labels()
    legend_ax.legend(handles, labels, loc="center", frameon=False)
    figure.suptitle(("Straight vessel" if toy == "line" else "Y-shaped vessel") + ": " + {
        "gap_position": "where the gap occurs", "gap_length": "how long the gap is",
        "radius": "local radius errors", "repair": "radius of the restored vessel",
    }.get(experiment, experiment.replace("_", " ")))
    return save(figure, directory, name)


def resolution(rows, directory, experiment, name):
    """Compare tissue or vessel resolution on the Y geometry."""
    figure, axes = plt.subplots(1, 2, figsize=(10, 4), layout="constrained")
    selection = [r for r in rows if r["experiment"] == experiment and r["toy"] == "y"]
    models = ("lateral", "terminal", "contact", "spread", "line_source")
    spacing = "tissue_spacing_mm" if experiment == "tissue_mesh" else "vessel_spacing_mm"
    for model in models:
        intact = sorted([r for r in selection if r["model"] == model and r["defect"] == "none"],
                        key=lambda r: r[spacing])
        broken = sorted([r for r in selection if r["model"] == model and r["defect"] == "gap"],
                        key=lambda r: r[spacing])
        axes[0].plot([r[spacing] for r in intact], [r["tumor_flow"] for r in intact],
                     "o-", color=COLORS[model], label=MODEL_LABELS[model])
        axes[1].plot([r[spacing] for r in broken], [r["remaining_percent"] for r in broken],
                     "o-", color=COLORS[model], label=MODEL_LABELS[model])
    for ax in axes:
        ax.set_xlabel("Tissue spacing (mm)" if experiment == "tissue_mesh" else "Largest vessel element (mm)")
        ax.invert_xaxis()
    axes[0].set(ylabel=r"Tumor flow ($\mathrm{mL}\,\mathrm{min}^{-1}$)", title="Intact Y-shaped vessel")
    axes[1].set(ylabel="Tumor flow remaining (percent)", title="One millimeter gap in the tumor branch")
    axes[0].legend(frameon=False)
    return save(figure, directory, name)


def parameter_plot(rows, directory, experiment, xkey, xlabel, name, models):
    """Show one parameter change on each toy geometry."""
    figure, axes = plt.subplots(1, 2, figsize=(10, 4), layout="constrained", sharey=True)
    for ax, toy in zip(axes, ("line", "y")):
        for model in models:
            data = sorted([r for r in rows if r["experiment"] == experiment and r["toy"] == toy
                           and r["model"] == model and r["defect"] == "gap"], key=lambda r: r[xkey])
            ax.plot([r[xkey] for r in data], [r["remaining_percent"] for r in data], "o-",
                    color=COLORS[model], label=MODEL_LABELS[model])
        ax.set(xlabel=xlabel, title="Straight vessel" if toy == "line" else "Y-shaped vessel")
    axes[0].set_ylabel("Tumor flow remaining (percent)")
    axes[0].legend(frameon=False)
    return save(figure, directory, name)


def extra_connection(rows, directory):
    """Show redistribution caused by a false connection."""
    figure, ax = plt.subplots(figsize=(7, 4), layout="constrained")
    for model in ("poiseuille", "lateral", "terminal", "contact", "spread", "line_source"):
        data = sorted([r for r in rows if r["experiment"] == "extra" and r["model"] == model],
                      key=lambda r: r["radius_factor"])
        ax.plot([r["radius_factor"] * 100 for r in data], [r["remaining_percent"] - 100 for r in data],
                "o-", color=COLORS[model], label=MODEL_LABELS[model])
    ax.axhline(0, color="#888888", ls=":", lw=1)
    ax.set(xlabel="False connection radius (percent of daughter radius)",
           ylabel=r"Change in tumor flow (percent)", title="An extra connection between the daughter vessels")
    ax.legend(frameon=False)
    return save(figure, directory, "extra_connection")


def pressure_profiles(solutions, directory):
    """Show the pressure lost across a missing and restored vessel."""
    figure, axes = plt.subplots(1, 2, figsize=(10, 4), layout="constrained", sharey=True)
    styles = {"Intact": ("#333333", "-"), "Gap": (COLORS["tumor"], "--"),
              "Restored": (COLORS["healthy"], ":")}
    for ax, model in zip(axes, ("poiseuille", "lateral")):
        for label, (color, style) in styles.items():
            solution = solutions[(model, label)]
            for dofs in solution.vessel.cell_dofs:
                x = (solution.vessel.coordinates[dofs, 0] - 0.003) * 1000
                y = solution.vessel_pressure[dofs]
                ax.plot(x, y, style, color=color)
            ax.plot([], [], style, color=color, label=label)
        ax.set(title=MODEL_LABELS[model], xlabel="Distance from inlet end (mm)")
        ax.legend(frameon=False)
    axes[0].set_ylabel("Vessel pressure (mmHg)")
    return save(figure, directory, "pressure_along_straight_vessel")


def equations(directory):
    """Summarize the compared mechanisms in familiar equations."""
    figure, ax = plt.subplots(figsize=(8, 4), layout="constrained")
    ax.axis("off")
    for i, model in enumerate(MODELS):
        ax.text(0.02, 0.93 - 0.135 * i, MODEL_LABELS[model], va="center")
        ax.text(0.43, 0.93 - 0.135 * i, MODEL_EQUATIONS[model], va="center", fontsize=14)
    return save(figure, directory, "model_equations")


def conservative_flux(rows, directory):
    """Compare cell flux and continuous pressure formulations."""
    figure, axes = plt.subplots(1, 2, figsize=(10, 4), layout="constrained")
    for ax, model in zip(axes, ("terminal", "spread")):
        for experiment, key, label, color in (
                ("tissue_mesh", model, "Continuous pressure", COLORS[model]),
                ("conservative_flux", "mixed_" + model, "Shared face flux", "#333333")):
            data = sorted([r for r in rows if r["experiment"] == experiment and r["toy"] == "y"
                           and r["model"] == key and r["defect"] == "none"],
                          key=lambda r: r["tissue_spacing_mm"])
            ax.plot([r["tissue_spacing_mm"] for r in data], [r["tumor_flow"] for r in data],
                    "o-", label=label, color=color)
        ax.set(xlabel="Tissue spacing (mm)", title=MODEL_LABELS[model])
        ax.invert_xaxis()
        ax.legend(frameon=False)
    axes[0].set_ylabel(r"Tumor flow ($\mathrm{mL}\,\mathrm{min}^{-1}$)")
    return save(figure, directory, "conservative_flux")


def gap_effect_map(rows, directory):
    """Color vessel positions by their effect on tumor delivery."""
    figure, axes = plt.subplots(1, 3, figsize=(11, 4), layout="constrained")
    toy = make_toy("y")
    xy = np.array(toy.points)[:, :2] * 1000
    image = None
    for ax, model in zip(axes, ("poiseuille", "lateral", "contact")):
        ax.add_patch(Circle(np.array(toy.tumor_center)[:2] * 1000, toy.tumor_radius * 1000,
                            color=COLORS["tumor"], alpha=0.12))
        for a, b in toy.branches:
            ax.plot(xy[[a, b], 0], xy[[a, b], 1], color="#cccccc", lw=4, zorder=1)
        for branch, (a, b) in enumerate(toy.branches):
            data = sorted([r for r in rows if r["experiment"] == "gap_position" and r["toy"] == "y"
                           and r["model"] == model and r["branch"] == branch], key=lambda r: r["position_percent"])
            t = np.array([r["position_percent"] / 100 for r in data])
            positions = xy[a] + t[:, None] * (xy[b] - xy[a])
            image = ax.scatter(positions[:, 0], positions[:, 1],
                               c=[r["remaining_percent"] for r in data], vmin=0, vmax=125,
                               cmap="viridis", s=40, zorder=3)
        ax.set(aspect="equal", title=MODEL_LABELS[model], xlim=(0, 32), ylim=(-12, 13),
               xlabel="Distance (mm)", ylabel="Distance (mm)")
    figure.colorbar(image, ax=axes, label="Tumor flow remaining after a gap (percent)", shrink=0.75)
    return save(figure, directory, "gap_effect_map")


def bifurcation_gap(rows, directory):
    """Compare flow after removing the center of the bifurcation."""
    figure, ax = plt.subplots(figsize=(10, 4), layout="constrained")
    data = [next(r for r in rows if r["experiment"] == "junction_gap" and r["model"] == model
                 and abs(r["gap_mm"] - 1) < 1e-8) for model in MODELS]
    ax.bar(np.arange(len(MODELS)), [r["remaining_percent"] for r in data],
           color=[COLORS[m] for m in MODELS])
    ax.set_xticks(np.arange(len(MODELS)), [MODEL_LABELS[m].replace(" ", "\n", 1) for m in MODELS])
    ax.set(ylabel="Tumor flow remaining (percent)", title="One millimeter gap centered on the bifurcation")
    ax.axhline(100, color="#888888", ls=":", lw=1)
    return save(figure, directory, "bifurcation_gap")


def make_figures(rows, directory, solutions):
    """Create all experiment PDFs with the shared theme."""
    set_theme()
    names = [toy_geometry(directory), intact_flows(rows, directory), equations(directory)]
    for toy in ("line", "y"):
        names.append(model_panels(rows, directory, "gap_position", toy, "position_percent",
                                  "Gap position along vessel (percent)", f"gap_position_{toy}"))
        names.append(model_panels(rows, directory, "gap_length", toy, "gap_mm",
                                  "Gap length (mm)", f"gap_length_{toy}"))
        names.append(model_panels(rows, directory, "radius", toy, "radius_percent",
                                  "Change in local radius (percent)", f"radius_change_{toy}"))
        names.append(model_panels(rows, directory, "repair", toy, "radius_factor",
                                  "Bridge radius relative to intact radius", f"bridge_radius_{toy}"))
    names += [resolution(rows, directory, "tissue_mesh", "tissue_resolution"),
              resolution(rows, directory, "vessel_mesh", "vessel_resolution"),
              parameter_plot(rows, directory, "terminal_support", "terminal_radius_mm",
                             "Terminal source radius (mm)", "terminal_source_radius",
                             ("lateral", "terminal", "contact", "spread", "line_source")),
              parameter_plot(rows, directory, "wall_strength", "wall_factor",
                             "Wall conductivity relative to default", "wall_conductivity",
                             ("lateral", "spread", "line_source")),
              parameter_plot(rows, directory, "spread_width", "spread_factor",
                             "Spread radius relative to vessel radius", "spread_radius", ("spread",)),
              extra_connection(rows, directory), pressure_profiles(solutions, directory),
              conservative_flux(rows, directory), gap_effect_map(rows, directory), bifurcation_gap(rows, directory)]
    return names
