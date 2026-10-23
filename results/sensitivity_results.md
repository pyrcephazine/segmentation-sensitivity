# Tumor flow sensitivity in a straight vessel and a Y-shaped vessel

Computed **1,270 cases** using DOLFINx 0.10.0 and networks_fenicsx. These are illustrative toy problems with a prescribed tumor region, not fitted animal data.

## Main findings

- At the default tissue spacing, a 1 mm gap halfway along the Y tumor branch leaves essentially zero delivery in the vessel-only and terminal-exchange models, 0.71% under wall exchange, 0.63% under spread exchange, and 77.1% under equal pressures.
- The equal-pressure prediction after that gap ranges from 47.1% to 110.1% across tissue meshes. Its large residual delivery is therefore strongly dependent on discretization, rather than evidence of robust vessel repair.
- Cutting the healthy branch instead raises terminal-model tumor flow by 22.2%, showing that off-target vessels can control the target flow.
- Restoring the missing interval with its original radius recovers the matched intact solve. The radius sweeps show how quickly this changes when the bridge is too narrow.

## Setup

Blood enters the upstream end at 80 mmHg. Venous pressure is 10 mmHg. There is no independent tissue source: vessel exchange supplies the blood bed, and distributed venous drainage removes blood. The outside of the box has zero normal flow. Only anatomical outlets exchange through terminal spheres; newly cut vessel ends are closed.

The straight vessel is 24 mm long. The Y has an 11 mm parent vessel and two 13.9 mm daughter vessels. A sphere of radius 4.5 mm around the straight vessel tip or the upper Y tip designates tumor. In the plots, the tumor region is shaded red.

$$J_T=\sum_a w_aQ_a+\int_\Lambda w_T(s)j(s)\,ds,\qquad Q_a=G_a(P_a-\overline p_a).$$

Here the weights are the fraction of each source deposited inside the tumor sphere. The 3D field is unresolved vascular-bed pressure. All percentages compare with the same model and parameters on the intact geometry, retaining the same vessel partitions around the defect so that splitting a vessel does not masquerade as damage. Signed exchange is retained.

[Toy geometry](figures/toy_geometry.pdf) · [Model equations](figures/model_equations.pdf)

| Parameter | Value |
|---|---:|
| Blood viscosity | 0.0035 Pa s |
| Straight vessel radius | 0.30 mm |
| Parent and daughter radii | 0.35 mm and 0.25 mm |
| Tissue mobility | $2\times10^{-9}$ m² Pa⁻¹ s⁻¹ |
| Venous drainage coefficient | $5\times10^{-6}$ Pa⁻¹ s⁻¹ |
| Wall conductivity | $5\times10^{-7}$ m Pa⁻¹ s⁻¹ |
| Terminal conductance | $2\times10^{-11}$ m³ Pa⁻¹ s⁻¹ |
| Terminal sphere radius | 2.5 mm |
| Default tissue spacing | 2.5 mm |

The fixed-permeability lumen model is matched to Poiseuille in the intact vessel; its permeability is held fixed when radius changes. Spread exchange uses a transverse uniform disk and the straight-cylinder pressure correction $\beta_\epsilon^{-1}=\beta^{-1}+[\log(\epsilon/r)-1/4]/(2\pi M)$. The default disk radius is five vessel radii. Equal pressures are enforced with one multiplier per fixed 4 mm wall interval, separated at a gap. The line-source model deposits on the centerline and samples pressure around the circumference; it is a resolution-dependent finite-element approximation of the singular source, not a smooth-source replacement.

## Intact flow

| Model | Straight vessel, mL min⁻¹ | Y-shaped vessel, mL min⁻¹ |
|---|---:|---:|
| Vessel flow | 7.328 | 5.909 |
| Fixed permeability | 7.328 | 5.909 |
| Wall exchange | 3.761 | 2.553 |
| Terminal exchange | 4.664 | 3.407 |
| Equal mean pressures | 2.578 | 1.529 |
| Spread exchange | 3.813 | 2.595 |
| Line source | 3.762 | 2.555 |

[Intact tumor flow](figures/intact_tumor_flow.pdf)

The direct vessel model discharges to a prescribed venous pressure; the coupled models also include resistance in the unresolved blood bed. The lower intact flow in the terminal model therefore reflects an added downstream load. Distributed exchange also changes where blood leaves the resolved vessels, so its tumor flow need not follow the change in total inlet flow.

## Gap location and length

Gaps are physical intervals removed before vessel meshing. Position is measured along the affected anatomical vessel, from its upstream end. All detached vessel pieces and original outlets remain in the solve. Location sweeps use a 1 mm gap; length sweeps use a gap centered halfway along the affected vessel.

| Model | Straight vessel, flow remaining | Y parent gap | Y tumor-branch gap | Y healthy-branch gap |
|---|---:|---:|---:|---:|
| Vessel flow | 0.0% | 0.0% | 0.0% | 107.4% |
| Fixed permeability | 0.0% | 0.0% | 0.0% | 107.4% |
| Wall exchange | 1.9% | 2.4% | 0.7% | 119.0% |
| Terminal exchange | 0.0% | 0.0% | 0.0% | 122.2% |
| Equal mean pressures | 53.0% | 60.3% | 77.1% | 107.7% |
| Spread exchange | 1.6% | 2.1% | 0.6% | 119.2% |
| Line source | 1.9% | 2.4% | 0.7% | 119.0% |

[Gap position, straight vessel](figures/gap_position_line.pdf) · [Gap position, Y-shaped vessel](figures/gap_position_y.pdf)

[Gap length, straight vessel](figures/gap_length_line.pdf) · [Gap length, Y-shaped vessel](figures/gap_length_y.pdf)

[Spatial map of gap effects](figures/gap_effect_map.pdf) · [Gap at the bifurcation](figures/bifurcation_gap.pdf)

A bifurcation gap removes the nearby half-gap length from all three incident vessels. It is a different topology change from removing an interval in only one branch.

Tissue-mediated flow after a gap is a different route from an intact arterial lumen. A nonzero tumor source therefore does not establish successful vessel reconnection. A healthy-branch defect can increase tumor delivery by changing competition for inflow.

## Radius errors and explicit repair

Radius errors affect a 3 mm interval. In the coupled models the changed radius also changes wall circumference and pressure-sampling locations. A bridge fills a 2 mm missing interval using an estimated radius. It restores the corresponding exchange terms as well as axial flow; the multiplier and line-source models retain their own exchange laws on the restored interval.

For a concrete comparison, halving radius in the Y tumor branch leaves 46.1% of intact tumor flow under Poiseuille resistance and 81.0% under fixed lumen permeability. Both laws were matched on the intact vessel, so this difference isolates radius dependence.

[Radius errors, straight vessel](figures/radius_change_line.pdf) · [Radius errors, Y-shaped vessel](figures/radius_change_y.pdf)

[Bridge radius, straight vessel](figures/bridge_radius_line.pdf) · [Bridge radius, Y-shaped vessel](figures/bridge_radius_y.pdf)

[Vessel pressure before and after repair](figures/pressure_along_straight_vessel.pdf)

## Mesh resolution

Tissue and vessel resolution are varied separately. Physical gap boundaries, radii, tumor location, source widths, and coefficients stay fixed. The charts show both absolute intact flow and the flow remaining after a gap; a stable percentage can otherwise hide a changing intact prediction.

[Tissue resolution](figures/tissue_resolution.pdf) · [Vessel resolution](figures/vessel_resolution.pdf)

| Model | Y intact flow at 5 mm tissue spacing | Y intact flow at 1.5 mm tissue spacing |
|---|---:|---:|
| Wall exchange | 2.576 | 2.552 |
| Terminal exchange | 3.471 | 3.375 |
| Equal mean pressures | 1.450 | 1.590 |
| Spread exchange | 2.627 | 2.576 |
| Line source | 2.577 | 2.560 |

The conservative volume-exchange calculations use shared face fluxes and cell pressures in the 3D bed. They are included for terminal and spread exchange rather than applying ordinary volume flux spaces to an exact surface source. The equal-pressure constraints use fixed physical wall intervals rather than one multiplier per vessel mesh element.

Despite the relatively small intact-flow change, equal-pressure flow after the gap varies from 47.1% to 110.1%. This is the clearest mesh sensitivity in the study. The finite-resistance exchange models consistently predict much less remaining tumor delivery for this gap.

[Conservative flux and pressure formulations](figures/conservative_flux.pdf)

## Exchange assumptions and false connections

Each source-width or conductivity comparison uses its own intact reference. Changing terminal source radius also changes overlap with the tumor when the source becomes wider than the tumor. The false connection joins the daughter vessels one quarter of the way beyond the bifurcation. For that experiment the healthy outlet conductance is half the tumor-outlet conductance, so a false connection can redistribute blood between branches with different loads.

Increasing wall conductivity from one tenth to ten times its default value raises remaining Y tumor flow after a gap from 0.007% to 13.7%. More exchange provides a stronger route through the blood bed, but does not restore the missing lumen. With the cylindrical pressure correction, increasing spread width from 1.5 to eight vessel radii changes remaining flow only from 0.70% to 0.62%.

For terminal exchange, increasing source radius from 4 to 6 mm reduces intact tumor-assigned flow from 3.506 to 1.670 mL min⁻¹. The wider source exceeds the 4.5 mm tumor radius and assigns more blood outside the target. A false connection as wide as a daughter vessel raises tumor flow by 1.3% in the direct vessel model but lowers it by 2.3% with wall exchange, which also adds exchange along the false vessel.

[Terminal source radius](figures/terminal_source_radius.pdf) · [Wall conductivity](figures/wall_conductivity.pdf) · [Spread radius](figures/spread_radius.pdf) · [Extra connection](figures/extra_connection.pdf)

## Data and model definitions

- [All numerical results](data/experiments.csv)
- [Parameters and package versions](data/run.json)
- [Model derivations](../notes/tumor_flow_models_review.md)
- [Implementation and reproduction](../README.md)

All plots are vector PDFs with LaTeX-rendered text and equations. Numerical flow values are in mL min⁻¹.

Network meshes use [networks_fenicsx](https://github.com/scientificcomputing/networks_fenicsx), by Jørgen Schartum Dokken and contributors. DOLFINx assembles the vessel and tissue finite elements; SciPy solves the coupled sparse systems.
