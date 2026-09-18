# Supplied-network simulation: pt002 inflow

Completed supplement, 11 September 2026. This run solves the intact inflow network from `oilblooddata.zip`, coupled to its supplied tissue mesh. A second solve checks the effect of halving the vessel cell size. Results are dimensionless.

## Completed mesh files

* [coupled_solution.vtu](coupled_solution.vtu): self-contained tissue and vessel mesh with both computed potentials and flux fields. This is the main deliverable and needs no companion file.
* [network_solution.vtu](network_solution.vtu): vessel-only solution, convenient for a Tube filter in ParaView.
* [tissue_solution.vtu](tissue_solution.vtu): full tissue solution.
* [metrics.csv](metrics.csv): the main and refinement measurements.
* [verification.json](verification.json): numerical and artifact checks.

The VTU coordinates are in the original millimetre coordinate frame. `potential_dimensionless` contains tissue u and vessel U. Cell field `compartment_label` is 1 for tissue and 2 for the inflow network; select label 2 with Threshold to expose the embedded network. The axial and tissue flux names end in `dimensionless`; wall exchange is per unit normalized arclength and is positive from vessel to tissue. Signed axial flow follows the source VTP edge orientation. Fields belonging to the other compartment contain zero placeholders. Source VTK point and edge IDs are zero-based; `source_exodus_node_id` preserves the one-based Exodus node IDs. Newly inserted vessel points have source point ID -1. Native XDMF/HDF5 files in the raw run directory retain normalized coordinates.

## Data and inlet selection

The input files are `oilblooddata/thermoembo_run/pt002/pt002_inflow_centerline.vtp` and `oilblooddata/thermoembo_run/pt002/mesh.exo`. The archive README identifies the centerline coordinates as millimetres and the Exodus coordinates as metres. They were converted to one coordinate system before coupling.

The original network has 797 vertices, 828 edges, 55 endpoints, one connected component, and cycle rank 32. All edges and loops were retained, including 83 edges marked `is_phantom` in the source. Edge subdivision preserves their geometry and connectivity. All stored point and edge radii are 0.1 mm; this supplied value was retained rather than treated as a measured radius distribution.

The inlet is original point **0**, at **[66.1982421875, 99.6728515625, 475.5] mm**. It is the unique degree-one vertex with the highest archived endpoint pressure, 836 mmHg. This supplies an explicit computational inlet selection. Archived pressures and flows come from a different resistance model and are used only to identify this endpoint; they are not calibration or validation targets for the new PDE solution.

## Model and normalization

The model is the same Laurino-Zunino circular-average 3D-1D weak formulation used in the synthetic pilot:

```text
Integral_Omega K3 grad(u).grad(v)
 + Integral_Lambda K1 A U_prime V_prime
 + Integral_Lambda kappa P (U - T_r u)(V - T_r v) = 0.
```

K3=K1=1, kappa=0.05, A=pi*r^2, and P=2*pi*r. U=1 at the selected inlet; other vessel endpoints have natural zero axial flux. The entire exterior tissue boundary has u=0. Tissue and vessel potentials are continuous P1; circular averages use DG1 with quadrature degree 20. Source terms are zero. The circle checks include every endpoint interpolation location and additional midpoint samples.

Lengths are normalized by the longest extent of the supplied tissue mesh, L=202.380554199 mm, with origin [-0.11389311981201172, 0.041460384368896486, 0.408371337890625] metres. Thus r=0.0004941186192 in the PDE. The imported tissue mesh is preserved: 126,569 tetrahedra and 27,947 vertices. Its original volume is 887.609130 mL. The main vessel discretization has 2,016 cells and 1,985 vertices, with maximum normalized cell length 0.005 (about 1.0119 mm). The coupled system has 29,932 pressure degrees of freedom and was solved on 4 MPI ranks using PETSc/MUMPS.

These are model potentials and dimensionless fluxes. The coefficients and inlet value have not been fitted to physiological measurements. The explicit use of the supplied radius and organ domain also means the numerical percentages from the earlier synthetic examples should not be transferred to this network.

## Results

| Quantity | Computed value |
| --- | --- |
| Conservative inlet flux | 1.23142857635e-05 |
| Integrated vessel-to-tissue exchange | 1.23142857635e-05 |
| Tissue exterior outflow | 1.23142857635e-05 |
| Mean tissue potential (volume-normalized) | 3.08670848106e-07 |
| Maximum tissue nodal potential | 2.95045856171e-05 |
| Vessel nodal potential range | 1.43228924415e-11 to 1 |
| Assembly / linear solve time | 115.91 s / 0.29 s |

## Validation

| Check | Result |
| --- | --- |
| Maximum relative flux imbalance | 1.665e-14 |
| Relative linear residual | 1.234e-19 |
| Relative energy-identity error | 8.254e-16 |
| Positive inlet value | 1 |
| Main coupling containment | 66,528 samples; zero outside |
| Halved maximum 1D cell size | 3,564 vessel cells; inlet-flux change 0.009984% |
| Existing toy-model regression | 4.530e-14 maximum relative difference |
| Mesh preservation and VTU read-back | Passed |

The refinement check concerns the one-dimensional discretization on the supplied tissue mesh. It does not establish three-dimensional mesh convergence. Input checksums and coordinate transformation are in `input_metadata.json`; source IDs are in the VTU files. Both runs' metrics are in `metrics.csv`. Extracted inputs and native solutions remain in the raw run directory on NOTS.

## Reproduce

From the study source directory on NOTS:

```bash
source environment.sh
python -m experiments.pt002_inflow.run --prepare
python3 submit.py supplied
# After that job completes:
python -m experiments.pt002_inflow.report
python -m experiments.collect
```

Reference: [Laurino and Zunino (2019), equations (11)-(12)](https://doi.org/10.1051/m2an/2019042).
