"""Shared sparse assembly for all pressure and exchange models."""

from dataclasses import dataclass

import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve

from config import Parameters
from fem import bulk_mesh, contact_average, lateral_operator, mixed_bulk_mesh, terminal_operator, vessel_mesh
from geometry import Defect


@dataclass
class Solution:
    """Flow outputs and optional spatial fields."""

    tumor_flow: float
    terminal_tumor_flow: float
    lateral_tumor_flow: float
    inlet_flow: float
    terminal_flows: np.ndarray
    vessel_pressure: np.ndarray
    tissue_pressure: np.ndarray
    vessel: object
    bulk: object
    dofs: int
    vessel_spacing: float
    tissue_spacing: float


def solve_system(matrix, fixed, values):
    """Solve with prescribed pressure or zero boundary flux values."""
    fixed = np.asarray(fixed, dtype=int)
    free = np.ones(matrix.shape[0], dtype=bool)
    free[fixed] = False
    result = np.zeros(matrix.shape[0])
    result[fixed] = values
    rhs = -matrix[free][:, fixed] @ result[fixed]
    result[free] = spsolve(matrix[free][:, free].tocsc(), rhs)
    return result


def solve(toy, model="lateral", parameters=Parameters(), defect=Defect()):
    """Solve one toy geometry and return flow in milliliters per minute."""
    if model == "bridge":
        model = "poiseuille"
        if defect.kind == "gap":
            defect = Defect("bridge", defect.branch, defect.position, defect.length, defect.radius_factor)
    vessel = vessel_mesh(toy, parameters.vessel_resolution, defect)
    Kv = vessel.stiffness(parameters, model)
    nv = vessel.size
    E = sparse.csr_matrix((np.ones(len(vessel.terminals)),
                           (np.arange(len(vessel.terminals)), vessel.terminals)),
                          shape=(len(vessel.terminals), nv))
    conductances = (parameters.terminal_conductance / parameters.conductance_scale
                    * np.where(vessel.tumor_terminals > 0, 1.0, parameters.healthy_terminal_factor))
    G = sparse.diags(conductances, format="csr")
    bulk, pressure, tissue, multipliers = None, None, np.empty(0), None
    lateral_flow = 0.0
    tumor_weights = vessel.tumor_terminals
    jump, weights, tumor_lateral = None, None, None
    if model in {"poiseuille", "darcy"}:
        matrix = Kv + E.T @ G @ E
        solution = solve_system(matrix, [vessel.root], [1.0])
        pressure = solution
        terminal_flow = G @ E @ pressure
        root = vessel.root
        tissue_spacing = np.nan
    else:
        mixed = model.startswith("mixed_")
        physical_model = model.removeprefix("mixed_")
        bulk = (mixed_bulk_mesh(toy, parameters.tissue_resolution, parameters) if mixed
                else bulk_mesh(toy, parameters.tissue_resolution))
        nt = bulk.size
        At, tumor_weights = terminal_operator(bulk, vessel, toy, parameters)
        D = sparse.hstack((-At, E), format="csr")
        Kt = (parameters.drainage / parameters.conductance_scale) * bulk.mass
        if not mixed:
            Kt = Kt + (parameters.tissue_mobility / parameters.conductance_scale) * bulk.stiffness
        base = sparse.block_diag((Kt, Kv), format="csr") + D.T @ G @ D
        if physical_model != "terminal":
            kind = "disk" if physical_model == "spread" else "circle"
            T, tumor_lateral = lateral_operator(bulk, vessel, toy, parameters, kind)
            jump = sparse.hstack((-T, vessel.N), format="csr")
            if physical_model == "contact":
                cell_average = contact_average(vessel, toy, parameters.multiplier_spacing)
                B = cell_average @ jump
                matrix = sparse.bmat([[base, B.T], [B, None]], format="csr")
            else:
                beta = 2 * np.pi * np.repeat(vessel.radii, 2) * parameters.wall_conductivity
                if physical_model == "spread":
                    correction = ((np.log(parameters.spread_radius_factor) - 0.25)
                                  / (2 * np.pi * parameters.tissue_mobility))
                    beta = 1 / (1 / beta + correction)
                weights = vessel.weights * beta / parameters.conductance_scale
                if physical_model == "line_source":
                    source, tumor_lateral = lateral_operator(bulk, vessel, toy, parameters, "line")
                    deposition = sparse.hstack((-source, vessel.N), format="csr")
                else:
                    deposition = jump
                matrix = base + deposition.T @ sparse.diags(weights) @ jump
        else:
            matrix = base
        if mixed:
            nf = bulk.flux_size
            coupling = sparse.vstack((bulk.B, sparse.csr_matrix((nv, nf))), format="csr")
            matrix = sparse.bmat([[-bulk.H, coupling.T], [coupling, matrix]], format="csr")
            root = nf + nt + vessel.root
            fixed = np.r_[bulk.boundary_flux_dofs, root]
            solution = solve_system(matrix, fixed, np.r_[np.zeros(len(fixed) - 1), 1.0])
            tissue = solution[nf:nf + nt]
            pressure = solution[nf + nt:]
            primary = solution[nf:]
        else:
            root = nt + vessel.root
            solution = solve_system(matrix, [root], [1.0])
            primary = solution[:nt + nv]
            tissue, pressure = primary[:nt], primary[nt:]
        terminal_flow = G @ D @ primary
        if physical_model == "contact":
            multipliers = solution[nt + nv:]
            lateral_flow = np.dot(np.asarray(cell_average @ tumor_lateral).ravel(), multipliers)
        elif jump is not None:
            lateral_flow = np.dot(tumor_lateral * weights, jump @ primary)
        tissue_spacing = bulk.spacing
    flow_factor = parameters.conductance_scale * (parameters.inlet_pressure - parameters.venous_pressure) * 6e7
    terminal_tumor = float(np.dot(tumor_weights, terminal_flow))
    inlet = float((matrix @ solution)[root])
    pressure_factor = (parameters.inlet_pressure - parameters.venous_pressure) / 133.322
    pressure_offset = parameters.venous_pressure / 133.322
    return Solution((terminal_tumor + lateral_flow) * flow_factor,
                    terminal_tumor * flow_factor, float(lateral_flow) * flow_factor,
                    inlet * flow_factor, np.asarray(terminal_flow) * flow_factor,
                    pressure_offset + pressure_factor * pressure,
                    pressure_offset + pressure_factor * tissue,
                    vessel, bulk, matrix.shape[0], float(vessel.lengths.max()), tissue_spacing)
