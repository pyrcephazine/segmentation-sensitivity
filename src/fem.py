"""Finite element meshes and interpolation operators."""

from functools import lru_cache

import basix.ufl
import networkx as nx
import numpy as np
import ufl
from dolfinx import fem, geometry, mesh
from dolfinx.fem import petsc as fem_petsc
from mpi4py import MPI
from networks_fenicsx import NetworkMesh
from scipy import sparse
from scipy.stats import qmc

from geometry import make_graph


def assemble(form):
    """Assemble a scalar finite element matrix."""
    matrix = fem_petsc.assemble_matrix(fem.form(form))
    matrix.assemble()
    indptr, indices, data = matrix.getValuesCSR()
    result = sparse.csr_matrix((data.copy(), indices.copy(), indptr.copy()), shape=matrix.getSize())
    matrix.destroy()
    return result


class Bulk:
    """A tetrahedral pressure mesh."""

    def __init__(self, toy, resolution):
        counts = [resolution, max(4, round(0.8 * resolution)), max(3, round(0.4 * resolution))]
        self.mesh = mesh.create_box(MPI.COMM_SELF, [np.array(toy.lower), np.array(toy.upper)],
                                    counts, cell_type=mesh.CellType.tetrahedron)
        self.space = fem.functionspace(self.mesh, ("Lagrange", 1))
        p, v = ufl.TrialFunction(self.space), ufl.TestFunction(self.space)
        self.stiffness = assemble(ufl.inner(ufl.grad(p), ufl.grad(v)) * ufl.dx)
        self.mass = assemble(p * v * ufl.dx)
        self.size = self.space.dofmap.index_map.size_local
        self.coordinates = self.space.tabulate_dof_coordinates()
        self.spacing = float(np.max((np.array(toy.upper) - toy.lower) / counts))
        self.tree = geometry.bb_tree(self.mesh, 3)
        nc = self.mesh.topology.index_map(3).size_local
        self.cell_dofs = np.array([self.space.dofmap.cell_dofs(i) for i in range(nc)])
        vertices = self.coordinates[self.cell_dofs]
        self.origins = vertices[:, 0]
        self.inverse = np.linalg.inv(np.transpose(vertices[:, 1:] - vertices[:, :1], (0, 2, 1)))

    def evaluation(self, points):
        """Interpolate bulk pressure at points inside the box."""
        points = np.asarray(points, dtype=np.float64).reshape(-1, 3)
        candidates = geometry.compute_collisions_points(self.tree, points)
        collisions = geometry.compute_colliding_cells(self.mesh, candidates, points)
        cells = np.array([collisions.links(i)[0] for i in range(len(points))], dtype=np.int32)
        reference = np.einsum("nij,nj->ni", self.inverse[cells], points - self.origins[cells])
        weights = np.column_stack((1 - reference.sum(axis=1), reference))
        return sparse.csr_matrix((weights.ravel(),
                                  (np.repeat(np.arange(len(points)), 4), self.cell_dofs[cells].ravel())),
                                 shape=(len(points), self.size))


@lru_cache(maxsize=8)
def bulk_mesh(toy, resolution):
    """Reuse a bulk mesh between parameter sweeps."""
    return Bulk(toy, resolution)


class Vessel:
    """A pressure space on a networks_fenicsx mesh."""

    def __init__(self, toy, resolution, defect):
        self.graph = make_graph(toy, defect)
        self.network = NetworkMesh(self.graph, N=resolution, comm=MPI.COMM_SELF)
        self.mesh = self.network.mesh
        self.space = fem.functionspace(self.mesh, ("Lagrange", 1))
        self.coordinates = self.space.tabulate_dof_coordinates()
        self.size = self.space.dofmap.index_map.size_local
        nc = self.mesh.topology.index_map(1).size_local
        self.cell_dofs = np.array([self.space.dofmap.cell_dofs(i) for i in range(nc)])
        self.endpoints = self.coordinates[self.cell_dofs]
        self.lengths = np.linalg.norm(self.endpoints[:, 1] - self.endpoints[:, 0], axis=1)
        self.tangents = (self.endpoints[:, 1] - self.endpoints[:, 0]) / self.lengths[:, None]
        markers = np.empty(nc, dtype=int)
        markers[self.network.subdomains.indices] = self.network.subdomains.values
        edges = list(self.graph.edges(data=True))
        self.radii = np.array([edges[i][2]["radius"] for i in markers])
        self.nominal_radii = np.array([edges[i][2]["nominal_radius"] for i in markers])
        self.branches = np.array([edges[i][2]["branch"] for i in markers])
        component = {node: j for j, nodes in enumerate(nx.connected_components(self.graph.to_undirected()))
                     for node in nodes}
        self.components = np.array([component[edges[i][0]] for i in markers])
        self.coefficient = fem.Function(fem.functionspace(self.mesh, ("DG", 0)))
        self.matrices = {}
        self.root = int(np.argmin(np.linalg.norm(self.coordinates - toy.points[0], axis=1)))
        terminal_nodes = [b for a, b in toy.branches if all(b != u for u, v in toy.branches)]
        self.terminals = np.array([np.argmin(np.linalg.norm(self.coordinates - toy.points[i], axis=1))
                                   for i in terminal_nodes], dtype=int)
        self.terminal_positions = self.coordinates[self.terminals]
        self.tumor_terminals = toy.in_tumor(self.terminal_positions).astype(float)
        abscissae = (np.array([-1, 1]) / np.sqrt(3) + 1) / 2
        self.quadrature_points = ((1 - abscissae[None, :, None]) * self.endpoints[:, :1]
                                  + abscissae[None, :, None] * self.endpoints[:, 1:]).reshape(-1, 3)
        self.weights = np.repeat(self.lengths / 2, 2)
        values = np.tile(np.column_stack((1 - abscissae, abscissae)), (nc, 1))
        self.N = sparse.csr_matrix((values.ravel(),
                                    (np.repeat(np.arange(2 * nc), 2), np.repeat(self.cell_dofs, 2, axis=0).ravel())),
                                   shape=(2 * nc, self.size))

    def stiffness(self, parameters, model):
        """Assemble axial Darcy or Poiseuille flow."""
        key = (model == "darcy", parameters.viscosity, parameters.conductance_scale)
        if key in self.matrices:
            return self.matrices[key]
        radii = self.radii
        if model == "darcy":
            c = np.pi * radii**2 * self.nominal_radii**2 / (8 * parameters.viscosity)
        else:
            c = np.pi * radii**4 / (8 * parameters.viscosity)
        for i, value in enumerate(c / parameters.conductance_scale):
            self.coefficient.x.array[self.coefficient.function_space.dofmap.cell_dofs(i)] = value
        p, v = ufl.TrialFunction(self.space), ufl.TestFunction(self.space)
        result = assemble(self.coefficient * ufl.inner(ufl.grad(p), ufl.grad(v)) * ufl.dx)
        self.matrices[key] = result
        return result


@lru_cache(maxsize=24)
def vessel_mesh(toy, resolution, defect):
    """Reuse vessel meshes between model comparisons."""
    return Vessel(toy, resolution, defect)


def normal_vectors(tangents):
    """Construct two perpendicular directions for each vessel cell."""
    guide = np.tile([0.0, 0.0, 1.0], (len(tangents), 1))
    guide[np.abs(tangents[:, 2]) > 0.9] = [1.0, 0.0, 0.0]
    first = np.cross(tangents, guide)
    first /= np.linalg.norm(first, axis=1)[:, None]
    return first, np.cross(tangents, first)


@lru_cache(maxsize=48)
def lateral_operator(bulk, vessel, toy, parameters, kind):
    """Average pressure and tumor membership around each vessel point."""
    points = vessel.quadrature_points
    if kind == "line":
        return bulk.evaluation(points), toy.in_tumor(points).astype(float)
    tangents = np.repeat(vessel.tangents, 2, axis=0)
    first, second = normal_vectors(tangents)
    angles = np.arange(parameters.angles) * 2 * np.pi / parameters.angles
    direction = (np.cos(angles)[None, :, None] * first[:, None]
                 + np.sin(angles)[None, :, None] * second[:, None])
    radii = np.repeat(vessel.radii, 2)
    if kind == "disk":
        radii = radii * parameters.spread_radius_factor
        nodes, weights = np.polynomial.legendre.leggauss(parameters.disk_rings)
        factors = np.sqrt((nodes + 1) / 2)
        weights = weights / 2
    else:
        factors, weights = np.array([1.0]), np.array([1.0])
    samples = (points[:, None, None, :] + radii[:, None, None, None]
               * factors[None, :, None, None] * direction[:, None])
    sample_weights = np.repeat(weights / parameters.angles, parameters.angles)
    n = len(sample_weights)
    average = sparse.csr_matrix((np.tile(sample_weights, len(points)),
                                (np.repeat(np.arange(len(points)), n), np.arange(len(points) * n))),
                               shape=(len(points), len(points) * n))
    values = bulk.evaluation(samples.reshape(-1, 3))
    tumor = toy.in_tumor(samples.reshape(-1, 3)).astype(float)
    return (average @ values).tocsr(), np.asarray(average @ tumor).ravel()


@lru_cache(maxsize=32)
def terminal_operator(bulk, vessel, toy, parameters):
    """Average over normalized terminal spheres."""
    sample = qmc.Sobol(3, scramble=False).random_base2(int(np.log2(parameters.terminal_samples)))
    radius = parameters.terminal_radius * sample[:, 0] ** (1 / 3)
    cosine = 2 * sample[:, 1] - 1
    angle = 2 * np.pi * sample[:, 2]
    direction = np.column_stack((np.sqrt(1 - cosine**2) * np.cos(angle),
                                 np.sqrt(1 - cosine**2) * np.sin(angle), cosine))
    offsets = radius[:, None] * direction
    rows, fractions = [], []
    lower, upper = np.array(toy.lower), np.array(toy.upper)
    for point in vessel.terminal_positions:
        points = point + offsets
        points = points[np.all((points >= lower) & (points <= upper), axis=1)]
        values = bulk.evaluation(points)
        rows.append(sparse.csr_matrix(np.asarray(values.mean(axis=0))))
        fractions.append(np.mean(toy.in_tumor(points)))
    return sparse.vstack(rows, format="csr"), np.array(fractions)


def contact_average(vessel, toy, spacing):
    """Use fixed physical wall intervals for pressure constraints."""
    rows, columns, weights, lookup = [], [], [], {}
    branches = np.repeat(vessel.branches, 2)
    components = np.repeat(vessel.components, 2)
    for i, (point, branch, component) in enumerate(zip(vessel.quadrature_points, branches, components)):
        if branch < len(toy.branches):
            origin = np.array(toy.points[toy.branches[branch][0]])
        else:
            edge = next((a, b) for a, b, data in vessel.graph.edges(data=True) if data["branch"] == branch)
            origin = np.array(vessel.graph.nodes[edge[0]]["pos"])
        interval = int(np.floor((np.linalg.norm(point - origin) + 1e-12) / spacing))
        key = (int(branch), int(component), interval)
        if key not in lookup:
            lookup[key] = len(lookup)
        rows.append(lookup[key])
        columns.append(i)
        weights.append(vessel.weights[i])
    matrix = sparse.csr_matrix((weights, (rows, columns)), shape=(len(lookup), len(branches)))
    return sparse.diags(1 / np.asarray(matrix.sum(axis=1)).ravel()) @ matrix


class MixedBulk:
    """Shared face fluxes and cell pressures for volume exchange."""

    def __init__(self, toy, resolution, parameters):
        base = bulk_mesh(toy, resolution)
        self.mesh, self.tree, self.spacing = base.mesh, base.tree, base.spacing
        self.space = fem.functionspace(self.mesh, ("DG", 0))
        self.size = self.space.dofmap.index_map.size_local
        self.cell_dofs = np.array([self.space.dofmap.cell_dofs(i)[0]
                                  for i in range(self.mesh.topology.index_map(3).size_local)])
        self.coordinates = self.space.tabulate_dof_coordinates()
        flux_space = fem.functionspace(self.mesh, basix.ufl.element("RT", self.mesh.basix_cell(), 1))
        self.flux_size = flux_space.dofmap.index_map.size_local
        self.mesh.topology.create_connectivity(2, 3)
        facets = mesh.exterior_facet_indices(self.mesh.topology)
        self.boundary_flux_dofs = fem.locate_dofs_topological(flux_space, 2, facets)
        u, w = ufl.TrialFunction(flux_space), ufl.TestFunction(flux_space)
        p, v = ufl.TrialFunction(self.space), ufl.TestFunction(self.space)
        self.H = assemble((parameters.conductance_scale / parameters.tissue_mobility)
                          * ufl.inner(u, w) * ufl.dx)
        self.B = assemble(ufl.div(u) * v * ufl.dx)
        self.mass = assemble(p * v * ufl.dx)

    def evaluation(self, points):
        """Read cell pressure at volume quadrature points."""
        points = np.asarray(points, dtype=np.float64).reshape(-1, 3)
        candidates = geometry.compute_collisions_points(self.tree, points)
        collisions = geometry.compute_colliding_cells(self.mesh, candidates, points)
        cells = np.array([collisions.links(i)[0] for i in range(len(points))])
        return sparse.csr_matrix((np.ones(len(points)), (np.arange(len(points)), self.cell_dofs[cells])),
                                 shape=(len(points), self.size))


@lru_cache(maxsize=4)
def mixed_bulk_mesh(toy, resolution, parameters):
    """Reuse mixed volume meshes."""
    return MixedBulk(toy, resolution, parameters)
