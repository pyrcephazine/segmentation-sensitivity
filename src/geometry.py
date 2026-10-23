"""Straight and bifurcating vessels with controlled defects."""

from dataclasses import dataclass

import networkx as nx
import numpy as np


@dataclass(frozen=True)
class Toy:
    """Vessel geometry and a spherical tumor."""

    name: str
    points: tuple
    branches: tuple
    radii: tuple
    tumor_center: tuple
    tumor_radius: float = 0.0045
    lower: tuple = (0.0, -0.012, -0.006)
    upper: tuple = (0.030, 0.012, 0.006)

    def in_tumor(self, points):
        """Return tumor membership for points."""
        return np.linalg.norm(np.asarray(points) - self.tumor_center, axis=-1) <= self.tumor_radius


@dataclass(frozen=True)
class Defect:
    """A local gap, radius change, bridge, or extra connection."""

    kind: str = "none"
    branch: int = 0
    position: float = 0.5
    length: float = 0.001
    radius_factor: float = 1.0


def make_toy(name):
    """Create one of the two toy geometries."""
    if name == "line":
        return Toy("line", ((0.003, 0, 0), (0.027, 0, 0)), ((0, 1),),
                   (0.00030,), (0.027, 0, 0))
    return Toy("y", ((0.003, 0, 0), (0.014, 0, 0), (0.026, 0.007, 0),
                     (0.026, -0.007, 0)), ((0, 1), (1, 2), (1, 3)),
               (0.00035, 0.00025, 0.00025), (0.026, 0.007, 0))


def make_graph(toy, defect=Defect()):
    """Split vessels at physical defect boundaries."""
    graph = nx.DiGraph()
    node_lookup = {}

    def node(point):
        key = tuple(np.round(point, 13))
        if key not in node_lookup:
            node_lookup[key] = len(node_lookup)
            graph.add_node(node_lookup[key], pos=list(point))
        return node_lookup[key]

    for branch, (start, end) in enumerate(toy.branches):
        p, q = np.array(toy.points[start]), np.array(toy.points[end])
        length = np.linalg.norm(q - p)
        limits = [0.0, 0.5, 1.0]
        lo, hi = -1.0, -1.0
        if branch == defect.branch and defect.kind in {"gap", "radius", "bridge"}:
            lo = max(0.001, defect.position - defect.length / (2 * length))
            hi = min(0.999, defect.position + defect.length / (2 * length))
            limits.extend([lo, hi])
        if defect.kind in {"junction", "junction_intact"} and toy.name == "y":
            lo, hi = ((1 - defect.length / (2 * length), 1.0) if branch == 0
                      else (0.0, defect.length / (2 * length)))
            limits.extend([lo, hi])
        if defect.kind in {"extra", "extra_intact"} and branch in {1, 2}:
            limits.append(defect.position)
        limits = sorted(set(np.round(limits, 12)))
        for a, b in zip(limits[:-1], limits[1:]):
            inside = lo < (a + b) / 2 < hi
            if inside and defect.kind in {"gap", "junction"}:
                continue
            factor = defect.radius_factor if inside and defect.kind in {"radius", "bridge"} else 1.0
            radius = toy.radii[branch] * factor
            graph.add_edge(node(p + a * (q - p)), node(p + b * (q - p)),
                           radius=radius, nominal_radius=toy.radii[branch],
                           branch=branch, bridge=inside and defect.kind == "bridge")
    if defect.kind == "extra" and toy.name == "y":
        starts = []
        for branch in (1, 2):
            a, b = toy.branches[branch]
            p, q = np.array(toy.points[a]), np.array(toy.points[b])
            starts.append(node(p + defect.position * (q - p)))
        graph.add_edge(*starts, radius=toy.radii[1] * defect.radius_factor,
                       nominal_radius=toy.radii[1], branch=3, bridge=True)
    return nx.convert_node_labels_to_integers(graph)


def branch_name(toy, branch):
    """Return a plain-language vessel name."""
    return "Straight vessel" if toy.name == "line" else ("Parent vessel", "Tumor branch", "Healthy branch")[branch]
