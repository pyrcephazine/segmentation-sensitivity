from __future__ import annotations

from dataclasses import asdict, dataclass
import math
import networkx as nx
import numpy as np


@dataclass(frozen=True)
class Break:
    branch: int
    fraction: float = 0.5
    width: float = 0.02


@dataclass
class Geometry:
    name: str
    points: np.ndarray
    branches: list[tuple[int, int]]
    inlet: int = 0

    @classmethod
    def line(cls):
        return cls('line', np.array([[.1,.5,.5],[.9,.5,.5]]), [(0,1)])

    @classmethod
    def y(cls):
        return cls('y', np.array([[.1,.5,.5],[.5,.5,.5],
                                 [.9,.75,.5],[.9,.25,.5]]), [(0,1),(1,2),(1,3)])

    @classmethod
    def tree(cls, levels=5):
        points = [[.08,.5,.5],[.18,.5,.5]]
        branches = [(0,1)]
        frontier = [(1,.15,.85,.15,.85)]
        for level in range(levels):
            following = []
            for parent, y0, y1, z0, z1 in frontier:
                if level % 2 == 0:
                    ym = (y0+y1)/2
                    boxes = [(y0,ym,z0,z1),(ym,y1,z0,z1)]
                else:
                    zm = (z0+z1)/2
                    boxes = [(y0,y1,z0,zm),(y0,y1,zm,z1)]
                for a,b,c,d in boxes:
                    child = len(points)
                    points.append([.18+.7*(level+1)/levels,(a+b)/2,(c+d)/2])
                    branches.append((parent,child))
                    following.append((child,a,b,c,d))
            frontier = following
        return cls('tree', np.array(points), branches)

    def discretize(self, max_h=.01, possible_breaks=()):
        """Put every possible gap boundary in the INTACT mesh as well.

        This makes all retained 1D cells identical across a comparison, and
        prevents a change in discretization from masquerading as break sensitivity.
        """
        graph = nx.DiGraph()
        for i, point in enumerate(self.points):
            graph.add_node(i, pos=point.tolist(), original_id=i)
        edge_id = 0
        for bid, (a,b) in enumerate(self.branches):
            p0,p1 = self.points[[a,b]]
            length = float(np.linalg.norm(p1-p0))
            knots = list(np.linspace(0,1,max(1,math.ceil(length/max_h))+1))
            for cut in possible_breaks:
                if cut.branch != bid:
                    continue
                lo,hi = cut.fraction-cut.width/(2*length),cut.fraction+cut.width/(2*length)
                if not 0 < lo < hi < 1:
                    raise ValueError(f'Gap {cut} would touch a junction or endpoint')
                knots.extend([lo,hi])
            knots = np.unique(np.round(knots,12))
            nodes = [a]
            for t in knots[1:-1]:
                idx = len(graph)
                graph.add_node(idx,pos=((1-t)*p0+t*p1).tolist(),original_id=idx)
                nodes.append(idx)
            nodes.append(b)
            for j,(u,v) in enumerate(zip(nodes[:-1],nodes[1:])):
                graph.add_edge(u,v,branch=bid,original_edge=edge_id,
                               fraction=float((knots[j]+knots[j+1])/2),
                               length=float((knots[j+1]-knots[j])*length))
                edge_id += 1
        return graph

    def perturb(self, intact, cuts):
        graph = intact.copy()
        removed = []
        for u,v,attr in list(graph.edges(data=True)):
            length = np.linalg.norm(self.points[self.branches[attr['branch']][1]]-
                                    self.points[self.branches[attr['branch']][0]])
            if any(c.branch == attr['branch'] and
                   abs(attr['fraction']-c.fraction) < c.width/(2*length)-1e-11
                   for c in cuts):
                removed.append(attr['original_edge'])
                graph.remove_edge(u,v)
        graph.remove_nodes_from(list(nx.isolates(graph)))
        if self.inlet not in graph or graph.degree(self.inlet) != 1:
            raise ValueError('The prescribed inlet must survive and remain an endpoint')
        mapping = {old:new for new,old in enumerate(sorted(graph.nodes))}
        inlet_new = mapping[self.inlet]
        graph = nx.relabel_nodes(graph,mapping,copy=True)
        components = list(nx.connected_components(graph.to_undirected()))
        supplied = next(c for c in components if inlet_new in c)
        retained = sum(d['length'] for _,_,d in graph.edges(data=True))
        detached = sum(d['length'] for u,_,d in graph.edges(data=True) if u not in supplied)
        original = sum(d['length'] for _,_,d in intact.edges(data=True))
        stats = dict(n_breaks=len(cuts),n_components=len(components),
                     retained_length=retained,removed_length=original-retained,
                     detached_length=detached,detached_fraction=detached/retained,
                     removed_cell_ids=sorted(removed),breaks=[asdict(c) for c in cuts])
        return graph, stats
