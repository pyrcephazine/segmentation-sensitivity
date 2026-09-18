"""Explicit partitioners avoiding the broken default ParMETIS ABI on NOTS.

The small network is split into balanced contiguous cell ranges. Ghost cells
are added on every shared vertex, including nonmanifold bifurcations. The 3D
mesh uses the installed parallel KaHIP partitioner. No PDE or mesh is changed.
"""
import numpy as np
from dolfinx import cpp, graph, mesh


def network_partitioner(comm, nparts, cell_types, cells):
    local=np.asarray(cells[0],dtype=np.int64).reshape(-1,2)
    counts=comm.allgather(len(local))
    total=sum(counts)
    first=sum(counts[:comm.rank])
    owners=np.minimum((np.arange(first,first+len(local))*nparts)//total,nparts-1)
    all_cells=np.vstack(comm.allgather(local))
    all_owners=np.concatenate(comm.allgather(owners))
    vertex_owners={}
    for cell,owner in zip(all_cells,all_owners):
        for v in cell: vertex_owners.setdefault(int(v),set()).add(int(owner))
    destinations=[]
    offsets=[0]
    for cell,owner in zip(local,owners):
        ghosts=(vertex_owners[int(cell[0])]|vertex_owners[int(cell[1])])-{int(owner)}
        destinations.extend([int(owner),*sorted(ghosts)])
        offsets.append(len(destinations))
    return cpp.graph.AdjacencyList_int32(np.asarray(destinations,dtype=np.int32),
                                         np.asarray(offsets,dtype=np.int32))


def tissue_partitioner():
    return mesh.create_cell_partitioner(graph.partitioner_kahip(seed=1),
                                        mesh.GhostMode.shared_facet)
