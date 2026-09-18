"""Read-only inventory of the provided archive; no clinical model is inferred."""
from .config import RESULTS, ARCHIVE
from pathlib import Path
import hashlib
import json
import zipfile
import xml.etree.ElementTree as ET
import numpy as np
import networkx as nx


def main():
    data=ARCHIVE
    root=RESULTS
    root.mkdir(parents=True,exist_ok=True)
    with zipfile.ZipFile(data) as z:
        infos=z.infolist()
        records=[]
        for info in infos:
            if '/thermoembo_run/' not in info.filename or not info.filename.endswith('_centerline.vtp'):
                continue
            doc=ET.fromstring(z.read(info))
            pts_node=doc.find('.//Points/DataArray')
            if pts_node.get('format')!='ascii':
                raise ValueError('Expected ASCII point coordinates in this archive')
            coords=np.fromstring(pts_node.text,sep=' ').reshape(-1,3)
            arrays={a.get('Name'):np.fromstring(a.text or '',sep=' ',dtype=np.int64)
                    for a in doc.findall('.//Lines/DataArray')}
            G=nx.Graph()
            G.add_nodes_from(range(len(coords)))
            start=0
            for stop in arrays['offsets']:
                conn=arrays['connectivity'][start:stop]
                G.add_edges_from(zip(conn[:-1],conn[1:])); start=stop
            components=list(nx.connected_components(G))
            records.append(dict(path=info.filename,n_points=len(coords),n_edges=G.number_of_edges(),
                           n_components=len(components),n_isolated_vertices=len(list(nx.isolates(G))),
                           n_endpoints=sum(d==1 for _,d in G.degree()),
                           cycle_rank=G.number_of_edges()-G.number_of_nodes()+len(components),
                           coordinate_min=coords.min(axis=0).tolist(),coordinate_max=coords.max(axis=0).tolist(),
                           point_fields=[a.get('Name') for a in doc.findall('.//PointData/DataArray')],
                           cell_fields=[a.get('Name') for a in doc.findall('.//CellData/DataArray')]))
        summary=dict(archive=str(data),archive_bytes=data.stat().st_size,
                     entries=len(infos),expanded_bytes=sum(i.file_size for i in infos),
                     label_meanings={'0':'background','1':'liver parenchyma','2':'inflow vessels',
                                     '3':'outflow vessels','4':'tumor'},
                     selected_for_simulation=False,
                     note='This initial study uses controlled synthetic geometry. Supplied centerlines are inventoried only.',
                     centerlines=records)
    sha=hashlib.sha256()
    with data.open('rb') as f:
        for block in iter(lambda:f.read(8*1024*1024),b''): sha.update(block)
    summary['sha256']=sha.hexdigest()
    (root/'data_audit.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(f'Inventoried {len(records)} supplied centerlines; label 2=inflow, 3=outflow.')


if __name__=='__main__': main()
