"""Generate publication-friendly figures, paired tables, and validation summaries."""
from .config import WORK
from pathlib import Path
import argparse
import csv
import json
import h5py
import numpy as np
import networkx as nx
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import FuncFormatter, NullFormatter
from mpl_toolkits.mplot3d.art3d import Line3DCollection
from scipy.spatial import cKDTree


def field(path,name):
    with h5py.File(path,'r') as h:
        return h[f'/Function/{name}/0'][:].squeeze()


def mesh_arrays(path):
    with h5py.File(path,'r') as h:
        return h['/Mesh/mesh/geometry'][:],h['/Mesh/mesh/topology'][:]


def load_rows(root):
    return json.loads((root/'metrics.json').read_text())


def validation_summary(base,runs):
    report={}
    keys=['inlet_flux','wall_exchange','tissue_mean_pressure','tissue_l2_norm']
    mpi=[]
    for name in ['line','y']:
        serial_path=base/'mpi1'/f'{name}_n8'
        if not (serial_path/'complete.json').exists(): continue
        serial={r['case']:r for r in load_rows(serial_path)}
        for ranks in [2,8]:
            other_path=base/f'mpi{ranks}'/f'{name}_n8'
            if not (other_path/'complete.json').exists(): continue
            for row in load_rows(other_path):
                if row['case'] in serial:
                    for key in keys:
                        error=abs(row[key]/serial[row['case']][key]-1)
                        mpi.append(dict(geometry=name,case=row['case'],ranks=ranks,
                                        metric=key,relative_difference=error))
    report['mpi_comparisons']=mpi
    report['maximum_serial_parallel_relative_difference']=max((r['relative_difference'] for r in mpi),default=None)
    analytic=base/'analytic/results.json'
    if analytic.exists(): report['analytic']=json.loads(analytic.read_text())
    refinement=[]
    for name in ['line','y']:
        for n in [12,24,36]:
            path=(runs if n==24 else base/'mesh')/f'{name}_n{n}'
            if not (path/'complete.json').exists(): continue
            for row in load_rows(path):
                if row['case'] not in ['intact','gap_at_0.5','gap_trunk','gap_upper']: continue
                refinement.append({k:row[k] for k in ['geometry','case','tissue_n','inlet_flux',
                                   'relative_inlet_flux_change','tissue_relative_l2_error','tissue_mean_pressure']})
    report['tissue_refinement']=refinement
    large_refinement=[]
    for n in [32,40,48]:
        path=(runs if n==40 else base/'large_mesh')/f'tree_n{n}'
        if not (path/'complete.json').exists(): continue
        for row in load_rows(path):
            if row['case'] in ['intact','seed_11_breaks_8','seed_11_breaks_32']:
                large_refinement.append({k:row[k] for k in ['geometry','case','tissue_n','inlet_flux',
                       'relative_inlet_flux_change','tissue_relative_l2_error','tissue_mean_pressure']})
    report['large_tree_refinement']=large_refinement
    checks=[]
    for test in ['network_refinement','quadrature']:
        for name in ['line','y']:
            path=base/test/f'{name}_n24'
            if not (path/'complete.json').exists(): continue
            standard={r['case']:r for r in load_rows(runs/f'{name}_n24')}
            for row in load_rows(path):
                original=standard[row['case']]
                checks.append(dict(test=test,geometry=name,case=row['case'],
                              relative_inlet_flux_difference=abs(row['inlet_flux']/original['inlet_flux']-1),
                              sensitivity_difference_percentage_points=100*abs(
                                  row['relative_inlet_flux_change']-original['relative_inlet_flux_change']),
                              tissue_l2_error_difference=abs(row['tissue_relative_l2_error']-
                                                             original['tissue_relative_l2_error'])))
    report['discretization_checks']=checks
    amplitude=[]
    for name in ['line','y']:
        path=base/'amplitude2'/f'{name}_n8'
        if not (path/'complete.json').exists(): continue
        unit={r['case']:r for r in load_rows(base/'mpi1'/f'{name}_n8')}
        for row in load_rows(path):
            for key in keys:
                amplitude.append(dict(geometry=name,case=row['case'],metric=key,
                    relative_scaling_error=abs(row[key]/(2*unit[row['case']][key])-1)))
    report['inlet_amplitude_scaling']=amplitude
    return report


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--base',default=str(WORK))
    parser.add_argument('--output',default=str(WORK/'analysis'))
    parser.add_argument('--toys-only',action='store_true',help='Plot completed line and Y cases while the tree study runs')
    args=parser.parse_args()
    base=Path(args.base); runs=base/'runs'; out=Path(args.output)
    out.mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,
                         'figure.dpi':130,'savefig.dpi':180,'axes.grid':True,'grid.alpha':.2})
    colors=['#2166ac','#b2182b','#ef8a62','#67a9cf','#1b7837','#762a83']
    pdf=PdfPages(out/'all_figures.pdf')
    def save(fig,name):
        fig.savefig(out/f'{name}.png',bbox_inches='tight')
        fig.savefig(out/f'{name}.pdf',bbox_inches='tight')
        pdf.savefig(fig,bbox_inches='tight'); plt.close(fig)
    names=['line','y'] if args.toys_only else ['line','y','tree']
    for name in names:
        path=runs/f'{name}_n{24 if name!="tree" else 40}'
        if not (path/'complete.json').exists():
            raise RuntimeError(f'Simulation suite is not complete: {path}')
    suites={name:load_rows(runs/f'{name}_n{24 if name!="tree" else 40}') for name in names}
    all_rows=[r for rows in suites.values() for r in rows]
    scalar=[{k:v for k,v in r.items() if not isinstance(v,(dict,list))} for r in all_rows]
    with (out/'all_cases.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(scalar[0]));writer.writeheader();writer.writerows(scalar)
    (out/'all_cases.json').write_text(json.dumps(all_rows,indent=2)+'\n')
    break_rows=[]
    peaks={}
    for name in names:
        path=runs/f'{name}_n{24 if name!="tree" else 40}'
        config=json.loads((path/'configuration.json').read_text())
        branch_points=np.asarray(config['branch_points'])
        for row in suites[name]:
            for gap in row['breaks']:
                p0,p1=branch_points[config['branches'][gap['branch']]]
                direction=(p1-p0)/np.linalg.norm(p1-p0)
                center=p0+gap['fraction']*(p1-p0)
                start=center-.5*gap['width']*direction
                end=center+.5*gap['width']*direction
                break_rows.append(dict(geometry=name,case=row['case'],seed=row['seed'],**gap,
                       **{f'start_{axis}':float(value) for axis,value in zip('xyz',start)},
                       **{f'end_{axis}':float(value) for axis,value in zip('xyz',end)}))
        peaks[name]=float(field(path/'intact/tissue.h5','tissue_pressure').max()/config['parameters']['inlet_value'])
    with (out/'break_intervals.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(break_rows[0]));writer.writeheader();writer.writerows(break_rows)
    (out/'baseline_tissue_peak_over_inlet.json').write_text(json.dumps(peaks,indent=2)+'\n')

    fig,axes=plt.subplots(1,2,figsize=(11,4),layout='constrained')
    line=suites['line']
    for i,row in enumerate(line):
        path=runs/'line_n24'/row['case']/'network.h5'
        points,cells=mesh_arrays(path); U=field(path,'vessel_pressure')
        segments=np.stack([points[cells,0],U[cells]],axis=2)
        label=row['case'].replace('gap_at_','gap at ').replace('gap_width_','width ')
        axes[0].add_collection(LineCollection(segments,colors=colors[i],linewidths=1.8,label=label))
    axes[0].set(xlim=(.08,.92),ylim=(-.03,1.03),xlabel='x',ylabel='Vessel potential U',title='Line: pressure on retained segments')
    axes[0].legend(fontsize=8)
    positions=[r for r in line if r['case'].startswith('gap_at_')]
    axes[1].bar([.2,.5,.8],[-100*r['relative_inlet_flux_change'] for r in positions],width=.13,color=colors[1:4])
    axes[1].set(xlabel='Gap position / line length',ylabel='Inlet flux reduction (%)',
                title='Identical 0.02-long gaps',xticks=[.2,.5,.8])
    save(fig,'line_sensitivity')

    fig,axes=plt.subplots(1,2,figsize=(11,4),layout='constrained')
    for i,row in enumerate(line):
        path=runs/'line_n24'/row['case']
        xyz,cells=mesh_arrays(path/'network.h5')
        data=np.genfromtxt(path/'network_cells.csv',delimiter=',',names=True)
        centers=np.column_stack([data['x'],data['y'],data['z']])
        distance,index=cKDTree(centers).query(xyz[cells].mean(axis=1))
        if distance.max()>1e-9:raise RuntimeError('Flow CSV and mesh coordinates differ')
        label=row['case'].replace('gap_at_','gap at ').replace('gap_width_','width ')
        for ax,key in zip(axes,['signed_axial_flow','wall_exchange_per_length']):
            values=data[key][index]
            segments=np.stack([xyz[cells,0],np.repeat(values[:,None],2,axis=1)],axis=2)
            ax.add_collection(LineCollection(segments,colors=colors[i],linewidths=1.8,label=label))
            ax.autoscale();ax.set(xlabel='x',xlim=(.08,.92))
            ax.ticklabel_format(axis='y',style='sci',scilimits=(0,0))
    axes[0].set(ylabel='Signed axial flux q',title='Line: cellwise axial flux')
    axes[0].legend(fontsize=8)
    axes[1].set(ylabel='Wall exchange per unit length',title='Line: mean exchange in each cell')
    save(fig,'line_flow')

    fig,axes=plt.subplots(1,2,figsize=(11,4),layout='constrained')
    y=suites['y']
    for ax,case,title in [(axes[0],'intact','Intact Y'),(axes[1],'gap_trunk','Y with trunk gap')]:
        path=runs/'y_n24'/case/'network.h5'
        points,cells=mesh_arrays(path); U=field(path,'vessel_pressure')
        lc=LineCollection(points[cells,:2],array=U[cells].mean(axis=1),cmap='viridis',linewidths=4,clim=(0,1))
        ax.add_collection(lc);ax.autoscale();ax.set_aspect('equal');ax.set(title=title,xlabel='x',ylabel='y')
        fig.colorbar(lc,ax=ax,label='Vessel potential U',shrink=.8)
    save(fig,'y_pressure')

    fig,axes=plt.subplots(1,2,figsize=(11,4),layout='constrained')
    perturbed=y[1:]; labels=[r['case'].replace('gap_','').replace('_',' ') for r in perturbed]
    axes[0].bar(labels,[-100*r['relative_inlet_flux_change'] for r in perturbed],color=colors[1:5])
    axes[1].bar(labels,[100*r['tissue_relative_l2_error'] for r in perturbed],color=colors[1:5])
    axes[0].set(ylabel='Inlet flux reduction (%)',title='Y: global exchange response')
    axes[1].set(ylabel='Relative tissue L2 error (%)',title='Y: change in tissue potential')
    for ax in axes: ax.tick_params(axis='x',rotation=20)
    save(fig,'y_sensitivity')

    fig,axes=plt.subplots(2,3,figsize=(12,7),layout='constrained')
    for row,(name,case) in enumerate([('line','gap_at_0.5'),('y','gap_trunk')]):
        p0=runs/f'{name}_n24'/'intact/tissue.h5'; p1=runs/f'{name}_n24'/case/'tissue.h5'
        points,_=mesh_arrays(p0); values=field(p0,'tissue_pressure')
        ppoints,_=mesh_arrays(p1); changed=field(p1,'tissue_pressure'); delta=field(p1,'tissue_pressure_difference')
        data=[(points,values,'intact'),(ppoints,changed,'broken'),(ppoints,delta,'difference')]
        for col,(xyz,value,title) in enumerate(data):
            select=np.isclose(xyz[:,2],.5)
            vmax=values.max() if col<2 else max(abs(delta.min()),abs(delta.max()))
            levels=np.linspace(0,vmax,30) if col<2 else np.linspace(-vmax,vmax,31)
            im=axes[row,col].tricontourf(xyz[select,0],xyz[select,1],value[select],levels=levels,
                                         cmap='viridis' if col<2 else 'RdBu_r',extend='both')
            axes[row,col].set(title=f'{name}: {title}, z=0.5',xlabel='x',ylabel='y',aspect='equal')
            fig.colorbar(im,ax=axes[row,col],shrink=.8,
                         format=FuncFormatter(lambda value,_: f'{0.0 if abs(value)<1e-16 else value:.1e}'))
    save(fig,'tissue_slices')

    if args.toys_only:
        pdf.close()
        print(f'Wrote {len(all_rows)} completed toy cases and figures to {out}')
        return

    tree=suites['tree']; counts=np.array([0,1,2,4,8,16,32]); seeds=[11,29,47,83,101]
    fig,axes=plt.subplots(2,2,figsize=(11,8),layout='constrained')
    specs=[('relative_inlet_flux_change',-100,'Inlet flux reduction (%)'),
           ('tissue_relative_l2_error',100,'Relative tissue L2 error (%)'),
           ('detached_fraction',100,'Disconnected retained length (%)')]
    aggregated=[]
    for ax,(key,scale,label) in zip(axes.flat,specs):
        vals=[]
        for i,seed in enumerate(seeds):
            rows=[tree[0]]+[r for r in tree if r['seed']==seed]
            ys=np.array([r[key]*scale for r in rows]); vals.append(ys)
            ax.plot(counts,ys,'o-',lw=1,alpha=.6,color=colors[i],label=f'seed {seed}')
        vals=np.array(vals)
        ax.fill_between(counts,vals.min(axis=0),vals.max(axis=0),color='gray',alpha=.12,label='seed range')
        ax.plot(counts,np.median(vals,axis=0),'k--',lw=2,label='median')
        ax.set(xlabel='Number of breaks',ylabel=label,xticks=counts,xscale='symlog')
        ax.set_xticks(counts,labels=[str(c) for c in counts])
        ax.set_xlim(-.15,36)
        if key=='relative_inlet_flux_change': ax.legend(fontsize=8,ncol=2)
    scatter=axes[1,1].scatter([100*r['detached_fraction'] for r in tree[1:]],
                            [-100*r['relative_inlet_flux_change'] for r in tree[1:]],
                            c=[r['n_breaks'] for r in tree[1:]],cmap='plasma',s=40)
    axes[1,1].set(xlabel='Disconnected retained length (%)',ylabel='Inlet flux reduction (%)',
                   title='Location and connectivity explain variation')
    fig.colorbar(scatter,ax=axes[1,1],label='Break count')
    save(fig,'random_break_sensitivity')
    for count in counts:
        rows=[r for r in tree if r['n_breaks']==count]
        loss=np.array([-100*r['relative_inlet_flux_change'] for r in rows])
        error=np.array([100*r['tissue_relative_l2_error'] for r in rows])
        aggregated.append(dict(breaks=int(count),replicates=len(rows),flux_reduction_median_pct=float(np.median(loss)),
                              flux_reduction_min_pct=float(loss.min()),flux_reduction_max_pct=float(loss.max()),
                              tissue_l2_error_median_pct=float(np.median(error))))
    with (out/'random_break_summary.csv').open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=list(aggregated[0]));writer.writeheader();writer.writerows(aggregated)
    (out/'random_break_summary.json').write_text(json.dumps(aggregated,indent=2)+'\n')

    worst=min((r for r in tree if r['n_breaks']==32),key=lambda r:r['inlet_flux'])
    fig=plt.figure(figsize=(12,5),layout='constrained')
    for j,row in enumerate([tree[0],worst]):
        ax=fig.add_subplot(1,2,j+1,projection='3d')
        path=runs/'tree_n40'/row['case']/'network.h5'
        xyz,cells=mesh_arrays(path); U=field(path,'vessel_pressure')
        lc=Line3DCollection(xyz[cells],array=U[cells].mean(axis=1),cmap='viridis',linewidths=2,clim=(0,1))
        ax.add_collection3d(lc);ax.set(xlim=(0,1),ylim=(0,1),zlim=(0,1),xlabel='x',ylabel='y',zlabel='z',
                                       title='Intact tree' if j==0 else f'32 breaks, seed {worst["seed"]}')
        ax.view_init(elev=22,azim=-60)
        fig.colorbar(lc,ax=ax,shrink=.65,label='Vessel potential U')
    save(fig,'tree_pressure')

    validation=validation_summary(base/'validation',runs)
    validation['nested_break_monotonicity_passed']=all(
        np.max(np.diff([tree[0]['inlet_flux']]+[r['inlet_flux'] for r in tree if r['seed']==seed]))<1e-10
        for seed in seeds)
    detached_checks=[]
    for row in all_rows:
        if row['n_components']<2: continue
        n=40 if row['geometry']=='tree' else 24
        casepath=runs/f'{row["geometry"]}_n{n}'/row['case']
        xyz,conn=mesh_arrays(casepath/'network.h5')
        G=nx.Graph();G.add_edges_from(conn)
        config=json.loads((casepath.parent/'configuration.json').read_text())
        inlet=np.asarray(config['branch_points'][0])
        inlet_id=int(np.argmin(np.linalg.norm(xyz-inlet,axis=1)))
        comp={v:i for i,nodes in enumerate(nx.connected_components(G)) for v in nodes}
        supplied=comp[inlet_id]
        csvdata=np.genfromtxt(casepath/'network_cells.csv',delimiter=',',names=True)
        centers=np.column_stack([csvdata['x'],csvdata['y'],csvdata['z']])
        distance,match=cKDTree(centers).query(xyz[conn].mean(axis=1))
        if distance.max()>1e-9: raise RuntimeError('Cell CSV and XDMF geometry do not agree')
        totals=np.zeros(row['n_components'])
        for edge,index in zip(conn,match):
            totals[comp[edge[0]]]+=csvdata['length'][index]*csvdata['wall_exchange_per_length'][index]
        detached=np.delete(totals,supplied)
        detached_checks.append(dict(geometry=row['geometry'],case=row['case'],
            maximum_detached_net_exchange_over_inlet=float(np.max(np.abs(detached))/row['inlet_flux'])))
    validation['detached_component_balances']=detached_checks
    validation['maximum_detached_net_exchange_over_inlet']=max(
        (r['maximum_detached_net_exchange_over_inlet'] for r in detached_checks),default=0.)
    (out/'validation_summary.json').write_text(json.dumps(validation,indent=2)+'\n')
    ref=validation['tissue_refinement']
    fig,axes=plt.subplots(1,2,figsize=(11,4),layout='constrained')
    for name,case in [('line','gap_at_0.5'),('y','gap_trunk'),('y','gap_upper')]:
        rows=sorted([r for r in ref if r['geometry']==name and r['case']==case],key=lambda r:r['tissue_n'])
        if rows:
            axes[0].plot([r['tissue_n'] for r in rows],[-100*r['relative_inlet_flux_change'] for r in rows],'o-',label=f'{name}: {case}')
    axes[0].set(xlabel='Tissue subdivisions per axis',ylabel='Inlet flux reduction (%)',title='Sensitivity versus tissue resolution')
    axes[0].legend(fontsize=8)
    if 'analytic' in validation:
        a=validation['analytic']['samples']
        h=np.array([r['h'] for r in a]);e=np.array([r['relative_flux_error'] for r in a])
        axes[1].loglog(h,e,'o-',label='computed error')
        axes[1].loglog(h,e[0]*(h/h[0])**2,'k--',label='second-order reference')
        axes[1].set(xlabel='Vessel cell size',ylabel='Relative analytic inlet-flux error',title='Closed-form limit check')
        axes[1].set_xticks(sorted(h),labels=[f'{value:g}' for value in sorted(h)])
        axes[1].xaxis.set_minor_formatter(NullFormatter())
        axes[1].legend()
    save(fig,'validation')
    pdf.close()
    print(f'Wrote {len(all_rows)} cases and figures to {out}')


if __name__=='__main__': main()
