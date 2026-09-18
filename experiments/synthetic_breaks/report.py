"""Build a concise, self-contained scientific report from completed simulations."""
from .config import WORK
from pathlib import Path
import html
import json
import textwrap
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from pypdf import PdfReader, PdfWriter


def markdown_table(headers,rows):
    return '\n'.join(['| '+' | '.join(headers)+' |','| '+' | '.join(['---']*len(headers))+' |']+
                     ['| '+' | '.join(str(v) for v in row)+' |' for row in rows])


def main():
    root=WORK/'analysis'
    rows=json.loads((root/'all_cases.json').read_text())
    v=json.loads((root/'validation_summary.json').read_text())
    summary=json.loads((root/'random_break_summary.json').read_text())
    peaks=json.loads((root/'baseline_tissue_peak_over_inlet.json').read_text())
    by={(r['geometry'],r['case']):r for r in rows}
    line=[by['line',f'gap_at_{f:.1f}'] for f in [.2,.5,.8]]
    y=[by['y',c] for c in ['gap_trunk','gap_upper','gap_lower','gap_both_daughters']]
    max_balance=max(r['conservation_relative_error'] for r in rows)
    max_residual=max(r['linear_relative_residual'] for r in rows)
    max_energy=max(r['energy_relative_error'] for r in rows)
    tree=by['tree','intact']
    final=summary[-1]
    mesh_changes=[]
    for geometry,case in [('line','gap_at_0.5'),('y','gap_trunk'),('y','gap_upper')]:
        selected={r['tissue_n']:r for r in v['tissue_refinement'] if r['geometry']==geometry and r['case']==case}
        if 24 in selected and 36 in selected:
            mesh_changes.append([geometry,case,100*abs(selected[24]['relative_inlet_flux_change']-
                                                      selected[36]['relative_inlet_flux_change']),
                                 100*abs(selected[24]['tissue_relative_l2_error']-
                                          selected[36]['tissue_relative_l2_error'])])
    large_changes=[]
    for case in ['intact','seed_11_breaks_8','seed_11_breaks_32']:
        selected={r['tissue_n']:r for r in v['large_tree_refinement'] if r['case']==case}
        if 40 in selected and 48 in selected:
            large_changes.append([case,abs(selected[40]['inlet_flux']/selected[48]['inlet_flux']-1)*100,
                                  100*abs(selected[40]['relative_inlet_flux_change']-
                                           selected[48]['relative_inlet_flux_change']),
                                  100*abs(selected[40]['tissue_relative_l2_error']-
                                           selected[48]['tissue_relative_l2_error'])])
    amplitude=max((r['relative_scaling_error'] for r in v['inlet_amplitude_scaling']),default=float('inf'))
    expected_comparisons=5*2*4
    checks=dict(expected_primary_cases=42,actual_primary_cases=len(rows),
                serial_parallel_comparisons=len(v['mpi_comparisons']),
                analytic=v.get('analytic',{}).get('status')=='passed',
                conservation=max_balance<1e-8,linear_residual=max_residual<1e-9,
                energy=max_energy<1e-8,mpi_agreement=v['maximum_serial_parallel_relative_difference'] is not None and v['maximum_serial_parallel_relative_difference']<1e-10,
                amplitude_scaling=amplitude<1e-10,
                nested_break_monotonicity=v['nested_break_monotonicity_passed'],
                detached_component_balance=v['maximum_detached_net_exchange_over_inlet']<1e-8,
                toy_mesh_comparisons=len(mesh_changes),large_mesh_comparisons=len(large_changes),
                discretization_checks=len(v['discretization_checks']))
    if (len(rows)!=42 or len(v['mpi_comparisons'])<expected_comparisons or len(mesh_changes)!=3
        or len(large_changes)!=3 or len(v['discretization_checks'])!=10):
        raise RuntimeError(f'Missing study/validation results: {checks}')
    if not all(checks[k] for k in ['analytic','conservation','linear_residual','energy','mpi_agreement',
                                   'amplitude_scaling','nested_break_monotonicity','detached_component_balance']):
        raise RuntimeError(f'A scientific validation failed: {checks}')
    (root/'verification.json').write_text(json.dumps(checks,indent=2)+'\n')
    parts=[]
    parts.append('# Vessel-break sensitivity in a coupled 3D--1D model\n\nCompleted controlled pilot, 11 September 2026. '
                 'All results are dimensionless and concern synthetic geometries.')
    lead=(f'Break location has a much larger effect than its small missing length alone suggests. '
          f'For identical gaps removing 2.5% of a straight vessel, inlet-flux reductions are '
          f'{-100*line[0]["relative_inlet_flux_change"]:.2f}%, '
          f'{-100*line[1]["relative_inlet_flux_change"]:.2f}%, and '
          f'{-100*line[2]["relative_inlet_flux_change"]:.2f}% at 20%, 50%, and 80% of its length. '
          f'A Y trunk break reduces inlet flux by {-100*y[0]["relative_inlet_flux_change"]:.2f}%; '
          f'a single daughter break reduces it by about {-100*y[1]["relative_inlet_flux_change"]:.2f}%. '
          f'In the larger tree, 32 breaks give a median reduction of {final["flux_reduction_median_pct"]:.2f}% '
          f'and a five-seed range of {final["flux_reduction_min_pct"]:.2f}--{final["flux_reduction_max_pct"]:.2f}%.')
    parts.append(lead)
    parts.append('## Scope and model\n\nThe solver follows the symmetric circular-average formulation in '
       'Laurino and Zunino (2019), equations (11)--(12), '
       '[DOI 10.1051/m2an/2019042](https://doi.org/10.1051/m2an/2019042). '
       'It replaces the upstream multihepatic physical model. In the unit tissue cube, '
       'the weak form is:\n\n'
       '```text\nIntegral_Omega K3 grad(u).grad(v)\n'
       ' + Integral_Lambda K1 A U_prime V_prime\n'
       ' + Integral_Lambda kappa P (U - T_r u)(V - T_r v) = 0.\n```\n\n'
       'The original inlet has U=1; the tissue exterior has u=0. All other vessel endpoints, '
       'including the new break ends, have zero axial flux. There are no source terms. '
       'K3=K1=1, radius r=0.025, area A=pi r^2, perimeter P=2 pi r, and kappa=0.05. '
       'Axial flux is q=-K1 A dU/ds. The reported outflow exits through the tissue exterior; '
       'distal vessel endpoints are sealed. This is a potential/diffusion model, without '
       'a calibration to dimensional blood-flow rates.\n\n'
       'Tissue and vessel potentials use continuous P1 finite elements. Circular averages use '
       'a separate DG1 space with branch-specific planes at junctions. Positive inlet data use '
       'explicit Dirichlet lifting. Fluxes are recovered from unconstrained-matrix reactions '
       'and independently checked against integrated wall exchange. Disconnected pieces remain '
       'coupled to tissue; neither their pressure nor their exchange is artificially set to zero.')
    parts.append('For this parameter choice, peak intact tissue potential is only '
       +', '.join(f'{100*peaks[name]:.3f}% of the inlet value in the {name}' for name in ['line','y','tree'])
       +'. Tissue feedback on the supplied vessel is therefore weak in this pilot. '
       'The response is close to a leaky one-dimensional cable, while detached vessels can still '
       'redistribute tissue potential locally. Sensitivity magnitudes are specific to these coefficients '
       'and boundary conditions; they are not universal percentages for anatomical vessels.')
    parts.append('## Controlled experiments\n\n'
       'A break deletes a finite 0.02-long interval, except in the explicit line gap-width sweep. '
       'The interval boundaries are included in the intact discretization, so retained cells and '
       'the tissue mesh are identical in each comparison. Breaks do not remove an inlet or a junction. '
       'The line and Y use 82,944 tetrahedra. The larger 63-branch, 32-leaf tree uses '
       f'{tree["tissue_cells"]:,} tetrahedra, {tree["network_cells"]:,} vessel cells and '
       f'{tree["total_dofs"]:,} total pressure unknowns.\n\n'
       'Five fixed seeds (11, 29, 47, 83, 101) each generate a uniform permutation of branch IDs. '
       'The 1, 2, 4, 8, 16 and 32-break cases use nested prefixes of that permutation. '
       'There are 42 primary solves: six line cases, five Y cases and 31 tree cases. '
       'The shared intact tree is solved once. All selected branches lose an interval at their midpoint. '
       'Individual trajectories and the seed range are shown; five seeds are exploratory and do not '
       'establish a population confidence interval. All three intact test networks are trees, so every '
       'break adds a disconnected component. Networks containing loops can retain alternative supply '
       'paths and require a separate sensitivity study.')
    toy_table=[]
    for r in rows:
        if r['geometry']=='tree' or r['case']=='intact':continue
        toy_table.append([r['geometry'],r['case'],f'{-100*r["relative_inlet_flux_change"]:.3f}',
                          f'{100*r["tissue_relative_l2_error"]:.3f}',f'{100*r["detached_fraction"]:.2f}'])
    parts.append('## Toy results\n\n'+markdown_table(['Geometry','Case','Inlet flux reduction (%)',
                 'Tissue relative L2 error (%)','Detached retained length (%)'],toy_table))
    parts.append('![Line sensitivity](line_sensitivity.png)\n\n![Line flux and exchange](line_flow.png)\n\n![Y sensitivity](y_sensitivity.png)\n\n'
                 '![Y potentials](y_pressure.png)\n\n![Tissue slices](tissue_slices.png)')
    parts.append('## Random-break results\n\n'+markdown_table(
         ['Breaks','Seeds / baselines','Median flux reduction (%)','Seed range (%)','Median tissue L2 error (%)'],
         [[r['breaks'],r['replicates'],f'{r["flux_reduction_median_pct"]:.3f}',
           f'{r["flux_reduction_min_pct"]:.3f} to {r["flux_reduction_max_pct"]:.3f}',
           f'{r["tissue_l2_error_median_pct"]:.3f}'] for r in summary]))
    parts.append('![Random-break sensitivity](random_break_sensitivity.png)\n\n![Tree potentials](tree_pressure.png)')
    parts.append('The inlet-flux response is monotone along every nested break sequence. '
       'Branch location and lost inlet connectivity explain why equal break counts can have very different '
       'effects. Detached components can still have nonzero potential and internal flow through tissue-mediated '
       'exchange, while their net wall exchange is zero. Vessel pressure and axial-flux L2 errors are '
       'computed on the retained vessel domain against the intact solution restricted to the same domain. '
       'Tissue L2 errors cover the full unchanged cube.')
    parts.append('## Numerical validation\n\n'+markdown_table(['Check','Result'],[
        ['Maximum inlet / wall / tissue-outflow imbalance',f'{max_balance:.3e} relative'],
        ['Maximum linear residual',f'{max_residual:.3e} relative'],
        ['Maximum energy-identity error',f'{max_energy:.3e} relative'],
        ['One-, two-, eight-rank comparison',f'{v["maximum_serial_parallel_relative_difference"]:.3e} maximum relative difference'],
        ['Doubling the positive inlet condition',f'{amplitude:.3e} maximum relative scaling error'],
        ['Maximum detached-component net exchange / inlet flux',f'{v["maximum_detached_net_exchange_over_inlet"]:.3e}'],
        ['Analytic-limit flux convergence',', '.join(f'{p:.3f}' for p in v['analytic']['observed_flux_orders'])+' observed orders']]))
    parts.append('The analytic limit takes K3=10^8, approaching zero tissue potential. Its vessel solution '
       'is U(s)=cosh((L-s)/ell)/cosh(L/ell), ell=sqrt(K1 A/(kappa P)). Halving vessel cell size '
       'from 0.04 to 0.005 gives the expected second-order flux convergence. '
       'Serial/MPI agreement compares both intact and broken cases, not only a communication test.')
    parts.append('Changes when increasing the toy tissue grid from 24 to 36 subdivisions per axis:\n\n'+
       markdown_table(['Geometry','Case','Flux-sensitivity change (percentage points)',
                       'Tissue-error change (percentage points)'],
                      [[a,b,f'{c:.5f}',f'{d:.5f}'] for a,b,c,d in mesh_changes]))
    parts.append('Changes in selected tree cases when increasing the grid from 40 to 48 subdivisions per axis:\n\n'+
       markdown_table(['Case','Absolute flux change (%)','Flux-sensitivity change (pp)','Tissue-error change (pp)'],
                      [[a,f'{b:.5f}',f'{c:.5f}',f'{d:.5f}'] for a,b,c,d in large_changes]))
    parts.append('The complete refinement tables include coarser grids, halved vessel cell size, and doubled '
       'circle quadrature degree in `validation_summary.json`. Integrated flux responses are more stable '
       'than spatial tissue-error measures; the latter retain visible discretization dependence. '
       'The large-tree refinement covers three selected cases, not every seed and break count.\n\n'
       '![Validation](validation.png)')
    parts.append('## Files and reproducibility\n\n'
       '* `all_cases.csv` and `random_break_summary.csv`: response tables.\n'
       '* `break_intervals.csv`: every removed interval with branch ID and 3D endpoint coordinates.\n'
       '* Raw runs stay under `STUDY_WORK/runs/` on NOTS, with full per-case solutions and break metadata.\n'
       '* Each `network.xdmf` or `tissue.xdmf` needs its companion `.h5` file. Open the XDMF file in ParaView.\n'
       '* All toy cases have complete tissue solutions. Tree tissue solutions are saved for intact and '
       'each seed at 32 breaks; the other tree cases retain complete vessel solutions and scalar metrics.\n'
       '* `all_figures.pdf`, individual PDF/PNG figures, `REPORT.pdf`, `REPORT.md`, and `REPORT.html`: shareable outputs.\n'
       '* `../src/`: solver and MPI mesh support; orchestration scripts and dependency records are in the repository root.\n\n'
       'Use a Tube filter for vessel visualization and color by vessel_pressure. Tissue flux vectors are '
       'cell data; apply Cell Data to Point Data before Stream Tracer if streamlines are desired. '
       'The native NOTS FEniCSx 0.10/PETSc/OpenMPI stack was used through a project-local Python venv. '
       'The main study used eight MPI ranks. An explicit KaHIP tissue partitioner and a shared-vertex '
       'network partitioner avoid the installed default ParMETIS failure.\n\n'
       'The input ZIP was inventoried and checksummed but not modified. Labels 2 and 3 identify inflow '
       'and outflow vessel classes. No supplied patient/phantom centerline is presented as a validated '
       'intact baseline in this synthetic pilot. The source inventory describes the 24 supplied centerlines. '
       'A later data study still needs a confirmed intact network, inlet selection, radii, units, and organ mesh.')
    markdown='\n\n'.join(parts)+'\n'
    (root/'REPORT.md').write_text(markdown)
    # A small self-contained renderer for the report's limited Markdown subset.
    blocks=[];intable=False;incode=False
    import re
    for line in markdown.splitlines():
        if line.startswith('```'):
            blocks.append('</pre>' if incode else '<pre>');incode=not incode;continue
        if incode: blocks.append(html.escape(line)+'\n');continue
        if line.startswith('|'):
            if not intable:blocks.append('<table>');intable=True
            cells=[c.strip() for c in line.strip('|').split('|')]
            if all(c=='---' for c in cells):continue
            blocks.append('<tr>'+''.join('<td>'+html.escape(c)+'</td>' for c in cells)+'</tr>');continue
        if intable:blocks.append('</table>');intable=False
        if not line: continue
        image=re.fullmatch(r'!\[(.*?)\]\((.*?)\)',line)
        if image:blocks.append(f'<figure><img src="{image[2]}" alt="{image[1]}"></figure>');continue
        if line.startswith('## '):blocks.append('<h2>'+html.escape(line[3:])+'</h2>')
        elif line.startswith('# '):blocks.append('<h1>'+html.escape(line[2:])+'</h1>')
        else:
            escaped=html.escape(line)
            escaped=re.sub(r'\[([^\]]+)\]\(([^)]+)\)',r'<a href="\2">\1</a>',escaped)
            blocks.append('<p>'+escaped+'</p>')
    style='body{max-width:1100px;margin:40px auto;padding:0 24px;font:16px/1.6 system-ui;color:#172b3a}h1,h2{line-height:1.2}table{border-collapse:collapse;font-size:14px;width:100%}td{border-bottom:1px solid #ddd;padding:8px}tr:first-child{font-weight:600;background:#eef3f7}img{width:100%}pre{background:#f4f6f8;padding:18px;overflow:auto}figure{margin:24px 0}'
    (root/'REPORT.html').write_text('<!doctype html><html><meta charset="utf-8"><title>Vessel-break sensitivity</title><style>'+style+'</style><body>'+''.join(blocks)+'</body></html>')
    # Narrative pages precede the complete vector figure collection in the PDF.
    narrative=[parts[0],lead,'MODEL: Unit cube, circular-average coupling, positive vessel inlet U=1; zero exterior tissue potential and sealed other vessel endpoints. K3=K1=1, radius=0.025, kappa=0.05. Axial flux q=-K1*pi*r^2*dU/ds. All quantities dimensionless. Retained disconnected pieces still exchange with tissue.',
               'REGIME: Intact peak tissue potential is less than 0.3% of the imposed inlet value, so feedback on the supplied vessel is weak here. The magnitudes depend on these coefficients and boundary conditions. All test geometries are trees; loops can provide alternative supply paths.',
               'TOY RESULTS:']
    narrative.extend(' | '.join(str(x) for x in row) for row in toy_table)
    narrative.append('Columns: geometry, case, inlet flux reduction (%), tissue L2 error (%), detached length (%).')
    narrative.append('RANDOM BREAKS:')
    narrative.extend(f'{r["breaks"]:2d} breaks: median flux loss {r["flux_reduction_median_pct"]:.3f}%; range {r["flux_reduction_min_pct"]:.3f}--{r["flux_reduction_max_pct"]:.3f}%; median tissue error {r["tissue_l2_error_median_pct"]:.3f}%.' for r in summary)
    narrative.extend([f'VALIDATION: Maximum conservation error {max_balance:.3e}; serial/MPI relative difference {v["maximum_serial_parallel_relative_difference"]:.3e}; analytic flux order approximately 2.',
                      'SCOPE: 42 primary synthetic cases. Five random seeds. This is a controlled dimensionless pilot, not a fitted data model. Read REPORT.html or REPORT.md for methods, refinement tables, complete qualifications, and file instructions.',
                      'REFERENCE: Laurino and Zunino (2019), DOI 10.1051/m2an/2019042.'])
    lines=[]
    for paragraph in narrative:lines.extend(textwrap.wrap(paragraph.replace('#',''),100)+[''])
    text_pdf=root/'_narrative.pdf'
    with PdfPages(text_pdf) as pdf:
        for start in range(0,len(lines),54):
            fig=plt.figure(figsize=(8.27,11.69));fig.text(.06,.96,'\n'.join(lines[start:start+54]),va='top',fontsize=9,family='DejaVu Sans',linespacing=1.55)
            pdf.savefig(fig);plt.close(fig)
    merged=PdfWriter()
    for file in [text_pdf,root/'all_figures.pdf']:
        for page in PdfReader(file).pages:merged.add_page(page)
    with (root/'REPORT.pdf').open('wb') as f:merged.write(f)
    text_pdf.unlink()
    print(lead)
    print('All required verification checks passed.')


if __name__=='__main__':main()
