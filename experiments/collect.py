"""Publish completed outputs under experiments/<name>/results/."""
from pathlib import Path
import argparse
import html
import json
import re
import shutil
import textwrap

from experiments.config import ROOT
from experiments.synthetic_breaks.config import RESULTS as SYNTHETIC_RESULTS, WORK as SYNTHETIC_WORK
from experiments.pt002_inflow.config import RESULTS as SUPPLIED_RESULTS, WORK as SUPPLIED_WORK


def inline(text):
    text=html.escape(text)
    text=re.sub(r'`([^`]+)`',r'<code>\1</code>',text)
    text=re.sub(r'\*\*([^*]+)\*\*',r'<strong>\1</strong>',text)
    return re.sub(r'\[([^\]]+)\]\(([^)]+)\)',r'<a href="\2">\1</a>',text)


def render(markdown):
    blocks=[];in_code=False;in_table=False;in_list=False
    for line in markdown.splitlines():
        if line.startswith('```'):
            blocks.append('</pre>' if in_code else '<pre>');in_code=not in_code;continue
        if in_code:blocks.append(html.escape(line)+'\n');continue
        if in_table and not line.startswith('|'):blocks.append('</table>');in_table=False
        item=line.startswith(('* ','- '))
        if in_list and not item:blocks.append('</ul>');in_list=False
        if line.startswith('|'):
            if not in_table:blocks.append('<table>');in_table=True
            cells=[x.strip() for x in line.strip('|').split('|')]
            if all(re.fullmatch(r':?-+:?',x) for x in cells):continue
            blocks.append('<tr>'+''.join('<td>'+inline(x)+'</td>' for x in cells)+'</tr>');continue
        if not line:continue
        image=re.fullmatch(r'!\[([^\]]*)\]\(([^)]+)\)',line)
        if image:blocks.append('<figure><img src="'+html.escape(image[2],quote=True)+'" alt="'+html.escape(image[1],quote=True)+'"></figure>');continue
        if item:
            if not in_list:blocks.append('<ul>');in_list=True
            blocks.append('<li>'+inline(line[2:])+'</li>');continue
        match=re.match(r'^(#{1,4}) (.*)',line)
        if match:
            n=len(match[1]);blocks.append(f'<h{n}>'+inline(match[2])+f'</h{n}>')
        else:blocks.append('<p>'+inline(line)+'</p>')
    if in_table:blocks.append('</table>')
    if in_list:blocks.append('</ul>')
    style='body{max-width:1100px;margin:40px auto;padding:0 24px;font:16px/1.6 system-ui;color:#172b3a}h1,h2{line-height:1.2}table{border-collapse:collapse;width:100%;font-size:14px}td{border-bottom:1px solid #ddd;padding:8px}tr:first-child{font-weight:600;background:#eef3f7}img{max-width:100%}pre{background:#f4f6f8;padding:16px;overflow:auto}code{font-size:.9em}figure{margin:24px 0}'
    return '<!doctype html><html><head><meta charset="utf-8"><title>Segmentation sensitivity study</title><style>'+style+'</style></head><body>'+''.join(blocks)+'</body></html>'


def write_report(folder, markdown):
    (folder/'REPORT.md').write_text(markdown)
    (folder/'REPORT.html').write_text(render(markdown))
    for target in re.findall(r'(?:href|src)="([^"]+)"',(folder/'REPORT.html').read_text()):
        if '://' not in target and not target.startswith('#'):
            assert (folder/target).exists(),target


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--from-results',type=Path,help='Import the historical archive layout (analysis/ and supplied/pt002_inflow/)')
    args=parser.parse_args()
    analysis=args.from_results/'analysis' if args.from_results else SYNTHETIC_WORK/'analysis'
    supplied=args.from_results/'supplied/pt002_inflow' if args.from_results else SUPPLIED_WORK
    assert json.loads((analysis/'verification.json').read_text())['actual_primary_cases']==42
    assert json.loads((supplied/'verification.json').read_text())['status']=='passed'
    for folder in [SYNTHETIC_RESULTS,SUPPLIED_RESULTS]:folder.mkdir(parents=True,exist_ok=True)
    for file in analysis.iterdir():
        if file.is_file() and not file.name.startswith('REPORT') and file.suffix in ['.png','.pdf','.csv','.json']:
            shutil.copy2(file,SYNTHETIC_RESULTS/file.name)
    shutil.copy2(analysis/'REPORT.pdf',SYNTHETIC_RESULTS/'REPORT.pdf')
    main=(analysis/'REPORT.md').read_text()
    main=main[:main.index('## Files and reproducibility')]+'''## Files and reproducibility

* This folder contains the synthetic-break study's figures, measurements, and verification records.
* `all_cases.csv`, `random_break_summary.csv`, and `break_intervals.csv` record responses and every gap definition.
* Experiment drivers are in the parent `synthetic_breaks/` directory; reusable model and mesh code is in `src/segmentation_sensitivity/` at the repository root.
* See [the repository README](../../../README.md) for environment setup, Slurm jobs, and reproduction commands.
* Full raw meshes remain in the historical NOTS runs and archive. New raw outputs use this experiment's `results/raw/` folder unless `STUDY_WORK` selects scratch storage.
* The supplied-network experiment has its own [report](../../pt002_inflow/results/REPORT.html).
'''
    for name in ['coupled_solution.vtu','network_solution.vtu','tissue_solution.vtu','metrics.csv','verification.json']:
        shutil.copy2(supplied/name,SUPPLIED_RESULTS/name)
    shutil.copy2(supplied/'input/metadata.json',SUPPLIED_RESULTS/'input_metadata.json')
    supplement=(supplied/'REPORT.md').read_text()
    for old,new in {
        'source scripts/environment.sh':'source environment.sh',
        'source src/environment.sh':'source environment.sh',
        'sbatch scripts/supplied.slurm':'python3 submit.py supplied',
        'sbatch src/supplied.slurm':'python3 submit.py supplied',
        'python -m study.supplied_results':'python -m experiments.pt002_inflow.report',
        'python -m src.supplied_results':'python -m experiments.pt002_inflow.report',
        'python -m study.supplied':'python -m experiments.pt002_inflow.run',
        'python -m src.supplied':'python -m experiments.pt002_inflow.run',
        'python -m study.package':'python package.py',
        'python -m src.package':'python package.py',
    }.items():supplement=supplement.replace(old,new)
    supplement=supplement.replace('python package.py','python -m experiments.collect').replace('sbatch supplied.slurm','python3 submit.py supplied')
    supplement=supplement.replace('The accompanying native XDMF/HDF5 files retain normalized coordinates.','Native XDMF/HDF5 files in the raw run directory retain normalized coordinates.')
    supplement=supplement.replace("Input checksums, coordinate transformation, source IDs, unmodified extracted inputs, native solutions, and both runs' metrics are included for reproducibility.","Input checksums and coordinate transformation are in `input_metadata.json`; source IDs are in the VTU files. Both runs' metrics are in `metrics.csv`. Extracted inputs and native solutions remain in the raw run directory on NOTS.")
    # Write both reports before validating their cross-links.
    (SUPPLIED_RESULTS/'REPORT.html').write_text(render(supplement))
    write_report(SYNTHETIC_RESULTS,main)
    write_report(SUPPLIED_RESULTS,supplement)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    lines=[]
    for line in supplement.splitlines():
        if line.startswith('```'):continue
        clean=re.sub(r'\[([^\]]+)\]\(([^)]+)\)',r'\1',line).replace('`','').replace('**','')
        lines.extend(textwrap.wrap(clean,100) or [''])
    with PdfPages(SUPPLIED_RESULTS/'REPORT.pdf') as pdf:
        for start in range(0,len(lines),52):
            fig=plt.figure(figsize=(8.27,11.69))
            fig.text(.06,.96,'\n'.join(lines[start:start+52]),va='top',fontsize=9,family='DejaVu Sans',linespacing=1.55)
            pdf.savefig(fig);plt.close(fig)
    print('Published reports, figures, and meshes in:')
    print(SYNTHETIC_RESULTS)
    print(SUPPLIED_RESULTS)


if __name__=='__main__':main()
