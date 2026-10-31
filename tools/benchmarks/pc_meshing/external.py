"""Mature solvers on frozen meshes; orchestration and existing diagnostics only."""
from pathlib import Path
import argparse,sys
from study import ROOT,BASE,run,docker
sys.path.insert(0,str(ROOT/'tools/benchmarks'))
import openfoam
from postprocess import foam_resources,fluent_fields

def container_path(path):
    return '/work/'+str(path.relative_to(ROOT)).replace(chr(92),'/')

def foam(names,budget,tol,tag=''):
    for name in names:
        prefix=BASE/'meshes'/name
        mesh=Path(str(prefix)+'.solver.cm2d')
        if not mesh.exists():
            print('SKIP rejected',name,flush=True)
            continue
        folder=(BASE/'openfoam'/(f'{name}-{tol}-entrypoint'+('-'+tag if tag else ''))).resolve()
        if not (folder/'comparison.json').exists():
            openfoam.prepare('annulus',Path(str(prefix)+'.foam'),mesh,folder,'upwind',6000)
            solution=folder/'system/fvSolution'
            solution.write_text(solution.read_text(encoding='utf-8').replace('p 1e-8; U 1e-8;',f'p {tol}; U {tol};'))
            control=folder/'system/controlDict'
            control.write_text(control.read_text(encoding='utf-8').replace('writeInterval 6000','writeInterval 200'))
        cmd=['docker','run','--rm','--network','none','--cpus','4','--memory','8g',
            '-v',str(ROOT)+':/work','-w','/work','--entrypoint','/bin/bash',
            'opencfd/openfoam-run:2606','/work/tools/benchmarks/pc_meshing/openfoam_runner.sh',
            container_path(folder),str(budget)]
        run(f'openfoam-{name}-{tol}-entrypoint'+('-'+tag if tag else ''),cmd,budget+40,kind='openfoam',mesh=name,
            tolerance=float(tol),prefix=str(folder.relative_to(ROOT)),threads=1,cpuLimit=4,gpuUsed=False)
        foam_resources(folder)

def firedrake(names,budget,linear):
    for name in names:
        mesh=BASE/'meshes'/f'{name}.msh'
        folder=BASE/'firedrake'/f'{name}-{linear}'
        args=['python3','tools/benchmarks/pc_meshing/firedrake_reference.py',
            '--mesh',container_path(mesh),'--output',container_path(folder),'--linear',linear]
        run(f'firedrake-{name}-{linear}',docker(args),budget,kind='firedrake',mesh=name,
            prefix=str(folder.relative_to(ROOT)),threads=1,cpuLimit=4,gpuUsed=False,boundaryModel='analytic-trace')

def fluent(names,budget,cores,tol,coupling):
    for name in names:
        mesh=BASE/'meshes'/f'{name}.solver.cm2d'
        folder=BASE/'fluent'/f'{name}-{cores}cores-{tol}-{coupling}'
        run(f'fluent-{name}-{cores}cores-{tol}-{coupling}',
            [sys.executable,ROOT/'tools/benchmarks/pc_meshing/fluent_reference.py',
             '--mesh',mesh,'--output',folder,'--cores',cores,'--tolerance',tol,'--coupling',coupling],
            budget,kind='fluent',mesh=name,prefix=str(folder.relative_to(ROOT)),
            threads=cores,gpuUsed=False,coupling=coupling)
        if (folder/'result.json').exists():
            fluent_fields(folder)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('solver',choices=['openfoam','firedrake','fluent'])
    p.add_argument('--names',nargs='+',required=True)
    p.add_argument('--budget',type=int,default=600)
    p.add_argument('--tolerance',default='1e-8')
    p.add_argument('--cores',type=int,default=2)
    p.add_argument('--linear',default='direct')
    p.add_argument('--tag',default='')
    p.add_argument('--coupling',choices=['simple','coupled'],default='simple')
    a=p.parse_args()
    if a.solver=='openfoam':foam(a.names,a.budget,a.tolerance,a.tag)
    elif a.solver=='firedrake':firedrake(a.names,a.budget,a.linear)
    else:fluent(a.names,a.budget,a.cores,a.tolerance,a.coupling)
