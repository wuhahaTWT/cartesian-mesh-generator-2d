#!/usr/bin/env python3
"""Actual Euler fields and failed accuracy gate; no synthetic flow rendering."""
import argparse,csv,json,math,sys
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'verification'))
import verify_euler as euler
import verify_native_flow as native

def main(root,destination):
    report=json.loads((root/'study.json').read_text())
    for name,digest in report['artifactSha256'].items():
        path=root/name
        if not path.resolve().is_relative_to(root.resolve()) or native.sha256_file(path)!=digest:raise ValueError('artifact changed: '+name)
    fig,axes=plt.subplots(3,2,figsize=(13,11),layout='constrained')
    xx=[i/1000 for i in range(1001)];exact=[euler.sod(x,.2) for x in xx]
    for ax,key,k in [(axes[0,0],'rho',0),(axes[0,1],'p',2)]:
        ax.plot(xx,[r[k] for r in exact],'k-',lw=1.1,label='Exact Riemann solution')
        for n in [100,200,400]:
            prefix=root/f'sod{n}';euler.audit(root/f'sod{n}.solver.cm2d',prefix)
            rows=list(csv.DictReader(Path(str(prefix)+'.cells.csv').open()));profile=rows[:n]
            ax.plot([float(r['x']) for r in profile],[float(r[key]) for r in profile],lw=.9,label=f'{n} x 4 cells')
        ax.set(xlabel='x',ylabel='Density' if key=='rho' else 'Absolute pressure',title=f'Sod shock tube: {key}, t = 0.2');ax.legend(fontsize=8);ax.grid(alpha=.2)
    for ax,name,meshname,key,title,limits in [
        (axes[1,0],'vortex160','vortex160','rho','Periodic vortex: actual density at t = 1',(6,15,6,15)),
        (axes[1,1],'circle','circle','mach','Mach 0.3 cylinder: inviscid startup at t = 0.4',(-2,3,-2,2))]:
        meshpath=root/f'{meshname}.solver.cm2d';prefix=root/name;euler.audit(meshpath,prefix)
        mesh=native.read_cm2d(meshpath);rows=list(csv.DictReader(Path(str(prefix)+'.cells.csv').open()))
        polygons=[[mesh.vertices[i] for i in cell.vertices] for cell in mesh.cells]
        coll=PolyCollection(polygons,array=[float(r[key]) for r in rows],cmap='viridis',edgecolors='#667078',linewidths=.11)
        ax.add_collection(coll);ax.set(xlim=limits[:2],ylim=limits[2:],aspect='equal',xlabel='x',ylabel='y',title=title);fig.colorbar(coll,ax=ax,shrink=.82,label=key)
    sod=[r for r in report['runs'] if r['name'].startswith('sod')]
    for field in ['rho','u','p']:
        axes[2,0].loglog([int(r['name'][3:]) for r in sod],[r['errors'][field] for r in sod],'o-',label=field)
    axes[2,0].set(xlabel='Streamwise cells',ylabel='Area-weighted L1 error',title='Sod refinement: errors decrease');axes[2,0].legend();axes[2,0].grid(alpha=.2)
    vortex=[r for r in report['runs'] if r['name'].startswith('vortex')]
    axes[2,1].plot([int(r['name'][6:]) for r in vortex],[100*r['densityPerturbationRelativeL2'] for r in vortex],'o-')
    axes[2,1].axhline(30,color='#bc382b',linestyle='--',label='Declared 30% fine-error target')
    axes[2,1].set(xlabel='Cells per side',ylabel='Density perturbation relative L2 [%]',title='Vortex accuracy study remains FAILED');axes[2,1].legend(fontsize=8);axes[2,1].grid(alpha=.2)
    fig.suptitle('Experimental native ideal-gas Euler: shared-face conservation, first-order accuracy limits',fontsize=14)
    fig.supxlabel('No viscosity or heat conduction. 128 grid misses 30% vortex target; 160 reaches 26.8% but misses original pair-ratio gate. App integration pending.',fontsize=9)
    destination.parent.mkdir(parents=True,exist_ok=True);fig.savefig(destination,dpi=145);plt.close(fig)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--study',type=Path,required=True);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args();main(args.study,args.output)
