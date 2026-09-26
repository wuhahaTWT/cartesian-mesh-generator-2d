#!/usr/bin/env python3
"""Plot actual native results and independently specified precision benchmarks."""
import argparse
import csv
import json
import math
from pathlib import Path
import sys
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.ticker import NullFormatter
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'verification'))
import verify_euler as e


def rows(prefix,suffix='.cells.csv'):
    with Path(str(prefix)+suffix).open() as f:return list(csv.DictReader(f))


def render(study,output):
    report=json.loads((study/'validation/verification.json').read_text());harmonic=json.loads((study/'harmonic.json').read_text())['harmonic'];app=json.loads((study/'circle-audit.json').read_text())
    root=Path(__file__).resolve().parents[2];mesh_path=root/app['mesh'];prefix=root/app['finalPrefix'];e.audit(mesh_path,prefix)
    mesh=e.geometry.read_cm2d(mesh_path);fields=rows(prefix);summary=json.loads(Path(str(prefix)+'.json').read_text())
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(2,3,figsize=(16,9.8),layout='constrained');colors={'linear':'#c77839','quadratic':'#00898d'}
    ax=axes[0,0];polygons=[[mesh.vertices[j] for j in c.vertices] for c in mesh.cells]
    collection=PolyCollection(polygons,array=[float(c['temperature']) for c in fields],cmap='inferno',edgecolors='#73889a',linewidths=.1)
    ax.add_collection(collection);ax.autoscale_view();ax.set_aspect('equal');ax.set(xlim=(-1.4,1.4),ylim=(-1.4,1.4),xlabel='x (m)',ylabel='y (m)',title=f'Actual App Cut-cell field | {len(fields)} cells')
    fig.colorbar(collection,ax=ax,shrink=.77,label='Gas temperature (K)')
    wall=max((r for r in rows(prefix,'.faces.csv') if r['kind']=='no-slip-wall'),key=lambda r:float(r['x']));wx,wy=float(wall['x']),float(wall['y'])
    detail=ax.inset_axes([.64,.62,.33,.31]);detail.add_collection(PolyCollection(polygons,array=[float(c['temperature']) for c in fields],cmap='inferno',clim=collection.get_clim(),edgecolors='#73889a',linewidths=.2));detail.set(xlim=(wx-.006,wx+.055),ylim=(wy-.05,wy+.05),xticks=[],yticks=[],title='Wall detail');detail.title.set_fontsize(8);detail.title.set_color('white')
    ax.text(.03,.025,f'Quadratic wall gradients | t = {summary["time"]:.3g} s\nShort transient, not steady qualification',transform=ax.transAxes,fontsize=8,bbox=dict(facecolor='white',alpha=.88,pad=3))
    ax=axes[0,1]
    for warped in (False,True):
        for scheme in ('linear','quadratic'):
            data=[r for r in harmonic if r['warped']==warped and r['scheme']==scheme]
            ax.loglog([r['n'] for r in data],[100*r['wallFluxRelativeL1'] for r in data],'-o' if warped else '--s',color=colors[scheme],mfc=colors[scheme] if warped else 'white',ms=5,label=scheme+' / '+('warped' if warped else 'regular'))
    ax.set(title='Non-polynomial steady heat: wall flux',xlabel='Cells per side',ylabel='Integrated wall heat L1 error (%)',xticks=[8,16,32],xticklabels=['8','16','32']);ax.xaxis.set_minor_formatter(NullFormatter());ax.grid(which='both',alpha=.15);ax.legend(fontsize=8)
    fine=[r for r in harmonic if r['n']==32 and r['warped']];ratio=fine[0]['wallFluxRelativeL1']/fine[1]['wallFluxRelativeL1']
    ax.text(.05,.05,f'32 x 32 warped: {100*fine[1]["wallFluxRelativeL1"]:.4f}%\n{ratio:.2f} times lower error',transform=ax.transAxes,fontsize=9)
    ax=axes[0,2];cold=report['coldRefinement'];n=[r['n'] for r in cold];error=[100*r['pressureRelativeLinf'] for r in cold]
    ax.loglog(n,error,'o-',color='#645b9a',label='Cold-start pressure vs continuum');ax.loglog(n,[error[0]*(n[0]/i)**2 for i in n],'--',color='#92999e',label='Second-order reference');ax.set(title=f'Whole-flow pressure | order {report["pressureOrder"]:.3f}',xlabel='Wall-normal cells',ylabel='Pressure Linf relative error (%)',xticks=n,xticklabels=list(map(str,n)));ax.xaxis.set_minor_formatter(NullFormatter());ax.legend(fontsize=8);ax.grid(which='both',alpha=.15)
    ax.text(.03,.05,'Pressure fixed by total mass and exact integral.\nPolynomial-exact temperature does not make\nthe entire flow machine-accurate.',transform=ax.transAxes,fontsize=8)
    ax=axes[1,0];yy=[i/200 for i in range(201)];ax.plot([.125*y*(1-y) for y in yy],yy,color='#273c4b',label='Analytic Couette heating')
    for level,color in [(8,'#c77839'),(16,'#645b9a'),(32,'#00898d')]:
        cells=rows(study/f'validation/cold-steady-{level}')[::2]
        ax.plot([float(c['temperature'])-1 for c in cells],[float(c['y']) for c in cells],'o',mfc='none',ms=4,color=color,label=f'{level} cells, cold start to t=60')
    ax.set(title='From rest and uniform cold gas to steady flow',xlabel='Temperature rise (benchmark units)',ylabel='y / H');ax.legend(fontsize=8);ax.grid(alpha=.15)
    ax=axes[1,1];time=report['temporal'];dt=[r['dt'] for r in time];err=[r['conservativeRmsError'] for r in time]
    ax.loglog(dt,err,'o-',color='#00898d',label='Conserved-state RMS error');ax.loglog(dt,[err[0]*(d/dt[0])**2 for d in dt],'--',color='#92999e',label='Second-order reference');ax.set(title=f'Time refinement, fixed mesh | order {report["temporalOrder"]:.3f}',xlabel='Maximum time step',ylabel='RMS difference from fine-step reference');ax.xaxis.set_minor_formatter(NullFormatter());ax.legend(fontsize=8);ax.grid(which='both',alpha=.15)
    ax.set_xticks([.00025,.0005,.001],labels=['2.5e-4','5e-4','1e-3']);ax.text(.04,.70,'Reference dt = 1.5625e-5\nCross-check dt = 3.125e-5\nNumerical time convergence',transform=ax.transAxes,fontsize=8)
    ax=axes[1,2];history=rows(study/'validation/cold-steady-16','.history.csv');tt=[];work=[];heat=[];energy=[];wi=qo=0
    for r in history:
        dt=float(r['dt']);wi-=dt*float(r['boundaryViscousWork']);qo+=dt*float(r['boundaryHeat']);tt.append(float(r['time']));work.append(wi);heat.append(qo);energy.append(float(r['totalEnergy'])-2.5)
    ax.plot(tt,work,color='#c77839',label='Mechanical energy input');ax.plot(tt,heat,'--',color='#00898d',label='Conducted heat out');ax.plot(tt,[w-q for w,q in zip(work,heat)],color='#645b9a',label='Work minus heat');stride=max(1,len(tt)//18);ax.plot(tt[::stride],energy[::stride],'o',mfc='none',ms=4,color='#273c4b',label='Computed total-energy gain')
    ax.set(title='Full cold-start energy budget',xlabel='Time (benchmark units)',ylabel='Energy per unit depth (benchmark units)');ax.legend(fontsize=8);ax.grid(alpha=.15)
    fig.suptitle('Native 2D wall accuracy: conservative heat, stress and independently checked error',fontsize=16,color='#243b4b')
    fig.supxlabel('Overall spatial/time order remains two. Constant-property ideal gas; no turbulence or general curved-wall engineering qualification.',fontsize=9)
    output.parent.mkdir(parents=True,exist_ok=True);fig.savefig(output,dpi=160);plt.close(fig)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--study',type=Path,required=True);parser.add_argument('--output',type=Path,required=True);args=parser.parse_args();render(args.study,args.output)
