#!/usr/bin/env python3
"""Render the audited experimental coupled SST channel, not a benchmark solution."""
import argparse,csv,json,sys
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'verification'))
import verify_native_flow as native
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--study',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
a=p.parse_args();study=json.loads(a.study.read_text());case=study['cases'][-1]
if not study['valid'] or not case['valid']:raise ValueError('independent coupled audit required')
mesh=native.read_cm2d(Path(case['mesh']));prefix=case['prefix']
def rows(suffix):
    with Path(prefix+suffix).open() as f:return list(csv.DictReader(f))
cells=rows('.cells.csv');history=rows('.history.csv')
poly=[[mesh.vertices[v] for v in c.vertices] for c in mesh.cells]
comparisons=study.get('comparisons',[])
fig,axes=plt.subplots(2,3 if comparisons else 2,figsize=(16 if comparisons else 11,9),layout='constrained')
for ax,key,title,label,cmap in zip(axes.flat,['speed','nuT','omega'],
        ['Coupled mean speed','Eddy / molecular viscosity','Specific dissipation'],
        ['Speed [m/s]','nu_t / nu','omega [1/s]'],['viridis','magma','cividis']):
    values=[float(r[key])/(.001 if key=='nuT' else 1) for r in cells]
    collection=PolyCollection(poly,array=values,cmap=cmap,edgecolors='#627280',linewidths=.15,rasterized=True)
    ax.add_collection(collection);ax.autoscale();ax.set_aspect('equal');ax.set(title=title,xlabel='x [m]',ylabel='y [m]')
    fig.colorbar(collection,ax=ax,shrink=.82,label=label)
ax=axes.flat[3];x=[int(r['iteration']) for r in history]
for key,label in [('momentumResidual','Momentum'),('kCellResidual','k'),('omegaCellResidual','omega')]:
    ax.semilogy(x,[max(float(r[key]),1e-18) for r in history],label=label)
ax.grid(alpha=.2);ax.legend();ax.set(title='Residuals of the current coupled fields',
    xlabel='SIMPLE iteration',ylabel='Diagonal-scaled residual (equation-specific units)')
if comparisons:
    x=list(range(len(comparisons)));width=.34
    for ax,metric,title,label in [
            (axes.flat[4],'seconds','Same final equation gates','Solver time [s]'),
            (axes.flat[5],'updates','Less repeated turbulence work','Total nonlinear SST updates')]:
        for strategy,offset,color in [('nested',-.5,'#8998a7'),('interleaved',.5,'#247b85')]:
            values=[c[strategy+'SolveSeconds'] if metric=='seconds' else next(
                q['turbulenceUpdates'] for q in study['cases'] if q['cells']==c['cells'] and
                q['turbulenceUpdatesPerIteration']==(500 if strategy=='nested' else 1)) for c in comparisons]
            bars=ax.bar([i+offset*width for i in x],values,width,label=strategy.capitalize(),color=color)
            ax.bar_label(bars,fmt='%.2f' if metric=='seconds' else '%g',padding=3,fontsize=9)
        ax.set_xticks(x,[f"{c['cells']:,} cells" for c in comparisons]);ax.set(title=title,ylabel=label)
        ax.margins(y=.2);ax.legend();ax.grid(axis='y',alpha=.15)
fig.suptitle(f'Experimental steady SST-2003m | {len(cells):,} actual cells',fontsize=16)
fig.supxlabel('Parabolic inlet, no-slip walls, pressure outlet; nu=0.001 m2/s, inlet k=0.001 m2/s2, omega=2 1/s.\n'
    'Independent discrete-equation audit passed. Physical accuracy, y+ and grid independence are not yet qualified.',fontsize=9)
fig.savefig(a.output,dpi=150);plt.close(fig)
