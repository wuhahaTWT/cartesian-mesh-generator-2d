#!/usr/bin/env python3
"""Plot computed laminar fields against independently defined continuum solutions."""
import argparse,csv,json,math,sys
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'verification'))
import verify_euler as e

def rows(prefix,suffix='.cells.csv'):
    with Path(str(prefix)+suffix).open() as f:return list(csv.DictReader(f))

def render(study,mesh_path,prefix,output):
    report=json.loads((study/'verification.json').read_text());e.audit(mesh_path,prefix)
    mesh=e.geometry.read_cm2d(mesh_path);cells=rows(prefix);s=json.loads(Path(str(prefix)+'.json').read_text())
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axs=plt.subplots(2,2,figsize=(13.6,9.4),layout='constrained')
    ax=axs[0,0];polygons=[[mesh.vertices[j] for j in c.vertices] for c in mesh.cells]
    col=PolyCollection(polygons,array=[float(c['temperature']) for c in cells],cmap='inferno',edgecolors='#67878e',linewidths=.12)
    ax.add_collection(col);ax.autoscale_view();ax.set_aspect('equal');fig.colorbar(col,ax=ax,shrink=.8,label='Gas temperature (K)')
    ax.set(title=f'Actual App: no-slip heated cylinder | {len(cells)} cells',xlabel='x (m)',ylabel='y (m)',xlim=(-1.4,1.4),ylim=(-1.4,1.4))
    wall=max((r for r in rows(prefix,'.faces.csv') if r['kind']=='no-slip-wall'),key=lambda r:float(r['x']));wx,wy=float(wall['x']),float(wall['y'])
    detail=ax.inset_axes([.65,.62,.32,.30]);detail.add_collection(PolyCollection(polygons,array=[float(c['temperature']) for c in cells],cmap='inferno',edgecolors='#67878e',linewidths=.35,clim=col.get_clim()))
    detail.set(xlim=(wx-.006,wx+.04),ylim=(wy-.055,wy+.055),xticks=[],yticks=[],title='Wall detail');detail.title.set_fontsize(8)
    ax.text(.02,.02,f't = {s["time"]:.3g} s | mu = {s["dynamicViscosity"]:g} Pa s | k = {s["thermalConductivity"]:g} W/(m K)',transform=ax.transAxes,fontsize=8,bbox=dict(facecolor='white',alpha=.86,pad=3))
    ax=axs[0,1];yy=[i/200 for i in range(201)]
    ax.plot([.125*y*(1-y) for y in yy],yy,color='#1f303f',label='Analytic steady Couette',lw=1.5)
    for n,color in [(8,'#df9632'),(16,'#8180ad'),(32,'#00878c')]:
        p=study/f'couette-{n}';e.audit(study/f'couette-{n}.solver.cm2d',p);data=rows(p)[::4]
        ax.plot([float(c['temperature'])-1 for c in data],[float(c['y']) for c in data],'o',mfc='none',ms=4,color=color,label=f'{n} wall-normal cells')
    ax.set(title=f'Viscous heating + Fourier conduction | order {report["couetteOrder"]:.3f}',xlabel='Temperature rise above wall (test units)',ylabel='y / H')
    ax.grid(alpha=.18);ax.legend(fontsize=8)
    ax=axs[1,0];xs=[i/250 for i in range(251)];decay=math.exp(-.1*(2*math.pi)**2*.1);amplitude=1e-4*math.sqrt(1.4)
    ax.plot(xs,[decay*math.cos(2*math.pi*x) for x in xs],color='#263547',label='Continuum shear decay')
    for n,color in [(16,'#df9632'),(32,'#8180ad'),(64,'#00878c')]:
        data=rows(study/f'shear-wave-{n}')[:n]
        ax.plot([float(c['x']) for c in data],[float(c['v'])/amplitude for c in data],'.',ms=3,color=color,label=f'{n} cells / wavelength')
    order=report['refinement']['shear-wave']['orders']['v']
    ax.set(title=f'Periodic transverse shear wave | order {order:.3f}',xlabel='x / wavelength',ylabel='v / initial amplitude',xlim=(0,1));ax.grid(alpha=.18);ax.legend(fontsize=8)
    ax=axs[1,1];data=rows(study/'couette-32','.history.csv');initial=Path(study/'couette-32.checkpoint').read_text().splitlines();index=next(i for i,v in enumerate(initial) if v.startswith('STATE '))
    g=e.geometry.measure(e.geometry.read_cm2d(study/'couette-32.solver.cm2d'),1e-11,1e-10)
    e0=math.fsum(a*float(row.split()[3]) for a,row in zip(g.areas,initial[index+1:-1]));tt=[0.];work=[0.];heat=[0.];energy=[0.]
    for r in data:
        dt=float(r['dt']);tt.append(float(r['time']));work.append(work[-1]-dt*float(r['boundaryViscousWork']));heat.append(heat[-1]+dt*float(r['boundaryHeat']));energy.append(float(r['totalEnergy'])-e0)
    net=[a-b for a,b in zip(work,heat)]
    ax.plot(tt,work,color='#cc8a2e',label='Wall mechanical input')
    ax.plot(tt,heat,'--',color='#007f88',label='Outward conducted heat')
    ax.plot(tt,net,color='#5e5399',label='Input work minus heat out')
    ax.plot(tt,energy,'o',mfc='none',ms=3,color='#1f303f',markevery=max(1,len(tt)//12),label='Computed total-energy change')
    ax.set(title=f'Moving-wall energy budget | error {report["couette"][-1]["energyBalanceRelative"]:.2e}',xlabel='Time (test units)',ylabel='Energy per unit depth (test units)');ax.grid(alpha=.18);ax.legend(fontsize=8)
    fig.suptitle('Native 2D compressible laminar transport: stress, wall work and heat',fontsize=15,color='#233945')
    fig.supxlabel('Constant mu and k, ideal gas, Stokes hypothesis. Couette starts from its analytic profile; short App evolution is not steady-state or engineering qualification.',fontsize=8)
    output.parent.mkdir(parents=True,exist_ok=True);fig.savefig(output,dpi=155);plt.close(fig)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--study',type=Path,required=True);p.add_argument('--mesh',type=Path,required=True);p.add_argument('--prefix',type=Path,required=True);p.add_argument('--output',type=Path,required=True);args=p.parse_args();render(args.study,args.mesh,args.prefix,args.output)
