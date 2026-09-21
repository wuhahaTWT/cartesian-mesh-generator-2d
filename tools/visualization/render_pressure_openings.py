#!/usr/bin/env python3
"""Render audited static-pressure-driven channel fields and refinement errors."""
import argparse
import json
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'verification'))
import verify_native_flow as native
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
import numpy as np


def render(root,output):
    report=json.loads((root/'study.json').read_text())
    if report['valid'] is not True:raise ValueError('study is not valid')
    for name,digest in report['artifactSha256'].items():
        if native.sha256_file(root/name)!=digest:raise ValueError('retained artifact changed: '+name)
    cases=report['cases'];fine=cases[-1];prefix=root/f"n{fine['n']}";path=Path(str(prefix)+'.solver.cm2d')
    audit=native.verify_case(path,prefix,'custom',.1,1,native.argument_parser().parse_args(['--max-iterations','10000']))
    if not audit['valid']:raise ValueError(audit['issues'])
    mesh=native.read_cm2d(path);measured=native.measure(mesh,1e-11,1e-9)
    rows=native.read_cells(Path(str(prefix)+'.cells.csv'),mesh,measured,'custom')
    polygons=[[mesh.vertices[i] for i in cell.vertices] for cell in mesh.cells]
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig=plt.figure(figsize=(12,8),layout='constrained');gs=fig.add_gridspec(3,2,height_ratios=[1.4,1,1])
    ax=fig.add_subplot(gs[0,:]);collection=PolyCollection(polygons,array=np.array([r['u'] for r in rows]),cmap='viridis',edgecolors=(.1,.1,.1,.23),linewidths=.18)
    ax.add_collection(collection);ax.autoscale_view();ax.set_aspect('equal');ax.set_title(f"Actual {fine['cells']:,}-cell velocity field and mesh | static p/rho: 4.8 -> 0")
    fig.colorbar(collection,ax=ax,label='u (m/s)',shrink=.75);ax.set_xlabel('x (m)');ax.set_ylabel('y (m)')
    ax=fig.add_subplot(gs[1,0]);ys=np.linspace(0,1,201);ax.plot(6*ys*(1-ys),ys,'k--',label='analytic')
    for case in cases:
        n=case['n'];p=root/f'n{n}';m=native.read_cm2d(Path(str(p)+'.solver.cm2d'));v=native.read_cells(Path(str(p)+'.cells.csv'),m,native.measure(m,1e-11,1e-9),'custom')
        x=min({r['x'] for r in v},key=lambda x:abs(x-2));profile=sorted([r for r in v if r['x']==x],key=lambda r:r['y'])
        ax.plot([r['u'] for r in profile],[r['y'] for r in profile],'.-',ms=3,label=f"{case['cells']} cells")
    ax.legend(fontsize=8);ax.set_xlabel('u (m/s)');ax.set_ylabel('y (m)');ax.set_title('Centre profile; finite mesh error remains')
    ax=fig.add_subplot(gs[1,1]);ax.plot([r['x'] for r in rows],[r['p'] for r in rows],'.',ms=1,label='actual cell pressure');ax.plot([0,4],[4.8,0],'k--',label='analytic');ax.legend(fontsize=8);ax.set_xlabel('x (m)');ax.set_ylabel('p/rho (m2/s2)');ax.set_title('Static pressure, not total pressure')
    ax=fig.add_subplot(gs[2,0]);h=[1/c['n'] for c in cases];u=[c['velocityL2'] for c in cases];q=[c['flowRelativeError'] for c in cases]
    ax.loglog(h,u,'o-',label='velocity L2 / 1 m/s');ax.loglog(h,q,'s-',label='flow rate relative error');ax.loglog(h,[u[-1]*(v/h[-1])**2 for v in h],'k--',label='second-order guide');ax.legend(fontsize=8);ax.set_xlabel('h (m)');ax.set_title('Same rectangular geometry, three flow grids')
    ax=fig.add_subplot(gs[2,1]);ax.axis('off');ax.text(0,.96,f"Fine velocity L2: {fine['velocityL2']:.6g} m/s\nFine flow rate: {fine['flowPerDepth']:.9f} m2/s (exact 1)\nFlow rate error: {fine['flowRelativeError']*100:.4f}%\nMax pressure error: {fine['pressureMaxError']:.3g} m2/s2\nVelocity observed order: {fine['observedVelocityOrder']:.4f}\nIndependent momentum / continuity: PASS\nRestart: byte-identical; equal pressure: rest\n\nReference scope: this rectangular laminar channel.\nDoes not qualify curved ports or general CFD accuracy.",va='top',linespacing=1.55)
    fig.suptitle('Pressure-driven incompressible flow | native solver 0.4.24',fontsize=15);output.parent.mkdir(parents=True,exist_ok=True);fig.savefig(output,dpi=170);plt.close(fig)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('root',type=Path);p.add_argument('output',type=Path);a=p.parse_args();render(a.root,a.output)
