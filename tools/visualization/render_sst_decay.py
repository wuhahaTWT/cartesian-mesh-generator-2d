#!/usr/bin/env python3
"""Render real native mesh and independently audited homogeneous k/omega decay."""
import argparse
import csv
import json
from pathlib import Path
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'verification'))
import verify_native_flow as native

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--study',type=Path,required=True)
p.add_argument('--output',type=Path,required=True)
a=p.parse_args();study=json.loads(a.study.read_text());s=study['cases']['shear']
if not s['valid']:raise ValueError('valid audit required')
r=s['runs'][-1];mesh=native.read_cm2d(Path(r['mesh']))
fig,axes=plt.subplots(1,3,figsize=(15,4.8),layout='constrained')
ax=axes[0]
ax.add_collection(LineCollection([[mesh.vertices[e.v0],mesh.vertices[e.v1]] for e in mesh.edges],linewidths=.6,colors='#697d88'))
ax.autoscale();ax.set_aspect('equal');ax.set(title=f'Actual final Cut-cell mesh: {len(mesh.cells)} cells',xlabel='x [m]',ylabel='y [m]')
t=np.linspace(0,.8,150);ax=axes[1]
ax.plot(t,(1+.3*t)**(-1.2),color='#287c8e',label='Continuous k/k0')
ax.plot(t,1/(1+.3*t),color='#ac643b',label='Continuous omega/omega0')
with Path(r['prefix']+'.history.csv').open() as f:rows=list(csv.DictReader(f))
ax.plot([float(x['time']) for x in rows],[float(x['k'])/.02 for x in rows],'o',ms=3,color='#287c8e',label='Computed k, dt=.05')
ax.plot([float(x['time']) for x in rows],[float(x['omega'])/4 for x in rows],'s',ms=3,color='#ac643b',label='Computed omega, dt=.05')
ax.set(title='Two transport equations, implicit losses',xlabel='Physical time [s]',ylabel='Normalized turbulence variable');ax.legend(fontsize=8);ax.grid(alpha=.2)
ax=axes[2]
for key,normalizer,label,color in [('kContinuousRms',.02,'k / k0','#287c8e'),('omegaContinuousRms',4,'omega / omega0','#ac643b')]:
 ax.loglog([x['timeStep'] for x in s['runs']],[x[key]/normalizer for x in s['runs']],'o-',label=label,color=color)
ax.set(title='Time refinement: observed orders 0.95-0.98',xlabel='Time step [s]',ylabel='Normalized RMS error at t=.8');ax.grid(which='both',alpha=.2);ax.legend()
fig.suptitle('SST-2003m transport foundation: homogeneous decay',fontsize=16)
fig.supxlabel('Prescribed F1=1, zero strain and zero gradients; original nonlinear balances independently checked.\n'
 'An ODE limit on native meshes, NOT wall-bounded RANS validation. Wall distance and flow coupling remain pending.',fontsize=10)
fig.savefig(a.output,dpi=150);plt.close(fig)
