#!/usr/bin/env python3
"""Render a clearly unaccepted SST diagnostic state and localized equation errors."""
import argparse
import csv
import hashlib
import json
import sys
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.colors import LogNorm
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'verification'))
import verify_native_flow as native

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--study',type=Path,required=True)
p.add_argument('--output',type=Path,required=True)
a=p.parse_args();study=json.loads(a.study.read_text())
run=study['runs']['fine-state-900'];audit=study['independentDiagnostic']
if study['accepted'] or audit['valid'] or not audit['diagnosticOnly']:
    raise ValueError('explicitly unaccepted diagnostic evidence required')
def checked(name):
    f=run['files'][name];path=Path(f['path'])
    if hashlib.sha256(path.read_bytes()).hexdigest()!=f['sha256']:raise ValueError('artifact hash mismatch')
    return path
with checked('flow.unconverged.cells.csv').open() as stream:cells=list(csv.DictReader(stream))
with checked('flow.history.csv').open() as stream:history=list(csv.DictReader(stream))
mesh_path=Path(run['mesh'])
if hashlib.sha256(mesh_path.read_bytes()).hexdigest()!=run['meshSha256']:raise ValueError('mesh hash mismatch')
mesh=native.read_cm2d(mesh_path)
if len(cells)!=len(mesh.cells):raise ValueError('cell count mismatch')
polygons=[[mesh.vertices[i] for i in cell.vertices] for cell in mesh.cells]
fig,ax=plt.subplots(2,2,figsize=(13,8),layout='constrained')
ax[0,0].add_collection(PolyCollection(polygons,facecolors='#f5f6f8',edgecolors='#718292',linewidths=.5))
m=audit['momentum']['cellResidual']
points=[('Momentum',m['worstCentre']['x'],m['worstCentre']['y'],'#c44947')]
for key,color in [('k','#327a67'),('omega','#7859a0')]:
    q=audit['scalar'][key]['worstCentre'];points.append((key,*q,color))
for name,x,y,color in points:ax[0,0].scatter(x,y,s=55,color=color,label=name,zorder=3)
ax[0,0].legend(fontsize=9);ax[0,0].set_title('Worst residual cells in the actual mesh')
field=PolyCollection(polygons,array=[float(c['omega']) for c in cells],cmap='viridis',norm=LogNorm(),edgecolors='none')
ax[0,1].add_collection(field);fig.colorbar(field,ax=ax[0,1],label='Omega [1/s]',shrink=.8)
ax[0,1].set_title('Unconverged field: diagnostic use only')
for axis in ax[0]:
    axis.set(xlim=(-.01,.035),ylim=(0,.0006),xlabel='x [m]',ylabel='y [m]')
    axis.axvline(0,color='#555',linestyle='--',linewidth=.6)
iterations=[int(h['iteration']) for h in history]
for key,label,limit,color in [('momentumResidual','Momentum',1e-7,'#c44947'),('kCellResidual','k local',1e-9,'#327a67'),('omegaCellResidual','omega local',1e-9,'#7859a0')]:
    ax[1,0].semilogy(iterations,[max(float(h[key])/limit,1e-300) for h in history],label=label,color=color)
ax[1,0].axhline(1,color='black',linestyle='--',linewidth=.8)
ax[1,0].set(xlabel='SIMPLE iteration',ylabel='Residual / configured limit',title='Global-looking progress still leaves local failures')
ax[1,0].legend(fontsize=9)
last=history[-1];labels=['Momentum','k local','omega local','Continuity'];ratios=[float(last[k])/t for k,t in [('momentumResidual',1e-7),('kCellResidual',1e-9),('omegaCellResidual',1e-9),('continuity',1e-8)]]
bars=ax[1,1].bar(labels,ratios,color=['#c44947','#327a67','#7859a0','#829aa5'])
ax[1,1].set_yscale('log');ax[1,1].axhline(1,color='black',linestyle='--',linewidth=.8)
ax[1,1].bar_label(bars,labels=[f'{r:.2g}' for r in ratios],padding=3)
ax[1,1].set(ylabel='Residual / configured limit',ylim=(1e-4,1e9),title='Iteration 900: values above 1 fail their gate')
for axis in ax[1]:axis.grid(axis='y',alpha=.2);axis.set_axisbelow(True)
fig.suptitle('12,800 cells | SST-2003m | Re_plate = 10 million | NOT CONVERGED',fontsize=14)
fig.supxlabel('Independent geometry, constitutive and face checks are consistent; original convergence checks still fail.\nVertical zoom. No physical-accuracy or engineering qualification is claimed.',fontsize=10)
a.output.parent.mkdir(parents=True,exist_ok=True);fig.savefig(a.output,dpi=150)
