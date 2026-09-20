#!/usr/bin/env python3
"""Plot measured SST isolation experiments without assigning solver acceptance."""
import argparse
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--study',type=Path,required=True)
p.add_argument('--output',type=Path,required=True)
a=p.parse_args();study=json.loads(a.study.read_text())
if study['accepted']:
    raise ValueError('this figure describes an explicitly unaccepted finest-grid study')
fig,axes=plt.subplots(1,3,figsize=(14,4.5),layout='constrained')
for name,label,color in [('frozen-half-long','Fixed carrier, relaxation 0.5','#466b93'),
                         ('frozen-full','Fixed carrier, relaxation 1.0','#ad6047')]:
    rows=study['runs'][name]['history']['samples'];x=[int(r['iteration']) for r in rows]
    for axis,key in zip(axes[:2],('kCell','omegaCell')):
        axis.plot(x,[float(r[key])/1e-9 for r in rows],label=label,color=color)
for axis,title in zip(axes[:2],('k local residual','Omega local residual')):
    axis.set(xlabel='Additional fixed-carrier updates',ylabel='Residual / original stopping limit',title=title)
    axis.set_yscale('log');axis.grid(alpha=.2);axis.legend(fontsize=8,loc='lower right')
run=study['runs']['initial-1e-4'];rows=run['history']['samples']
for key,label,limit,color in [('momentumResidual','Momentum',1e-7,'#ad6047'),
                             ('kCellResidual','k local',1e-9,'#327a67'),
                             ('omegaCellResidual','Omega local',1e-9,'#7859a0')]:
    axes[2].semilogy([int(r['iteration']) for r in rows],[float(r[key])/limit for r in rows],label=label,color=color)
axes[2].axhline(1,color='black',linestyle='--',linewidth=.8)
axes[2].set(xlabel='Coupled SIMPLE iteration',ylabel='Residual / original stopping limit',title='Higher initial k; identical inlet conditions')
axes[2].legend(fontsize=8,loc='upper right');axes[2].grid(alpha=.2)
fig.suptitle('12,800 cells | SST-2003m | Same physical inputs | NOT CONVERGED',fontsize=14)
fig.supxlabel('Fixed-carrier tests start from the saved iteration-900 state. Ratios above 1 fail their original gate.\nNeither stronger updates nor changed initial turbulence resolved this failure; no engineering qualification claimed.',fontsize=10)
a.output.parent.mkdir(parents=True,exist_ok=True)
fig.savefig(a.output,dpi=150)
