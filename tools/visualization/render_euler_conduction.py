#!/usr/bin/env python3
"""Render actual coupled heat fields, energy accounting and continuum references."""
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
import verify_euler as euler
import verify_native_flow as geometry
from verify_heat_conduction import linear_euler_fourier_mode


def rows(path):
    with path.open() as file:return list(csv.DictReader(file))


def render(study,mesh_path,prefix,destination):
    evidence=json.loads((study/'conduction-validation.json').read_text())
    euler.audit(mesh_path,prefix)
    mesh=geometry.read_cm2d(mesh_path);data=rows(Path(str(prefix)+'.cells.csv'))
    summary=json.loads(Path(str(prefix)+'.json').read_text())
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(2,2,figsize=(13.4,9.5),layout='constrained')
    ax=axes[0,0]
    polygons=[[mesh.vertices[i] for i in c.vertices] for c in mesh.cells]
    collection=PolyCollection(polygons,array=[float(r['temperature']) for r in data],cmap='inferno',edgecolors='#7d9aa2',linewidths=.28)
    ax.add_collection(collection);ax.autoscale_view();ax.set_aspect('equal')
    ax.set(xlabel='x (m)',ylabel='y (m)',title=f'Actual App result | {len(mesh.cells)} cells | t = {summary["time"]:.5g} s')
    ax.text(.02,.04,f'Wall: {summary["wallValue"]:g} K; k = {summary["thermalConductivity"]:g} W/(m K)',transform=ax.transAxes,color='white',fontsize=9,bbox=dict(facecolor='#202833',alpha=.78,pad=4))
    fig.colorbar(collection,ax=ax,shrink=.76,label='Gas temperature (K)')

    ax=axes[0,1];history=rows(study/'flux.history.csv');euler.audit(study/'box.solver.cm2d',study/'flux')
    tt=[0.];energy=[0.];heat=[0.];energy0=.1/(1.4-1)
    for row in history:
        tt.append(float(row['time']));energy.append(float(row['totalEnergy'])-energy0)
        heat.append(heat[-1]-float(row['dt'])*float(row['boundaryHeat']))
    ax.plot(tt,heat,color='#243544',lw=1.8,label='Integrated inward wall heat')
    ax.plot(tt,energy,'o',mfc='none',mec='#cf782c',ms=4,markevery=max(1,len(tt)//13),label='Computed total-energy increase')
    balance=evidence['walls']['flux']['historyEnergyBalanceRelative']
    ax.set(title=f'Sealed heat-flux box | budget error {balance:.2e}',xlabel='Time (test units)',ylabel='Energy gain per unit depth (test units)')
    ax.grid(alpha=.17);ax.legend(fontsize=9)

    ax=axes[1,0];mode=evidence['continuumMode'];amplitude=mode['amplitude'];end=mode['time']
    a,b,c=linear_euler_fourier_mode(end,mode['conductivity'])
    xx=[i/300 for i in range(301)]
    data=rows(study/'wave64.cells.csv')[:64];euler.audit(study/'wave64.solver.cm2d',study/'wave64')
    for key,reference,color,label,offset in [
        ('temperature',lambda x:c*math.cos(2*math.pi*x),'#007f88','Temperature perturbation',1),
        ('rho',lambda x:a*math.cos(2*math.pi*x),'#7869a9','Density perturbation',1),
        ('u',lambda x:b*math.sin(2*math.pi*x),'#d77b31','Generated acoustic velocity',0)]:
        ax.plot(xx,[reference(x)/amplitude for x in xx],color=color,label=label,lw=1.2)
        ax.plot([float(r['x']) for r in data],[(float(r[key])-offset)/amplitude for r in data],'.',color=color,ms=3)
    ax.set(title='Coupled Fourier mode | 64 x 4 cells | t = 0.08',xlabel='x / wavelength',ylabel='Perturbation / initial amplitude',xlim=(0,1))
    ax.text(.5,.05,'Lines: independent continuum linearization; dots: solver',ha='center',transform=ax.transAxes,fontsize=8)
    ax.grid(alpha=.17);ax.legend(fontsize=8,loc='upper right')

    ax=axes[1,1];runs=evidence['refinement'];nn=[r['nx'] for r in runs]
    for key,color in [('temperature','#007f88'),('rho','#7869a9'),('u','#d77b31')]:
        ax.loglog(nn,[r['normalizedL1'][key] for r in runs],'o-',color=color,label=f'{key}: order {runs[-1]["observedOrder"][key]:.2f}')
    anchor=runs[1]['normalizedL1']['temperature']
    ax.loglog(nn,[anchor*(32/n)**2 for n in nn],'--',color='#9da6ac',lw=1,label='Second-order slope')
    ax.set(title='Coupled-mode spatial refinement | amplitude 1e-5',xlabel='Cells per wavelength',ylabel='Area-weighted L1 / initial amplitude',xticks=nn)
    ax.set_xticklabels([str(n) for n in nn]);ax.xaxis.set_minor_formatter(NullFormatter());ax.grid(alpha=.17,which='both');ax.legend(fontsize=9)
    fig.suptitle('Native 2D Euler-Fourier: total-energy coupling, conservation and refinement',fontsize=15,color='#20333a')
    fig.supxlabel('Constant conductivity, ideal gas, slip walls; no viscous stress or solid coupling. Measured development evidence, not general engineering qualification.',fontsize=9)
    destination.parent.mkdir(parents=True,exist_ok=True);fig.savefig(destination,dpi=155);plt.close(fig)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--study',type=Path,required=True)
    parser.add_argument('--mesh',type=Path,required=True);parser.add_argument('--prefix',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args();render(args.study,args.mesh,args.prefix,args.output)
