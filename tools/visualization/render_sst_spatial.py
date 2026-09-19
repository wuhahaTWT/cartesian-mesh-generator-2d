#!/usr/bin/env python3
"""Plot actual native Cut-cell wall distances and reconstructed SST coefficients."""
import argparse,csv,json,sys
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'verification'))
import verify_native_flow as native
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--study',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
a=p.parse_args();study=json.loads(a.study.read_text());r=study['cases']['shear'][-1]
if not study['valid']:raise ValueError('audit required')
mesh=native.read_cm2d(Path(r['mesh']))
with Path(r['prefix']+'.cells.csv').open() as f:rows=list(csv.DictReader(f))
polygons=[[mesh.vertices[i] for i in c.vertices] for c in mesh.cells]
fig,axes=plt.subplots(1,3,figsize=(14,5),layout='constrained')
for ax,key,title,cmap,label in zip(axes,['distance','F1','nuT'],['Nearest selected wall segment','Reconstructed SST blending function','Spatial turbulent viscosity'],['viridis','cividis','magma'],['Distance [m]','F1','nu_t [m2/s]']):
 values=[float(row[key]) for row in rows]
 collection=PolyCollection(polygons,array=values,cmap=cmap,edgecolors='#697681',linewidths=.13,rasterized=True)
 ax.add_collection(collection);ax.autoscale();ax.set_aspect('equal');ax.set(title=title,xlabel='x [m]',ylabel='y [m]');fig.colorbar(collection,ax=ax,shrink=.75,label=label)
axes[0].plot([0,1],[0,0],lw=3,color='#bd6a3b',label='Selected wall');axes[0].legend(fontsize=8)
fig.suptitle(f'SST spatial foundation on {len(rows):,} actual Cut-cells',fontsize=16)
fig.supxlabel('Affine prescribed k, omega and U=(y,0); actual segment distance and reconstructed tensor strain.\n'
               'Six frozen transport cases independently audited. This is not a physical wall-flow or converged RANS solution.',fontsize=10)
fig.savefig(a.output,dpi=150);plt.close(fig)
