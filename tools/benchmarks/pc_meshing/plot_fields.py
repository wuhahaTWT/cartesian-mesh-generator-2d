"""Scientific plot of actual accepted fields, retaining all pressure outliers."""
from pathlib import Path
import csv,json,math,sys
from study import ROOT,BASE,read
sys.path.insert(0,str(ROOT/'tools/flow'));sys.path.insert(0,str(ROOT/'tools/benchmarks'))
import native_mesh as nm
from annulus_cross import exact,errors
import matplotlib;matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.colors import SymLogNorm
import numpy as np

def main():
    out=ROOT/'artifacts/current';items=[]
    selections=[
        ('Native shape160','shape160',BASE/'solutions/shape160-coupled-adaptive-aggregation-1e-8-connected.cells.csv'),
        ('Native graded128','graded128',BASE/'solutions/graded128-coupled-adaptive-ic0-1e-8.cells.csv'),
        ('Fluent shape160','shape160',BASE/'fluent/shape160-1cores-1e-6/comparison.cells.csv'),
        ('Fluent Gmsh65k','gmsh65k',BASE/'fluent/gmsh65k-2cores-1e-8-coupled/comparison.cells.csv'),
        ('Fluent graded105k','graded-indexed512',BASE/'fluent/graded-indexed512-2cores-1e-8-coupled/comparison.cells.csv')]
    for label,name,file in selections:
        if not file.exists():continue
        mesh=nm.read_cm2d(BASE/'meshes'/f'{name}.solver.cm2d')
        with file.open(encoding='utf-8-sig') as f:rows=[{k:float(v) for k,v in row.items()} for row in csv.DictReader(f)]
        rows.sort(key=lambda r:r['cell']);assert len(rows)==len(mesh.cells)
        area=math.fsum(r['area'] for r in rows)
        gauge=math.fsum(r['area']*(r['p']-exact(r['x'],r['y'])[2]) for r in rows)/area
        span=-exact(.5,0)[2]
        values=np.array([(r['p']-exact(r['x'],r['y'])[2]-gauge)/span for r in rows])
        items.append((label,mesh,values,errors(rows)))
    if not items:return
    limit=max(float(np.max(np.abs(x[2]))) for x in items)
    fig,axes=plt.subplots(1,len(items),figsize=(4.1*len(items),4.5),layout='constrained',squeeze=False)
    for ax,(label,mesh,values,diag) in zip(axes[0],items):
        polygons=[[mesh.vertices[v] for v in c.vertices] for c in mesh.cells]
        image=PolyCollection(polygons,array=values,cmap='RdBu_r',norm=SymLogNorm(.01,vmin=-limit,vmax=limit),edgecolors='none')
        ax.add_collection(image);ax.set_aspect('equal');ax.set(xlim=(-1.02,1.02),ylim=(-1.02,1.02),xlabel='x / m',ylabel='y / m')
        ax.set_title(f"{label}\nPressure diagnostic RMS = {100*diag['pressureL2OverAnalyticPressureSpan']:.3g}%")
    fig.colorbar(image,ax=list(axes[0]),label='Pressure deviation / reference span\n(best constant removed; symlog)')
    fig.savefig(out/'pc-mesh-study-fields.png',dpi=135);plt.close(fig)

if __name__=='__main__':main()
