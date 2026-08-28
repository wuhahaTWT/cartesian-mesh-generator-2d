#!/usr/bin/env python3
"""Render real before/after CM2D solver meshes, including the Q1 micro face."""
import argparse
import json
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"verification"))
from generate_q0_baselines import read_cm2d
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection, PolyCollection

def main():
    p=argparse.ArgumentParser()
    p.add_argument("--before",type=Path,default=Path("build-q2/before"))
    p.add_argument("--after",type=Path,default=Path("build-q2/after"))
    p.add_argument("--output",type=Path,default=Path("artifacts/q2/superellipse-solver.png"))
    a=p.parse_args()
    before=read_cm2d(a.before/"superellipse.hybrid.solver.cm2d")
    after=read_cm2d(a.after/"superellipse.hybrid.solver.cm2d")
    report=json.loads((a.after/"superellipse.hybrid.quality-contract.json").read_text())
    min_face=report["ordinary_metrics"]["face_length_over_local_background_h"]["worst"]
    fig,axes=plt.subplots(1,3,figsize=(15,5),gridspec_kw={"width_ratios":[1.8,1,1]})
    fig.subplots_adjust(top=.78,bottom=.13,left=.05,right=.98,wspace=.3)
    polys=[[after.vertices[i] for i in c.vertices] for c in after.cells]
    axes[0].add_collection(PolyCollection(polys,facecolors="#e5f0f6",edgecolors="#547484",linewidths=.35))
    axes[0].autoscale();axes[0].set_aspect("equal")
    axes[0].set_title(f"Actual solver mesh: {len(after.cells)} cells\nmin(face/local_h) = {min_face:.7g} > 0.01")
    axes[0].set_xlabel("x");axes[0].set_ylabel("y")
    origin=(-2.326100423965795,0.)
    for ax,mesh,title in zip(axes[1:],[before,after],["Before: 9.79753e-9 internal face","After: one canonical grid-line vertex"]):
        lines=[]
        for e in mesh.edges:
            pts=[mesh.vertices[e.v0],mesh.vertices[e.v1]]
            if any(abs(p[0]-origin[0])<.01 and abs(p[1])<.01 for p in pts):
                lines.append([((p[0]-origin[0])*1.e9,p[1]*1.e9) for p in pts])
        ax.add_collection(LineCollection(lines,colors="#32546a",linewidths=1.3))
        local=[((p[0]-origin[0])*1.e9,p[1]*1.e9) for p in mesh.vertices
               if abs(p[0]-origin[0])<2.e-8 and abs(p[1])<3.e-8]
        ax.scatter([p[0] for p in local],[p[1] for p in local],color="#c74631",s=40,zorder=4)
        ax.set_xlim(-14,14);ax.set_ylim(-12,24);ax.set_aspect("equal")
        ax.set_title(title,fontsize=10)
        ax.set_xlabel("x offset / 1e-9");ax.set_ylabel("y / 1e-9")
        ax.grid(alpha=.2)
    fig.suptitle("Q2 superellipse — original input retained; transition resampling before Cut-cell construction",fontsize=13,y=.97)
    a.output.parent.mkdir(parents=True,exist_ok=True)
    fig.savefig(a.output,dpi=180)

if __name__=="__main__": main()
