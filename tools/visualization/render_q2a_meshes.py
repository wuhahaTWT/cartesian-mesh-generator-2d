#!/usr/bin/env python3
"""Render the actual five Q2-A solver artifacts, no synthetic mesh illustration."""
import argparse
import json
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"verification"))
from generate_q0_baselines import read_cm2d
from generate_q1_baselines import CASES
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

p=argparse.ArgumentParser()
p.add_argument("--evidence",type=Path,default=Path("build-q2a/evidence"))
p.add_argument("--output",type=Path,default=Path("artifacts/q2a/solver-meshes.png"))
a=p.parse_args()
fig,axes=plt.subplots(2,3,figsize=(15,9))
for ax,name in zip(axes.flat,CASES):
    mesh=read_cm2d(a.evidence/f"shared-{name}.hybrid.solver.cm2d")
    report=json.loads((a.evidence/f"shared-{name}.hybrid.quality-contract.json").read_text())
    value=report["ordinary_metrics"]["face_length_over_local_background_h"]["worst"]
    ax.add_collection(LineCollection([[mesh.vertices[e.v0],mesh.vertices[e.v1]] for e in mesh.edges],
                                    colors="#50748a",linewidths=.3))
    ax.add_collection(LineCollection([[mesh.vertices[e.v0],mesh.vertices[e.v1]] for e in mesh.edges if e.patch==1],
                                    colors="#d47b24",linewidths=1.4))
    ax.autoscale();ax.set_aspect("equal")
    status="PASS" if value>=.01 else "FAIL (existing; Q2-B)"
    ax.set_title(f"{name} | {len(mesh.cells)} solver cells\nmin(face/local_h)={value:.7g} | {status}",fontsize=10)
    ax.set_xlabel("x");ax.set_ylabel("y")
axes[1,2].axis("off")
axes[1,2].text(.05,.9,"Q2-A: shared construction\n\nActual exported solver meshes\nOrange: embedded solid wall\nBlue: fluid cell edges\n\nQ1 face/local_h hard limit: 0.01\n\nNarrow gap / sharp trailing edge\nstill need constrained local repair.\n\nThis is not full Q2 acceptance.",
               transform=axes[1,2].transAxes,va="top",fontsize=12,linespacing=1.6)
fig.suptitle("Native 2D — Q2-A shared intersection / common edge partition",fontsize=15)
fig.tight_layout(rect=(0,0,1,.96))
a.output.parent.mkdir(parents=True,exist_ok=True)
fig.savefig(a.output,dpi=170)
