"""Plot saved topology-gradient measurements."""
import numpy as np
from pathlib import Path

def plot_report(report,path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    objectives=list(dict.fromkeys(c["objective"] for c in report["curves"]))
    fig, axes = plt.subplots(len(objectives),2,figsize=(11,4*len(objectives)),constrained_layout=True,squeeze=False)
    for ax, (objective, beta) in zip(axes.flat, ((o,b) for o in objectives for b in (0,3))):
        for c in report["curves"]:
            if c["objective"] == objective and c["beta"] == beta and len(c["errors"]) == len(report["steps"]):
                color = plt.cm.viridis(c["reynolds"]/200)
                ax.loglog(report["steps"], np.maximum(c["errors"],1e-16), color=color,
                          alpha=.8, lw=1, label=f"Re={c['reynolds']}" if c["sample"] == 0 else None)
        ax.axhline(report["acceptance"], color="red", ls="--", lw=.8)
        ax.set(title=f"{objective}, projection beta={beta}", xlabel="Central difference step h", ylabel="Relative directional derivative error")
        ax.grid(alpha=.2); ax.legend(fontsize=8)
    fig.suptitle("MAC Navier-Stokes-Brinkman discrete adjoint | 3 random designs per Re\nEquation/derivative consistency only; no physical-accuracy claim")
    fig.savefig(path, dpi=145)
    plt.close(fig)

