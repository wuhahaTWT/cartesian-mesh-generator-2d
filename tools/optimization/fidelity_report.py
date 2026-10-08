#!/usr/bin/env python3
"""Rebuild the compact T01 evidence, all-pair ranking and scientific figures."""
import argparse
from collections import Counter
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess

import numpy as np

from fidelity_metrics import ranking
from optimize_flow import write_json
from topology_cases import CASES
from topology_artifacts import contours,signed_area

ROOT=Path(__file__).resolve().parents[2]
STUDY=ROOT/"outputs/topology-fidelity"
SOURCE_PATHS=[".github/workflows/flow-research.yml","apps/cartmesh2d_flow_cli.cpp",
    "include/cartmesh2d/fv/Incompressible2D.hpp","include/cartmesh2d/fv/FlowCheckpoint2D.hpp",
    "src/fv/Incompressible2D.cpp","src/fv/FlowBoundaryIO2D.cpp","tests/flow_boundary_test.cpp",
    "tests/flow_topology_test.py","tests/navier_stokes_topology_test.py","tests/topology_fidelity_test.py",
    "tools/flow/native_mesh.py"]


def load(path):
    return json.loads(path.read_text())


def evidence(path):
    return dict(path=str(path.relative_to(ROOT)),sha256=hashlib.sha256(path.read_bytes()).hexdigest(),bytes=path.stat().st_size)


def compact(row):
    keys=("id","case","reynolds","darcyPort","alphaMax","filterScale","initialization","seed","path","status",
          "issue","optimizerStatus","optimizationConverged","finalParameters","porousObjective","leakage","grid","gap","gapInterval",
          "nativeStatus","nativeIssue","native","nativeAttempts","nativeRecovery","engineeringReferences","engineeringBest",
          "samePorousArea","porousReductionVsEngineering","densitySha256")
    result={k:row[k] for k in keys if k in row}
    result["projectedKkt"]=row.get("final",{}).get("projectedKkt")
    result["problem"]=row.get("problem")
    attempts=row.get("nativeAttempts",row.get("native",[]))
    result["allRecordedNativeAttemptsAccepted"]=bool(attempts) and all(a.get("metrics") is not None for a in attempts)
    if "sharpProjection" in row:
        p=row["sharpProjection"]
        result["sharpProjection"]={k:p[k] for k in ("matching","sourceFieldSha256")}
        result["sharpProjection"]["extraction"]={k:p["extraction"].get(k) for k in
                ("components","fluidArea","holes","boundarySha256","portConnectivity")}
    return result


def counts(rows):
    gaps=[r["gap"] for r in rows if "gap" in r]
    stats=dict(designs=len(rows),optimizationStatuses=dict(Counter(r["status"] for r in rows)),
        nativeMethods=dict(Counter(r.get("grid",{}).get("method","extraction-failed") for r in rows)),
        withTrueObjective=sum(r.get("grid",{}).get("J_T") is not None for r in rows),
        withAtLeastTwoGrids=sum(len(r.get("grid",{}).get("acceptedLevels",[]))>=2 for r in rows),
        componentAttemptStatuses=dict(Counter(s for r in rows for a in r.get("nativeAttempts",r.get("native",[]))
                                             for s in a.get("attemptStatuses",[a.get("status","unknown")]))),
        gapRange=[min(gaps),max(gaps)] if gaps else None)
    for name in ("solidFluxFraction","solidDissipationFraction"):
        values=[r["leakage"][name] for r in rows if "leakage" in r]
        stats[name+"Range"]=[min(values),max(values)] if values else None
    stats["failedOrIncompleteIds"]=[r["id"] for r in rows if len(r.get("grid",{}).get("acceptedLevels",[]))<2]
    return stats


def engineering_comparison(matrix):
    comparisons=[]
    candidates={r["id"]:r for r in matrix["rows"]}
    for baseline in matrix["baselines"]:
        candidate=candidates[baseline["id"].removesuffix("-engineering")]
        pair=dict(candidate=candidate["id"],baseline=baseline["id"],pairedLevels=[],
                  porousReduction=1-candidate["porousObjective"]/baseline["porousObjective"])
        bm={r["level"]:r.get("metrics") for r in baseline.get("native",[])}
        for a in candidate.get("native",[]):
            b=bm.get(a["level"])
            if a.get("metrics") is not None and b is not None:
                ja,jb=a["metrics"]["dimensionlessTotalPower"],b["dimensionlessTotalPower"]
                pair["pairedLevels"].append(dict(level=a["level"],candidate=ja,baseline=jb,trueReduction=1-ja/jb))
        comparisons.append(pair)
    return comparisons


def plotting():
    os.environ.setdefault("MPLCONFIGDIR",str(STUDY/"plot-cache"))
    os.environ.setdefault("XDG_CACHE_HOME",str(STUDY/"plot-cache"))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size":9,"axes.spines.top":False,"axes.spines.right":False})
    return plt


def matching_audit(matrix):
    flow=[];area=[];orders=[];gci=[]
    for row in matrix["rows"]+matrix["baselines"]:
        p=row["problem"];n=2 if p["case"] in ("double-pipe","four-terminal") else 1
        expected=n*(2/3)*.02*p["port_width"]
        for grid in row.get("native",[]):
            m=grid.get("metrics")
            if m is not None:
                flow.append(abs(m["inletFlux"]/expected-1))
                area.append(abs(m["area"]-p["volume_fraction"]*p["width"]*p["height"]))
        if row.get("grid",{}).get("method")=="three-grid-GCI":
            orders.append(row["grid"]["observedOrder"]);gci.append(row["grid"]["gciRelative"])
    return dict(maxRelativeInletFlowMismatch=max(flow,default=None),maxAbsoluteAreaMismatch=max(area,default=None),
                observedOrderRange=[min(orders),max(orders)] if orders else None,
                relativeGciRange=[min(gci),max(gci)] if gci else None)


def poiseuille_sanity(matrix):
    row=next(r for r in matrix["rows"] if r["id"]=="double-pipe-re-0-da-0.0001")
    p=row["problem"];path=ROOT/row["nativePath"]/"sharp/final.npz"
    with np.load(path) as field:
        groups=contours(field["rho"],float(field["width"]),float(field["height"]))
    valid=len(groups)==2;dimensions=[]
    for group in groups:
        valid &= len(group)==1
        xy=group[0];lo=xy.min(0);hi=xy.max(0);tolerance=1e-11+1e-9*max(hi-lo)
        on_boundary=np.minimum(abs(xy-lo),abs(xy-hi)).min(1)
        valid &= bool(np.all(on_boundary<=tolerance) and abs(abs(signed_area(xy))-np.prod(hi-lo))<=tolerance)
        dimensions.append((hi-lo).tolist())
    if not valid:return dict(status="not-rectangular-no-analytic-substitution")
    expected=sum((8*p["viscosity"]*length/width**2)*(2*width/3) for length,width in dimensions)
    return dict(id=row["id"],actualExtractedRectangles=True,dimensions=dimensions,analyticDimensionlessPower=expected,
        rows=[dict(level=x["level"],J_T=x["metrics"]["dimensionlessTotalPower"],relativeError=x["metrics"]["dimensionlessTotalPower"]/expected-1)
              for x in row["native"] if x.get("metrics")],
        scope="Independent plane Poiseuille sanity check of actual rectangles; no substitution for failed curved designs or general physical accuracy.")


def scatter(report,path):
    plt=plotting()
    fig,axes=plt.subplots(2,2,figsize=(12,10),constrained_layout=True)
    for ax,case in zip(axes.flat,CASES):
        bounds=[];notes=[]
        for re,color,marker in ((0,"#1766a6","o"),(50,"#b94632","s")):
            groups=[g for g in report["groups"] if g["case"]==case and g["reynolds"]==re]
            if not groups:continue
            g=groups[0];r=ranking(g["rows"])
            tau="n/a" if r["kendallTau"] is None else f"{r['kendallTau']:.3f}"
            notes.append(f"Re={re}: tau={tau}, available {r['available']}/{r['attempted']}, GCI {r['gciAvailable']}")
            labelled=False
            for row in g["rows"]:
                grid=row.get("grid",{});jt=grid.get("J_T")
                if jt is None:continue
                jb=row["porousObjective"];u=grid.get("uncertainty")
                label=f"Re={re}" if not labelled else None;labelled=True
                ax.errorbar(jb,jt,yerr=u,fmt=marker,color=color,mfc=color if u is not None else "white",
                            ms=5,capsize=3,lw=1,label=label)
                bounds.extend([jb,jt]+([jt-u,jt+u] if u is not None else []))
        if bounds:
            lo=min(0,min(bounds));hi=max(bounds)*1.08
            ax.plot([lo,hi],[lo,hi],"--",color="0.45",lw=1,label="y=x")
            ax.set(xlim=(lo,hi),ylim=(lo,hi))
        ax.set(title=case,xlabel="Porous total-pressure power, J_B",ylabel="Sharp-wall total-pressure power, J_T")
        ax.grid(alpha=.2);ax.legend(loc="upper left",fontsize=8)
        ax.text(.02,.02,"\n".join(notes),transform=ax.transAxes,fontsize=8,va="bottom",
                bbox=dict(facecolor="white",alpha=.85,edgecolor="none"))
    fig.suptitle("Porous predictions vs native Cut-cell results | all available candidates\nBars: conditional three-grid GCI; hollow: unavailable GCI; failed candidates retained in JSON",fontsize=12)
    fig.savefig(path,dpi=145);plt.close(fig)


def gap_plot(matrix,path):
    plt=plotting()
    fig,axes=plt.subplots(4,2,figsize=(12,14),constrained_layout=True)
    re_values=[0,10,50,100];darcy=[1e-3,1e-4,1e-5]
    for index,case in enumerate(CASES):
        rows=[r for r in matrix["rows"] if r["case"]==case]
        by={(r["reynolds"],r["darcyPort"]):r for r in rows}
        for col,series in enumerate((darcy,re_values)):
            ax=axes[index,col]
            for k,value in enumerate(series):
                data=[by.get((re,value) if col==0 else (value,da)) for re,da in
                      ([(re,0) for re in re_values] if col==0 else [(0,da) for da in darcy])]
                x=[r["reynolds"] if col==0 else r["alphaMax"] for r in data]
                y=[r.get("gap",np.nan) for r in data]
                line=ax.plot(x,y,"-",lw=1,label=f"Da_w={value:g}" if col==0 else f"Re={value}")[0]
                for xx,yy,row in zip(x,y,data):
                    if not np.isfinite(yy):continue
                    interval=row.get("gapInterval")
                    error=np.array([[yy-interval[0]],[interval[1]-yy]]) if interval else None
                    ax.errorbar(xx,yy,yerr=error,fmt="o",ms=4,capsize=2,color=line.get_color(),
                                mfc=line.get_color() if interval else "white")
            ax.axhline(0,color="0.45",ls="--",lw=1);ax.grid(alpha=.2)
            ax.set(title=case,xlabel="Mean-port Reynolds number" if col==0 else "Maximum Brinkman drag, alpha_max",
                   ylabel="gap = (J_B - J_T) / J_T")
            if col==1:ax.set_xscale("log")
            ax.legend(fontsize=8,ncol=2)
            missing=sum("gap" not in r for r in rows)
            single=sum(len(r.get("grid",{}).get("acceptedLevels",[]))==1 for r in rows)
            ax.text(.02,.02,f"{missing}/12 missing J_T; {single}/12 single-grid only",transform=ax.transAxes,fontsize=8,
                    bbox=dict(facecolor="white",alpha=.8,edgecolor="none"))
    fig.suptitle("Fidelity gap across Re and resistance | fixed area per case\nBars: propagated conditional GCI only; hollow: uncertainty unavailable. MAC discretization + extraction also contribute.",fontsize=12)
    fig.savefig(path,dpi=140);plt.close(fig)


def report(matrix_path,ranking_path,artifact):
    a=load(STUDY/"a/gradients/summary.json");b=load(STUDY/"b-complete.json")
    extra=STUDY/"a/total-pressure-gradients/summary.json"
    if extra.exists():
        a["curves"]+=load(extra)["curves"]
    plotting()
    import sys
    sys.path.insert(0, str(ROOT/'tools/visualization'))
    from plot_topology_gradients import plot_report
    plot_report(a,artifact.with_name("native-topology-fidelity-gradients.png"))
    c=load(matrix_path);d=load(ranking_path)
    scatter(d,artifact.with_name("native-topology-fidelity-ranking.png"))
    gap_plot(c,artifact.with_name("native-topology-fidelity-gap.png"))
    result=dict(schema="cartmesh2d-topology-fidelity-v1",physicalAccuracyQualified=False,
        stageA=dict(curves=len(a["curves"]),passed=sum(x["passed"] for x in a["curves"]),
                    interiorMinima=sum(x["interiorMinimum"] for x in a["curves"]),
                    worstMinimumError=max(x["minimumError"] for x in a["curves"]),
                    maximumEquationResidual=max(x["residual"] for x in a["curves"]),
                    exactLegacyComparison=load(STUDY/"a/re0-legacy.json")),
        stageB=dict(cases=CASES,rows=[compact(r) for r in b["rows"]],
                    optimizationStatuses=dict(Counter(r["status"] for r in b["rows"])),
                    initialAndRetryAttempts=len(b["attempts"]),limits=b["limits"]),
        stageC=dict(status=c["status"],controls=c["controls"],recoveryRule=c.get("recoveryRule"),statistics=counts(c["rows"]),
                    rows=[compact(r) for r in c["rows"]],baselines=[compact(r) for r in c["baselines"]],
                    baselineStatistics=counts(c["baselines"]),engineeringComparisons=engineering_comparison(c)),
        stageD=dict(status=d["status"],rules=d["rules"],statistics=counts([r for g in d["groups"] for r in g["rows"]]),
                    groups=[dict(case=g["case"],reynolds=g["reynolds"],
                    rows=[compact(r) for r in g["rows"]],distinctness=g.get("distinctness"),ranking=ranking(g["rows"])) for g in d["groups"]]),
        definitions=dict(objective="Port flux-weighted total pressure power. Native kinematic power/(nu*Upeak^2/mu), per density and depth; Re=0 uses exact Stokes (zero kinetic/inertia coefficient).",
            darcy="Da_w=mu/(alpha_max*w^2). C uses 1e-3,1e-4,1e-5. Da_H=Da_w*(w/H)^2.",
            uncertainty="GCI=1.25*abs(J_f-J_m)/(r_fm^p-1), h=sqrt(area/cells); p from three-grid nonuniform-ratio Richardson equation. Conditional on asymptotic power law; two-grid change is not GCI.",
            reversal="Strict opposite signs of J_B difference and J_T difference, resolved only if abs(J_T difference)>GCI_a+GCI_b.",
            leakage="Absolute internal-face throughput fraction through mean rho<0.5, counting repeated crossings; not a distinct inlet-particle fraction."),
        limitations=["Benchmark boundary adaptations prevent exact literature objective/crossover reproduction.",
            "Fixed MAC analysis grids have no spatial convergence qualification; gap includes their error and equal-area sharp extraction.",
            "Budget-exhausted optimization, extraction, mesh-quality and native convergence failures remain visible.",
            "GCI is a conditional discretization estimate; missing/oscillatory grids do not receive invented uncertainty.",
            "Grid-resolved reversal is evidence within these equations and discretizations, not a physically qualified engineering ranking.",
            "No external checkMesh, full native CTest, frontend/App or platform CI qualification in this bounded CLI study."],
        sources=[evidence(p) for p in (STUDY/"a/gradients/summary.json",STUDY/"a/re0-legacy.json",STUDY/"b-complete.json",matrix_path,ranking_path)],
        resources=dict(freeBytes=shutil.disk_usage(ROOT).free,studyBytes=sum(p.stat().st_size for p in STUDY.rglob("*") if p.is_file())))
    result["figures"]=[evidence(p) for p in sorted(artifact.parent.glob("native-topology-fidelity-*.png"))]
    result["workingTreeProvenance"]=dict(
        head=subprocess.check_output(["git","rev-parse","HEAD"],cwd=ROOT,text=True).strip(),
        branch=subprocess.check_output(["git","branch","--show-current"],cwd=ROOT,text=True).strip(),
        newCommitHashes=[],commitIssue="Sandbox refused .git/index.lock; planned per-stage paths/messages in outputs/codex_tasks/T01_progress.md. No push.",
        sourceSnapshotRole="Final working tree. Per-run source/executable hashes and executed commands remain in raw ledgers; this snapshot does not replace their historical provenance.",
        files=[evidence(ROOT/p) for p in SOURCE_PATHS]+[evidence(p) for p in sorted((ROOT/"tools/optimization").glob("*.py"))],
        executables=[evidence(ROOT/"build"/p) for p in ("cartmesh2d_cli","cartmesh2d_flow_cli")])
    if extra.exists():result["sources"].append(evidence(extra))
    for ledger in (c.get("initialLedger"),d.get("initialRankingLedger")):
        if ledger:result["sources"].append(evidence(ROOT/ledger))
    result["stageC"]["matchingAudit"]=matching_audit(c)
    result["stageC"]["poiseuilleSanity"]=poiseuille_sanity(c)
    for group in result["stageD"]["groups"]:
        selected=set((group.get("distinctness") or {}).get("representativeIds",[]))
        group["rankingDistinctMaterialRepresentatives"]=ranking([r for r in group["rows"] if r["id"] in selected])
        group["duplicateRule"]="First predeclared representative within density Linf 1e-6; full all-candidate ranking is also retained. This is not a claim of distinct sharp-wall connectivity."
    write_json(artifact,result)
    if artifact.stat().st_size>=2_000_000 or any(f["bytes"]>=2_000_000 for f in result["figures"]):
        raise RuntimeError("artifact exceeds the task's 2 MB per-file ceiling")
    print(json.dumps(dict(stageC=result["stageC"]["statistics"],stageD=[dict(case=g["case"],reynolds=g["reynolds"],
          distinct=g["distinctness"]["distinctMaterialFields"],**{k:g["ranking"][k] for k in ("attempted","available","gciAvailable","kendallTau")},
          resolvedReversals=len(g["ranking"]["resolvedReversals"])) for g in result["stageD"]["groups"]]),indent=2))
    return result


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--matrix",type=Path,default=STUDY/"c-recovery/summary.json")
    p.add_argument("--ranking",type=Path,default=STUDY/"d-uncertainty/summary.json")
    p.add_argument("--artifact",type=Path,default=ROOT/"artifacts/current/native-topology-fidelity.json")
    args=p.parse_args();report(args.matrix.resolve(),args.ranking.resolve(),args.artifact.resolve())
