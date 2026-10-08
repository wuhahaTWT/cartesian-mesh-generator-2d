#!/usr/bin/env python3
"""Bounded, resumable T01 experiment driver; all attempts remain in the ledger."""
import argparse
from contextlib import redirect_stdout
from dataclasses import asdict
import hashlib
import json
import os
from pathlib import Path
import time
import copy

import numpy as np

from brinkman import Problem
from navier_stokes_brinkman import NavierStokesBrinkman
from topology_cases import CASES,problem_parameters
from engineering_baselines import best_reference
from optimize_flow import parser as optimizer_parser,run as optimize,snapshot,write_json,metrics
from topology_artifacts import extract,contours,port_connectivity
from compare_sharp_designs import prepare_design,flow_metrics
import native_flow as bridge
from fidelity_metrics import leakage,grid_uncertainty,ranking

ROOT=Path(__file__).resolve().parents[2]
STUDY=ROOT/"outputs/topology-fidelity"


def analyze_case(output,case,reynolds,da=1e-4,filter_scale=1.,initialization="uniform",iterations=(100,180,800),seed=0):
    output=output.resolve()

    spec=problem_parameters(case)
    alpha=1/(da*spec["port_width"]**2)
    args=optimizer_parser().parse_args(["--output",str(output),"--case",case,
        "--nx",str(spec["nx"]),"--ny",str(spec["ny"]),"--width",str(spec["width"]),
        "--port-width",str(spec["port_width"]),"--volume",str(spec["volume_fraction"]),
        "--filter-radius",str(.055*filter_scale),"--alpha-max",str(alpha),"--reynolds",str(reynolds),
        "--objective","total-pressure-power","--initialization",initialization,"--seed",str(seed),
        "--iterations",*map(str,iterations),"--max-seconds","300","--snapshot-every","1000","--no-plot"])
    output.parent.mkdir(parents=True,exist_ok=True)
    start=time.monotonic()
    row=dict(case=case,reynolds=reynolds,darcyPort=da,alphaMax=alpha,filterScale=filter_scale,
             initialization=initialization,seed=seed,path=str(output.relative_to(ROOT)),status="running",native=[])
    with output.with_suffix(".log").open("w") as log,redirect_stdout(log):
        try:
            report=optimize(args)
            row.update(status="budget-exhausted" if report["status"] in ("iteration-limit","time-limit") else report["status"],
                       optimizerStatus=report["status"],optimizationConverged=report["optimizationConverged"],
                       final=report["final"],finalParameters=report["finalParameters"],problem=report["problem"])
            model=NavierStokesBrinkman(Problem(**report["problem"]))
            beta,q=report["finalParameters"]["beta"],report["finalParameters"]["q"]
            with np.load(output/"final.npz") as f:
                x=f["design"].ravel()
            evaluation=model.evaluate(x,q,beta,args.objective)
            row["porousObjective"]=evaluation.objective
            row["leakage"]=leakage(model,evaluation,q)
            row["extraction"]=extract(output)
            best,baselines=best_reference(model,beta,q,args.objective)
            label,bx,be=best
            snapshot(output/"reference-engineering.npz",model,bx,be,q,beta)
            row["engineeringReferences"]=baselines
            row["engineeringBest"]=dict(label=label,**metrics(model,bx,be))
            row["porousReductionVsEngineering"]=1-evaluation.objective/be.objective
            row["samePorousArea"]=dict(candidate=evaluation.volume*model.cells*model.area,
                                      baseline=be.volume*model.cells*model.area)
        except (ValueError,ArithmeticError,RuntimeError) as exc:
            row.update(status="experiment-failed",issue=str(exc))
    row["elapsedSeconds"]=time.monotonic()-start
    write_json(output.with_suffix(".json"),row)
    return row


def benchmarks(output):
    output.mkdir(parents=True,exist_ok=False)
    report=dict(schema="cartmesh2d-topology-benchmarks-v1",cases=CASES,rows=[],physicalAccuracyQualified=False,
        limits=["Diffuser is an adapted case: exact original full-text parameters unavailable.",
                "Four-terminal omits literature leads and prescribes outlet velocity; no exact crossover/objective reproduction.",
                "MAC analysis grid is a research starting resolution, not mesh-independent.",
                "Stationarity uses the existing projectedKkt <=1e-3; iteration budgets are not convergence."])
    write_json(output/"summary.json",report)
    for case in CASES:
        for re in (0,50):
            row=analyze_case(output/f"{case}-re-{re}",case,re)
            report["rows"].append(row)
            write_json(output/"summary.json",report)
            print(case,re,row["status"],row.get("final",{}).get("projectedKkt"),flush=True)

    write_json(output/"summary.json",report)
    plot_benchmarks(report,output/"benchmarks.png")
    return report


def native_design(row,output,levels=(4,5,6),source_name="final.npz",iterations=6000,timeout=90):
    """Retain each level independently, including meshing and flow failures."""
    if "problem" not in row:
        row["nativeStatus"]="no-accepted-analysis";return row
    source=ROOT/row["path"]
    model=NavierStokesBrinkman(Problem(**row["problem"]))
    target=model.problem.volume_fraction*model.problem.width*model.problem.height
    output.mkdir(parents=True,exist_ok=False)
    row["nativePath"]=str(output.relative_to(ROOT));row["native"]=[]
    try:
        design=prepare_design(output/"sharp",source/source_name,model,target,row.get("id","candidate"))
        row["sharpProjection"]=design
    except (ArithmeticError,ValueError,RuntimeError) as exc:
        row.update(nativeStatus="extraction-failed",nativeIssue=str(exc))
        write_json(output/"summary.json",row);return row
    for level in levels:

        item=dict(level=level,metrics=None,path=str((output/f"level-{level}").relative_to(ROOT)),
                  meshPaddingFraction=1/30,meshAggregationFraction=.25)
        row["native"].append(item)
        args=argparse.Namespace(directory=output/"sharp",output=output/f"level-{level}",
            mesh_cli=ROOT/"build/cartmesh2d_cli",flow_cli=ROOT/"build/cartmesh2d_flow_cli",
            levels=[level],timeout=timeout,iterations=iterations,speed=.02,nu=1.,tolerance=1e-8,
            small_alpha=.25,mesh_only=False,component=None,reynolds=row["reynolds"],
            convergence="strict",velocity_relaxation=.2,steady_acceleration="anderson")
        with (output/f"level-{level}.log").open("w") as log,redirect_stdout(log):
            try:
                report=bridge.run(args)
                item["metrics"]=flow_metrics(report)
                item["status"]="accepted" if item["metrics"] is not None else "rejected"
                leaves=report.get("components",[report])
                item["attemptStatuses"]=[c["status"] for leaf in leaves for c in leaf["cases"]]
            except (ArithmeticError,ValueError,RuntimeError) as exc:
                item.update(status="bridge-failed",issue=str(exc))
        write_json(output/"summary.json",row)
    row["grid"]=grid_uncertainty(row["native"])
    jt=row["grid"].get("J_T")
    if jt is not None and jt>0:
        row["gap"]=(row["porousObjective"]-jt)/jt
        uncertainty=row["grid"].get("uncertainty")
        if uncertainty is not None and jt>uncertainty:
            row["gapInterval"]=[row["porousObjective"]/(jt+uncertainty)-1,row["porousObjective"]/(jt-uncertainty)-1]
    row["nativeStatus"]="evaluated"
    write_json(output/"summary.json",row)
    return row


def matrix_study(output,benchmarks_path):
    output.mkdir(parents=True,exist_ok=False)
    b=json.loads(benchmarks_path.read_text())
    reusable={(r["case"],r["reynolds"]):r for r in b["rows"]}
    report=dict(schema="cartmesh2d-topology-fidelity-matrix-v1",rows=[],baselines=[],
                controls=dict(reynolds=[0,10,50,100],darcyPort=[1e-3,1e-4,1e-5],levels=[4,5,6],
                    tolerance=1e-8,iterations=6000,secondsPerNativeAttempt=90,speed=.02,
                    areaProjection="Interior density offset for equal sharp area; no smoothing, hole filling or component deletion"),
                status="running",physicalAccuracyQualified=False)
    write_json(output/"summary.json",report)
    try:
        for case in CASES:
            for re in (0,10,50,100):
                for da in (1e-3,1e-4,1e-5):
                    identity=f"{case}-re-{re}-da-{da:g}"
                    if da==1e-4 and (case,re) in reusable:
                        row=copy.deepcopy(reusable[case,re]);row["analysisReusedFromB"]=True
                        model=NavierStokesBrinkman(Problem(**row["problem"]))
                        with np.load(ROOT/row["path"]/"final.npz") as f:
                            x=f["design"].ravel();q=float(f["q"]);beta=float(f["beta"])
                        e=model.evaluate(x,q,beta,"total-pressure-power")
                        row["leakage"]=leakage(model,e,q)
                    else:
                        row=analyze_case(output/"analysis"/identity,case,re,da)
                    row["id"]=identity
                    report["rows"].append(row)
                    write_json(output/"summary.json",report)
                    native_design(row,output/"native"/identity)
                    write_json(output/"summary.json",report)
                    print(identity,row["status"],row.get("grid",{}).get("method"),row.get("gap"),flush=True)
                    if da==1e-4 and "engineeringBest" in row:
                        baseline=copy.deepcopy(row)
                        baseline.update(id=identity+"-engineering",porousObjective=row["engineeringBest"]["objective"],
                                        status="engineering-reference",native=[])
                        for key in ("grid","gap","gapInterval","leakage","nativeStatus","nativePath","sharpProjection"):
                            baseline.pop(key,None)
                        report["baselines"].append(baseline)
                        native_design(baseline,output/"native"/(identity+"-engineering"),source_name="reference-engineering.npz")
                        write_json(output/"summary.json",report)
        report["status"]="completed-attempts"
    except KeyboardInterrupt as exc:
        report.update(status="interrupted",issue=str(exc))

    write_json(output/"summary.json",report)
    return report


def refresh_grid(row):
    for key in ("gap","gapInterval"):
        row.pop(key,None)
    row["grid"]=grid_uncertainty(row.get("native",[]))
    jt=row["grid"].get("J_T")
    if jt is not None and jt>0:
        jb=row["porousObjective"]
        row["gap"]=(jb-jt)/jt
        uncertainty=row["grid"].get("uncertainty")
        if uncertainty is not None and jt>uncertainty:
            row["gapInterval"]=[jb/(jt+uncertainty)-1,jb/(jt-uncertainty)-1]


def retry_native(row,output,levels,source_name="final.npz"):
    """Never erase an attempt or choose a grid using its objective value."""
    attempt=copy.deepcopy(row)
    for key in ("nativeAttempts","grid","gap","gapInterval"):
        attempt.pop(key,None)
    native_design(attempt,output,levels=levels,source_name=source_name)
    row.setdefault("nativeAttempts",copy.deepcopy(row.get("native",[])))
    row["nativeAttempts"].extend(copy.deepcopy(attempt.get("native",[])))
    selected={r["level"]:r for r in row.get("native",[])}
    for item in attempt.get("native",[]):
        previous=selected.get(item["level"])
        if previous is None or (previous.get("metrics") is None and item.get("metrics") is not None):
            selected[item["level"]]=item
    row["native"]=sorted(selected.values(),key=lambda r:r["level"])
    row.setdefault("nativeRecovery",[]).append(dict(path=str(output.relative_to(ROOT)),levels=list(levels),
        status=attempt.get("nativeStatus"),issue=attempt.get("nativeIssue")))
    refresh_grid(row)


def recover_matrix(output,source):
    output.mkdir(parents=True,exist_ok=False)
    report=copy.deepcopy(json.loads(source.read_text()))
    report.update(status="recovering",initialLedger=str(source.relative_to(ROOT)),
        recoveryRule="Retry all boundary-template failures with an actual-geometry template; if fewer than two grids accepted, try level 7 then level 3. Keep every rejection; first accepted attempt per level; finest accepted J_T. No objective-based selection.")
    write_json(output/"summary.json",report)
    try:
        for section in ("rows","baselines"):
            for row in report[section]:
                source_name="reference-engineering.npz" if section=="baselines" else "final.npz"
                levels=[r["level"] for r in row.get("native",[]) if "boundary-rejected" in r.get("attemptStatuses",[])]
                if levels:
                    retry_native(row,output/(row["id"]+"-boundary"),levels,source_name)
                    write_json(output/"summary.json",report)
                if row.get("nativeStatus")=="extraction-failed" or "problem" not in row:
                    continue
                for level in (7,3):
                    if len(row.get("grid",{}).get("acceptedLevels",[]))>=2:
                        break
                    retry_native(row,output/(row["id"]+f"-level-{level}"),[level],source_name)
                    write_json(output/"summary.json",report)
                print(row["id"],row.get("grid",{}).get("acceptedLevels"),row.get("grid",{}).get("method"),flush=True)
        report["status"]="completed-attempts"
    except KeyboardInterrupt as exc:
        report.update(status="interrupted",issue=str(exc))

    write_json(output/"summary.json",report)
    return report


def distinctness(rows):
    """A diagnostic of different material fields, not a geometry quality gate.

    Absolute density 1e-6 removes last-bit differences when counting duplicate
    optimized designs; it is recorded explicitly and is not a CFD tolerance.
    Full fields and every candidate remain available regardless of this count.
    """
    fields=[];representatives=[];pairs=[]
    for row in rows:
        path=ROOT/row["path"]/"final.npz"
        if not path.exists():
            continue
        with np.load(path) as f:
            rho=f["rho"].copy()
        row["densitySha256"]=hashlib.sha256(rho.tobytes()).hexdigest()
        fields.append((row["id"],rho))
        if not any(np.max(abs(rho-other))<=1e-6 for _,other in representatives):
            representatives.append((row["id"],rho))
    for i,(a,x) in enumerate(fields):
        for b,y in fields[i+1:]:
            pairs.append(dict(a=a,b=b,densityLinf=float(np.max(abs(x-y)))))
    return dict(materialFields=len(fields),distinctMaterialFields=len(representatives),densityAbsoluteTolerance=1e-6,
                representativeIds=[x[0] for x in representatives],allPairDistances=pairs,
                note="Material-field diagnostic only; distinct density fields need not have distinct connectivity or manufacturable wall separation.")


def rank_study(output,source):
    output.mkdir(parents=True,exist_ok=False)
    matrix=json.loads(source.read_text())
    variants=[dict(darcy=1e-3,filter=0.75,initialization="random",seed=1101),
              dict(darcy=1e-3,filter=1.25,initialization="random",seed=1102),
              dict(darcy=1e-4,filter=1.5,initialization="geometric",seed=0),
              dict(darcy=1e-3,filter=1.75,initialization="uniform",seed=0),
              dict(darcy=1e-3,filter=2.25,initialization="random",seed=1104)]
    extras=[dict(darcy=3e-3,filter=1.5,initialization="random",seed=1105),
            dict(darcy=3e-4,filter=2.,initialization="uniform",seed=0),
            dict(darcy=3e-3,filter=2.,initialization="random",seed=1107)]
    report=dict(schema="cartmesh2d-topology-ranking-v1",status="running",sourceMatrix=str(source.relative_to(ROOT)),
        groups=[],variants=variants,extraVariantsForDuplicateDesigns=extras,physicalAccuracyQualified=False,
        rules=["Three Darcy candidates reused from C plus five predeclared starts/filter radii per group.",
               "Extra candidates only if fewer than eight distinct density fields; no selection by objective or reversal.",
               "Native levels 4/5/6; if fewer than two accepted, try 3 then 7 to bound disk use. Every failed attempt retained.",
               "Resolved reversal requires abs(J_Ta-J_Tb) > GCI_a+GCI_b; two-grid changes never replace GCI.",
               "Tau includes all available values; separate >=2-grid and GCI-only subsets expose missing uncertainty."])
    write_json(output/"summary.json",report)
    try:
        for case in CASES:
            for re in (0,50):
                group=dict(case=case,reynolds=re,rows=[copy.deepcopy(r) for r in matrix["rows"] if r["case"]==case and r["reynolds"]==re])
                for row in group["rows"]:row["analysisReusedFromC"]=True
                report["groups"].append(group)
                for i,v in enumerate(variants+extras):
                    group["distinctness"]=distinctness(group["rows"])
                    if i>=len(variants) and group["distinctness"]["distinctMaterialFields"]>=8:
                        break
                    initialization="merged" if case=="double-pipe" and v["initialization"]=="geometric" else v["initialization"]
                    identity=f"{case}-re-{re}-variant-{i+1}"
                    row=analyze_case(output/"analysis"/identity,case,re,v["darcy"],v["filter"],initialization,seed=v["seed"])
                    row["id"]=identity;group["rows"].append(row)
                    write_json(output/"summary.json",report)
                    native_design(row,output/"native"/identity)
                    write_json(output/"summary.json",report)
                    if row.get("nativeStatus")!="extraction-failed" and "problem" in row:
                        for level in (3,7):
                            if len(row.get("grid",{}).get("acceptedLevels",[]))>=2:break
                            retry_native(row,output/"native"/(identity+f"-level-{level}"),[level])
                            write_json(output/"summary.json",report)
                    group["distinctness"]=distinctness(group["rows"])
                    group["ranking"]=ranking(group["rows"])
                    write_json(output/"summary.json",report)
                    print(identity,row["status"],row.get("grid",{}).get("acceptedLevels"),
                          "distinct",group["distinctness"]["distinctMaterialFields"],flush=True)
                group["distinctness"]=distinctness(group["rows"])
                group["ranking"]=ranking(group["rows"])
                write_json(output/"summary.json",report)
        report["status"]="completed-attempts"
    except KeyboardInterrupt as exc:
        report.update(status="interrupted",issue=str(exc))

    write_json(output/"summary.json",report)
    return report


def rank_uncertainty(output,source):
    """Budgeted third-grid attempt for EVERY eligible two-grid candidate."""
    output.mkdir(parents=True,exist_ok=False)
    report=copy.deepcopy(json.loads(source.read_text()))
    report.update(status="running",initialRankingLedger=str(source.relative_to(ROOT)))
    report["rules"].append("Uniform uncertainty supplement: for exactly two accepted grids and no previous level-3 attempt, try level 3 once; independent of objective, order or reversal. All original attempts remain.")
    write_json(output/"summary.json",report)
    try:
        for group in report["groups"]:
            for row in group["rows"]:
                attempted={r["level"] for r in row.get("nativeAttempts",row.get("native",[]))}
                if len(row.get("grid",{}).get("acceptedLevels",[]))==2 and 3 not in attempted:
                    retry_native(row,output/row["id"],[3])
                    print(row["id"],row["grid"]["acceptedLevels"],row["grid"]["method"],flush=True)
                group["ranking"]=ranking(group["rows"])
                write_json(output/"summary.json",report)
        report["status"]="completed-attempts"
    except KeyboardInterrupt as exc:
        report.update(status="interrupted",issue=str(exc))

    write_json(output/"summary.json",report)
    return report


def plot_benchmarks(report,path):
    os.environ.setdefault("MPLCONFIGDIR",str(STUDY/"plot-cache"))
    os.environ.setdefault("XDG_CACHE_HOME",str(STUDY/"plot-cache"))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(4,4,figsize=(14,12),constrained_layout=True)
    for row_index,case in enumerate(CASES):
        # Schematic of reported connectivity, not copied/digitized geometry.
        reference=axes[row_index,0]
        if case=="double-pipe":
            for y in (.25,.75):
                reference.plot([0,.35,.65,1],[y,.5,.5,y],color="0.35")
            note="Borrvall-Petersson: shared trunk\nSource geometry not digitized"
        elif case=="elbow":
            reference.plot([0,.8],[.8,0],color="0.35",label="low Re shortcut")
            angle=np.linspace(np.pi/2,0,65)
            reference.plot(.1+.7*np.cos(angle),.1+.7*np.sin(angle),"--",color="#a65628",label="inertia rounds turn")
            reference.legend(fontsize=7)
            note="B-P / Gersborg-Hansen\nDifferent fluid drag and Re scales"
        elif case=="four-terminal":
            angle=np.linspace(np.pi/2,-np.pi/2,65)
            for x in (.3*np.cos(angle),1-.3*np.cos(angle)):
                reference.plot(x,.5+.2*np.sin(angle),color="0.35")
            for y in (.3,.7):
                reference.plot([0,1],[y,y],"--",color="#a65628")
            note="Olesen: U-returns (low Re)\nParallel paths (high Re); leads omitted here"
        else:
            reference.text(.5,.5,"Original numerical geometry\nnot verified",ha="center",va="center",fontsize=9)
            note="B-P diffuser-inspired adaptation\nNo literature-shape equivalence claim"
        reference.set(xlim=(0,1),ylim=(0,1),aspect="equal")
        reference.set_title("Literature topology (schematic)\n"+note,fontsize=8)
        rows=[r for r in report["rows"] if r["case"]==case]
        for col,re in enumerate((0,50)):
            r=next(r for r in rows if r["reynolds"]==re)
            ax=axes[row_index,col+1]
            if "porousObjective" not in r:
                ax.text(.5,.5,r.get("issue",r["status"]),ha="center",wrap=True,fontsize=7)
                continue
            with np.load(ROOT/r["path"]/"final.npz") as f:
                groups=contours(f["rho"],float(f["width"]),float(f["height"]))
            for group in groups:
                for loop in group:
                    closed=np.vstack((loop,loop[0]));ax.plot(closed[:,0],closed[:,1],lw=.8)
            ax.set(title=f"{case}, Re={re}\n{r['status']}; KKT={r['final']['projectedKkt']:.2g}",
                   xlim=(0,r["problem"]["width"]),ylim=(0,1),aspect="equal")
        r=rows[0];ax=axes[row_index,3]
        if "engineeringBest" in r:
            with np.load(ROOT/r["path"]/"reference-engineering.npz") as f:
                groups=contours(f["rho"],float(f["width"]),float(f["height"]))
            for group in groups:
                for loop in group:
                    closed=np.vstack((loop,loop[0]));ax.plot(closed[:,0],closed[:,1],color="#a65628",lw=.8)
            ax.set(title=f"Engineering reference (Re=0)\n{r['engineeringBest']['label']}",
                   xlim=(0,r["problem"]["width"]),ylim=(0,1),aspect="equal")
    fig.suptitle("Actual unsmoothed contours | 4 literature-inspired cases\nAdapted boundary conditions; no exact literature reproduction or physical qualification",fontsize=12)
    for ax in axes[:,1:].flat:
        ax.title.set_fontsize(9)
    fig.savefig(path,dpi=135);plt.close(fig)


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("stage",choices=["b","c","c-recovery","d","d-uncertainty"])
    p.add_argument("--output",type=Path,required=True)
    p.add_argument("--source",type=Path,help="Previous-stage ledger; omitted uses this task's canonical preserved evidence")
    args=p.parse_args()
    if args.stage=="b":benchmarks(args.output.resolve())
    elif args.stage=="c":matrix_study(args.output.resolve(),args.source.resolve() if args.source else STUDY/"b-complete.json")
    elif args.stage=="c-recovery":recover_matrix(args.output.resolve(),args.source.resolve() if args.source else STUDY/"c-matrix/summary.json")
    elif args.stage=="d":rank_study(args.output.resolve(),args.source.resolve() if args.source else STUDY/"c-recovery/summary.json")
    else:rank_uncertainty(args.output.resolve(),args.source.resolve() if args.source else STUDY/"d-ranking/summary.json")
