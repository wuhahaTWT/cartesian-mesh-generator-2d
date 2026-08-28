#!/usr/bin/env python3
"""Reproduce Q2 evidence; never confuse partial superellipse repair with full Q2."""
from __future__ import annotations
import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import subprocess
import sys

from generate_q1_baselines import CASES
from generate_q0_baselines import read_cm2d

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def run(command, log=None):
    result = subprocess.run([str(x) for x in command], text=True,
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=600)
    if log:
        log.write_text(result.stdout, encoding="utf-8")
    if result.returncode:
        raise RuntimeError(f"{command}: {result.stdout}")
    return result.stdout

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--executable", type=Path, default=Path("build-q2/cartmesh2d_hybrid_cli"))
    parser.add_argument("--before", type=Path, default=Path("build-q2/before"))
    parser.add_argument("--after", type=Path, default=Path("build-q2/after"))
    parser.add_argument("--output", type=Path, default=Path("artifacts/q2/comparison.json"))
    parser.add_argument("--collect-only", action="store_true")
    parser.add_argument("--openfoam", action="store_true")
    parser.add_argument("--repeat", action="store_true")
    parser.add_argument("--require-full-acceptance", action="store_true")
    args=parser.parse_args()
    args.after.mkdir(parents=True, exist_ok=True)
    repo=Path.cwd()
    result={"format_version":"cartmesh2d-q2-comparison-v1",
            "base_commit":run(["git","rev-parse","556bb90"]).strip(),
            "source_diff_sha256":hashlib.sha256(run(["git","diff","HEAD"]).encode()).hexdigest(),
            "image":"opencfd/openfoam-run:2606", "cases":{}}
    source_files=[Path("CMakeLists.txt")]
    for directory in ("src","include","apps"):
        source_files.extend(p for p in Path(directory).rglob("*") if p.is_file())
    result["source_files_sha256"]={str(p):sha(p) for p in sorted(source_files)}
    short="face_length_over_local_background_h"
    for name,(boundary,params) in CASES.items():
        prefix=args.after/name
        case=args.after/f"{name}-case"
        command=[args.executable,boundary,prefix,*map(str,params),case,"0.01"]
        if not args.collect_only:
            stdout=run(command,args.after/f"{name}.generation.log")
            if "hybrid_status=success" not in stdout:
                raise AssertionError(f"{name}: fallback is not Q2 acceptance: {stdout}")
        post=json.loads(Path(f"{prefix}.hybrid.quality-contract.json").read_text())
        pre=json.loads((args.before/f"{name}.hybrid.quality-contract.json").read_text())
        geometry=json.loads(Path(f"{prefix}.hybrid.json").read_text())
        registry=json.loads(Path(f"{prefix}.hybrid.intersections.json").read_text())
        # The H4-2 reader requires closed two-valent layer fronts and rejects
        # the unchanged Q1 H4-3 fronts, whose ends terminate on the wall.
        # Use the full independent OpenFOAM owner/neighbour/volume reader for
        # all five; additionally use the specialized H4-2 reader where valid.
        if name in ("circle","superellipse"):
            run([sys.executable,"tools/verification/check_hybrid_mesh2d.py",
                 f"{prefix}.hybrid.vtk",f"{prefix}.hybrid.json","--output",
                 args.after/f"{name}.independent-hybrid.json"])
        run([sys.executable,"tools/verification/check_openfoam2d.py",case,
             "--report",args.after/f"{name}.independent-foam.json"])
        solver=read_cm2d(Path(f"{prefix}.hybrid.solver.cm2d"))
        for i,record in enumerate(registry["records"]):
            assert record["id"]==i and record["local_h"]>0
            delta=sum((a-b)**2 for a,b in zip(record["original_position"],
                record["canonical_vertex"]["position"]))**.5
            assert abs(delta-record["displacement"])<=record["local_h"]*1.e-12
            if record["solver_vertex_id"] is not None:
                vertex=solver.vertices[record["solver_vertex_id"]]
                assert sum((a-b)**2 for a,b in zip(vertex,record["canonical_vertex"]["position"]))**.5<=record["local_h"]*1.e-9
        deterministic=None
        if args.repeat:
            repeat=args.after/f"repeat-{name}"
            run([args.executable,boundary,repeat,*map(str,params)])
            suffixes=[".hybrid.cm2d",".hybrid.solver.cm2d",".hybrid.quality-contract.json",".hybrid.intersections.json"]
            deterministic=all(sha(Path(f"{prefix}{s}"))==sha(Path(f"{repeat}{s}")) for s in suffixes)
            assert deterministic, f"{name}: repeated geometry/report changed"
        foam_log=args.after/f"{name}.checkMesh.log"
        if args.openfoam:
            run(["docker","run","--rm","-v",f"{repo}:/home/openfoam/workingDir",
                 "-w","/home/openfoam/workingDir",result["image"],"checkMesh","-case",
                 "/home/openfoam/workingDir/"+str(case),"-writeAllFields"],foam_log)
        foam_status="NOT RUN"
        if foam_log.exists():
            log=foam_log.read_text()
            assert "Mesh OK" in log and "Failed " not in log, f"{name}: checkMesh failed"
            foam_status="PASS"
        def hard(report):
            return dict(sorted(Counter(i["metric"] for i in report["issues"] if i["level"]=="hard").items()))
        before_hard,after_hard=hard(pre),hard(post)
        metric_regressions={}
        for section in ("ordinary_metrics","boundary_layer_metrics"):
            for metric,summary in pre[section].items():
                old,new=summary["worst"],post[section][metric]["worst"]
                delta=new-old
                lower=summary["worst_direction"]=="min"
                if (delta<0 if lower else delta>0) and abs(delta)>1.e-8*max(1.,abs(old)):
                    metric_regressions[metric]={"before":old,"after":new,"worst_direction":summary["worst_direction"]}
        entry={"input":boundary,"input_sha256":sha(Path(boundary)),
            "generation_command":list(map(str,command)),
            "min_face_over_local_h_before":pre["ordinary_metrics"][short]["worst"],
            "min_face_over_local_h_after":post["ordinary_metrics"][short]["worst"],
            "min_absolute_face_before":pre["legacy_hard_safety"]["min_face_length_absolute"],
            "min_absolute_face_after":post["legacy_hard_safety"]["min_face_length_absolute"],
            "hard_issue_counts_before":before_hard,"hard_issue_counts_after":after_hard,
            "new_or_increased_hard_counts":{k:v for k,v in after_hard.items() if v>before_hard.get(k,0)},
            "worst_metric_regressions":metric_regressions,
            "legacy_hard_safety_before":pre["legacy_hard_safety"],
            "legacy_hard_safety_after":post["legacy_hard_safety"],
            "area_error":geometry["area_error"],"solver_cells":len(solver.cells),
            "record_count":len(registry["records"]),"deterministic":deterministic,
            "independent_reader":"check_openfoam2d.py",
            "specialized_hybrid_reader":"PASS" if name in ("circle","superellipse") else "H4-2-only; not applicable to terminated H4-3 fronts",
            "checkMesh":foam_status,"quality_contract_status":post["status"],
            "solver_sha256":sha(Path(f"{prefix}.hybrid.solver.cm2d")),
            "registry_sha256":sha(Path(f"{prefix}.hybrid.intersections.json")),
            "checkMesh_log_sha256":sha(foam_log) if foam_log.exists() else None}
        assert post["legacy_hard_safety"]["valid"]
        if name=="superellipse":
            assert entry["min_face_over_local_h_after"]>=.01
            assert entry["min_absolute_face_before"]<1.e-8
            assert entry["min_absolute_face_after"]>1.e-8
        result["cases"][name]=entry
        print(name,entry["min_face_over_local_h_before"],"->",entry["min_face_over_local_h_after"],foam_status,flush=True)
    result["all_cases_face_contract_pass"]=all(c["min_face_over_local_h_after"]>=.01 for c in result["cases"].values())
    result["q2_status"]="PASS" if result["all_cases_face_contract_pass"] and all(
        c["checkMesh"]=="PASS" and c["deterministic"] is True and not c["new_or_increased_hard_counts"] and not c["worst_metric_regressions"]
        for c in result["cases"].values()) else "PARTIAL_NOT_ACCEPTED"
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
    if args.require_full_acceptance and result["q2_status"]!="PASS":
        raise SystemExit("Q2 NOT ACCEPTED: see comparison.json; partial superellipse repair is not full Q2")

if __name__=="__main__":
    main()
