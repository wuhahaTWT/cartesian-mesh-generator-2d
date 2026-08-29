#!/usr/bin/env python3
"""Run the existing real OpenFOAM acceptance command on five Q2-A meshes."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
from generate_q1_baselines import CASES
from verify_q2a_construction import run

p=argparse.ArgumentParser()
p.add_argument("--evidence",type=Path,default=Path("build-q2a/evidence"))
p.add_argument("--output",type=Path,default=Path("artifacts/q2a/openfoam.json"))
a=p.parse_args()
repo=Path.cwd().resolve()
image="opencfd/openfoam-run:2606"
result={"image":image,"image_id":subprocess.check_output(
    ["docker","image","inspect",image,"--format","{{.Id}}"],text=True).strip(),"cases":{}}
for name in CASES:
    case=(a.evidence/f"{name}-case").resolve().relative_to(repo)
    command=["docker","run","--rm","-v",f"{repo}:/home/openfoam/workingDir",
             "-w","/home/openfoam/workingDir",image,"checkMesh","-case",
             "/home/openfoam/workingDir/"+str(case),"-writeAllFields"]
    log=a.evidence/f"{name}.checkMesh.log"
    seconds,text=run(command,log)
    assert "Mesh OK" in text and "Failed " not in text, f"{name}: checkMesh failed; see {log}"
    result["cases"][name]={"command":command,"status":"PASS","wall_seconds":seconds,
                           "log_sha256":hashlib.sha256(log.read_bytes()).hexdigest()}
    print(name,"Mesh OK",flush=True)
a.output.parent.mkdir(parents=True,exist_ok=True)
a.output.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n")
