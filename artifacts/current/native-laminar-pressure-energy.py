#!/usr/bin/env python3
"""Run native operator restrictions and analyze their exported matrices.

Python performs linear algebra and bookkeeping only. Every matrix coefficient
and face momentum flux comes from the native product or its native research
adapter; this is not an independently reconstructed flow solver.
"""
import hashlib
import json
import os
from pathlib import Path
import subprocess

os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import numpy as np
from scipy.linalg import eigh, eigvals
from scipy.sparse import coo_matrix

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "outputs/cloud-laminar/viscous-energy"
CYL = ROOT / "outputs/cloud-laminar/cylinder-joint"
MESHES = {"bad": CYL/"far-20-fixed128-grid1.solver.cm2d",
          "coarse": CYL/"far-20-fixed128-grid0.solver.cm2d",
          "repaired": CYL/"far-20-fixed128-grid1-volume-075.solver.cm2d",
          "fine": CYL/"far-20.solver.cm2d"}


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def relative(p):
    return str(p.relative_to(ROOT))


def run(command, output=None):
    if output:
        with output.open("w") as stream:
            subprocess.run(command, cwd=ROOT, stdout=stream, check=True)
    else:
        subprocess.run(command, cwd=ROOT, check=True)


def spectrum(name, mode):
    path = OUT/(name+".matrix.csv")
    cells_path = OUT/(name+".cells.csv")
    a = np.loadtxt(path, delimiter=",", skiprows=1)
    cells = np.loadtxt(cells_path, delimiter=",", skiprows=1)
    n = len(cells)*(1 if mode == "pressure" else 2)
    if int(a[:, 0].max()) != n-1 or int(a[:, 1].max()) != n-1:
        raise RuntimeError("incomplete native matrix")
    b = coo_matrix((a[:, 2], (a[:, 0].astype(int), a[:, 1].astype(int))), shape=(n,n)).toarray()
    h = (b+b.T)/2
    w, v = eigh(h)
    e = eigvals(b)
    critical = 0 if mode == "pressure" else -1
    eigenvector = v[:, critical]
    np.savetxt(OUT/(name+".critical-vector.csv"),eigenvector,delimiter=",")
    return {"matrix":relative(path),"matrixSha256":sha(path),"cellsSha256":sha(cells_path),
            "supportCells":len(cells),"dimension":n,"symmetricMinimum":float(w[0]),"symmetricMaximum":float(w[-1]),
            "eigenvalueMinimumRealPart":float(e.real.min()),"eigenvalueMaximumRealPart":float(e.real.max()),
            "criticalSymmetricEigenResidual":float(np.linalg.norm(h@eigenvector-w[critical]*eigenvector)),
            "interpretation":"restricted geometric pressure-flux block at uniform rAU=1 s, not full saddle Jacobian" if mode=="pressure" else "restricted viscous perturbation generator at nu=1 m2/s; eigenvalues in 1/s"}


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    for name in ("viscous-energy","branch-energy"):
        run(["g++","-std=c++20","-O2","-Wall","-Wextra","-Wpedantic","-Werror","-Iinclude",f"artifacts/current/native-laminar-{name}.cpp","build/libcartmesh2d_fv.a","build/libcartmesh2d.a","-o",f"build/native-laminar-{name}"])
    evidence={"schema":1,"sourceBase":"37d8bbade12d6caf5472a1aa6bc16c2fa8cb655e","sourceBaseTree":"d212321ff9fc7e436d814fdd8df371fe62ee98a7",
              "scope":"Mechanism research; no change to product equations, defaults, gates, boundaries or candidate selection",
              "environment":{"platform":"Linux x86_64","compiler":"GCC 13.3.0"},
              "supportDefinition":"all embedded-wall cells plus complete graph rings; omitted perturbation DOFs are zero, physical cells and mesh are not removed",
              "meshSha256":{name:sha(path) for name,path in MESHES.items()},"spectra":{},"powerComparisons":{}}
    cases=[("bad","symmetric",1), ("bad","symmetric",3),("bad","laplacian",1),("bad","orthogonal",1),
           ("coarse","symmetric",1),("repaired","symmetric",1),("fine","symmetric",1),
           ("bad","pressure",3),("coarse","pressure",3),("repaired","pressure",3)]
    for mesh,mode,rings in cases:
        name=f"{mesh}-{mode}-ring{rings}"
        run(["build/native-laminar-viscous-energy",relative(MESHES[mesh]),relative(OUT/name),str(rings),mode],OUT/(name+".log"))
        evidence["spectra"][name]=spectrum(name,mode)
    comparisons=[("bad-branches","bad","far-20-fixed128-grid1-face-limited-linear","far-20-fixed128-grid1-short-nu01-from-nu1"),
                 ("coarse-paths","coarse","far-20-fixed128-grid0-face-limited-linear","far-20-fixed128-grid0-short-nu01-from-nu1"),
                 ("bad-anderson","bad","far-20-fixed128-grid1-anderson-flat","far-20-fixed128-grid1-short-nu01-from-nu1")]
    for name,mesh,left,right in comparisons:
        output=OUT/(name+".json")
        run(["build/native-laminar-branch-energy",relative(MESHES[mesh]),relative(CYL/left),relative(CYL/right),relative(OUT/(name+".csv"))],output)
        result=json.loads(output.read_text())
        result["rawSha256"]={str(CYL/prefix)+suffix:sha(CYL/(prefix+suffix)) for prefix in (left,right) for suffix in (".cells.csv",".faces.csv")}
        # Store repository-relative provenance; retain every native cell.
        result["rawSha256"]={str(Path(p).relative_to(ROOT)):digest for p,digest in result["rawSha256"].items()}
        evidence["powerComparisons"][name]=result
    run(["build/native-laminar-viscous-energy",relative(MESHES["bad"]),relative(OUT/"bad"),"0","affine"],OUT/"affine.json")
    evidence["affinePressurePatch"]=json.loads((OUT/"affine.json").read_text())
    ablation=ROOT/"outputs/cloud-laminar/pressure-pair-ablation/runs.json"
    if ablation.exists():
        record=json.loads(ablation.read_text())
        evidence["centralPairAblation"]={"scope":record["scope"],"binarySha256":record["binarySha256"],"runs":[]}
        for r in record["runs"]:
            s=r.get("summary",{})
            prefix=Path(r["command"][-1])
            r["rawSha256"]={str(prefix.relative_to(ROOT))+suffix:sha(Path(str(prefix)+suffix)) for suffix in (".json",".cells.csv",".faces.csv",".residuals.csv")}
            evidence["centralPairAblation"]["runs"].append({"rawSha256":r["rawSha256"],"exitCode":r["exitCode"],"seconds":r["seconds"],"meshSha256":r["meshSha256"],"summary":{k:s.get(k) for k in ("status","converged","cells","coupledEvaluations","maximumSpeedRatio","pressureRangeRatio","wallForceX","pressureForceX","wallViscousForceX","globalRelativeImbalance")}})
    evidence["qualification"]=[
        "Negative largest restricted viscous energy eigenvalues exclude compact perturbation viscous anti-dissipation on tested supports only, not whole-domain incompressible stability.",
        "The pressure-flux block is not positive semidefinite on both abnormal and normal meshes; its extremal eigenvalue alone is not a quality gate.",
        "Native same-equation branch differences close their complete cell momentum-power balance. Pressure power sustains the difference against viscous and convective removal; this is a numerical compatibility finding, not proof that either branch is the physical solution.",
        "The central-pair ablation is deliberately not affine exact on unequal cut cells and cannot be promoted as a general fix even if field amplitudes improve.",
        "No cells, local extrema, boundary atoms or solver gates were altered or hidden. No new macOS, frontend, Electron or packaged-app validation."]
    evidence["researchBinarySha256"]={name:sha(ROOT/"build"/name) for name in ("native-laminar-viscous-energy","native-laminar-branch-energy")}
    (ROOT/"artifacts/current/native-laminar-pressure-energy.json").write_text(json.dumps(evidence,indent=2)+"\n")


if __name__=="__main__":
    main()
