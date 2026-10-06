#!/usr/bin/env python3
"""Research-only central adjoint pressure/divergence ablation.

Build one replacement native object in outputs; do not change product sources
or the normal binary. This deliberately removes distance/skew linear exactness
and is NOT a qualified discretization or proposed default. It isolates the
pressure/velocity pairing found in the native branch-energy experiment.
"""
import hashlib
import json
from pathlib import Path
import subprocess
import time

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "outputs/cloud-laminar/pressure-pair-ablation"


def replace_once(text, old, new):
    if text.count(old) != 1:
        raise RuntimeError(f"source anchor changed: {old[:90]}")
    return text.replace(old, new)


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    header_path = ROOT / "include/cartmesh2d/fv/detail/FlowFaceOperators2D.hpp"
    source_path = ROOT / "src/fv/Incompressible2D.cpp"
    header = header_path.read_text()
    # Only the pressure shared-face routine is replaced. Viscous gradients,
    # convection, pressure compact normal derivative and all gates are intact.
    start = header.index("inline std::vector<double> pressureFaceValues(")
    end = header.index("inline std::vector<Vector2D> conservativePressureGradient(", start)
    header = header[:start] + '''inline std::vector<double> pressureFaceValues(
    const FvMesh2D& mesh, const std::vector<double>& p,
    const std::vector<Vector2D>&, const std::vector<double>& boundary,
    const std::vector<bool>& fixed) {
    std::vector<double> result(mesh.faces.size());
    for(std::size_t id=0;id<mesh.faces.size();++id) {
        const auto& f=mesh.faces[id];
        result[id]=f.neighbour?.5*(p[f.owner]+p[*f.neighbour]):(fixed[id]?boundary[id]:p[f.owner]);
    }
    return result;
}

''' + header[end:]
    overlay = OUT / "include/cartmesh2d/fv/detail/FlowFaceOperators2D.hpp"
    overlay.parent.mkdir(parents=True, exist_ok=True)
    overlay.write_text(header)
    source = source_path.read_text()
    source = replace_once(source,
        "const double uf=interpolate(f,r.u)+dot(interpolateGradient(f,gup),skew),vf=interpolate(f,r.v)+dot(interpolateGradient(f,gvp),skew);",
        "const double uf=.5*(r.u[i]+r.u[j]),vf=.5*(r.v[i]+r.v[j]);")
    # Both physical-time and steady under-relaxation defects must use exactly
    # the same face interpolation as the predictor; retain their existing math.
    for suffix in (",skew);",):
        old_u = "const double oldUf=interpolate(f,oldU)+dot(interpolateGradient(f,gu)" + suffix
        old_v = "const double oldVf=interpolate(f,oldV)+dot(interpolateGradient(f,gv)" + suffix
        if source.count(old_u) != 2 or source.count(old_v) != 2:
            raise RuntimeError("old velocity interpolation anchors changed")
        source = source.replace(old_u,"const double oldUf=.5*(oldU[i]+oldU[j]);")
        source = source.replace(old_v,"const double oldVf=.5*(oldV[i]+oldV[j]);")
    copied = OUT / "Incompressible2D.cpp"
    copied.write_text(source)
    commands = [
        ["g++", "-std=c++20", "-O3", "-DNDEBUG", "-I"+str(OUT/"include"), "-Iinclude", "-c", str(copied), "-o", str(OUT/"Incompressible2D.o")],
        ["g++", "-O3", "build/CMakeFiles/cartmesh2d_flow_cli.dir/apps/cartmesh2d_flow_cli.cpp.o", str(OUT/"Incompressible2D.o"), "build/libcartmesh2d_fv.a", "build/libcartmesh2d.a", "-o", str(OUT/"flow-cli-central-pair")],
    ]
    with (OUT/"build.log").open("w") as log:
        for command in commands:
            subprocess.run(command,cwd=ROOT,stdout=log,stderr=subprocess.STDOUT,check=True)
    records={"scope":"research ablation; sacrifices non-midpoint affine consistency; not product acceptance", "productSourceSha256":{str(p.relative_to(ROOT)):sha(p) for p in (header_path,source_path)},"binarySha256":sha(OUT/"flow-cli-central-pair"),"commands":commands,"runs":[]}
    for grid in (0,1):
        mesh=f"outputs/cloud-laminar/cylinder-joint/far-20-fixed128-grid{grid}.solver.cm2d"
        prefix=OUT/f"grid{grid}"
        command=[str(OUT/"flow-cli-central-pair"),"--mesh",mesh,"--case","external","--nu",".1","--speed","1","--convection","face-limited-linear","--tolerance","1e-8","--max-iterations","2500","--output",str(prefix)]
        started=time.monotonic()
        with prefix.with_suffix(".log").open("w") as log:
            result=subprocess.run(command,cwd=ROOT,stdout=log,stderr=subprocess.STDOUT)
        record={"command":command,"exitCode":result.returncode,"seconds":time.monotonic()-started,"meshSha256":sha(ROOT/mesh)}
        if prefix.with_suffix(".json").exists():record["summary"]=json.loads(prefix.with_suffix(".json").read_text())
        records["runs"].append(record)
        (OUT/"runs.json").write_text(json.dumps(records,indent=2)+"\n")
        print(grid,result.returncode,record["seconds"],flush=True)


if __name__ == "__main__":
    main()
