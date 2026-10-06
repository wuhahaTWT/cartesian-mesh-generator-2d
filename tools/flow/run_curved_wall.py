#!/usr/bin/env python3
"""Prepare real annular cut cells and invoke native spatial/coupled solvers.

No Python PDE reconstruction. Radii .5/1 m, walls 2/2.2 K, R=1 and
k=.37 are verification parameters, not an air application recommendation.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path
import subprocess
import time
from run_compressible_channel import compress_histories


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--build", type=Path, default=Path("build"))
    p.add_argument("--segments", type=int, default=128)
    p.add_argument("--phase", type=float, default=.375)
    p.add_argument("--levels", type=int, nargs="+", default=[4, 5, 6])
    p.add_argument("--include-quadratic", action="store_true",
                   help="diagnose optional quadratic walls, including known unstable coarse cases")
    p.add_argument("--skip-coupled", action="store_true")
    a = p.parse_args()
    if not 16 <= a.segments <= 1024 or not math.isfinite(a.phase) or any(n not in range(4, 8) for n in a.levels):
        p.error("unsupported geometry/levels")
    root = a.output.resolve();root.mkdir(parents=True, exist_ok=False)
    build = a.build.resolve()
    binaries = {n: hashlib.sha256((build/n).read_bytes()).hexdigest()
                for n in ["cartmesh2d_cli", "cartmesh2d_curved_heat_benchmark", "cartmesh2d_euler_cli"]}
    report = {"scope": "NEW native annulus research; spatial isolation is not physical qualification",
              "binarySha256": binaries, "commands": [], "cases": {}, "failures": []}
    def save():
        (root/"run.json").write_text(json.dumps(report, indent=2)+"\n")
    def run(name, cmd):
        cmd = list(map(str, cmd));start = time.perf_counter()
        with (root/(name+".log")).open("w") as f:
            result = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT)
        report["commands"].append({"name": name, "argv": cmd, "returncode": result.returncode,
                                   "completeProcessSeconds": time.perf_counter()-start})
        if result.returncode: report["failures"].append(name)
        save();return result.returncode
    try:
        for level in a.levels:
            d = root/f"L{level}";d.mkdir();xy = d/"annulus.xy"
            xy.write_text("\n\n".join("\n".join(
                f"{r*math.cos(2*math.pi*(i+a.phase)/a.segments):.17g} {r*math.sin(2*math.pi*(i+a.phase)/a.segments):.17g}"
                for i in range(a.segments)) for r in [1., .5])+"\n")
            if run(f"L{level}-mesh", [build/"cartmesh2d_cli", xy, d/"mesh", level, .125, .1,
                    "interior", d/"foam", level-1, 0]): continue
            for scheme in (["linear", "quadratic"] if a.include_quadratic else ["linear"]):
                for mode in ["trace", "isothermal"]:
                    name = f"L{level}-{scheme}-{mode}";prefix = d/f"{scheme}-{mode}"
                    if not run(name, [build/"cartmesh2d_curved_heat_benchmark", d/"mesh.solver.cm2d", prefix, scheme, mode]):
                        report["cases"][name] = json.loads(prefix.with_suffix(".json").read_text());save()
            if level == min(a.levels) and not a.skip_coupled:
                run(f"L{level}-coupled", [build/"cartmesh2d_euler_cli", "--mesh", d/"mesh.solver.cm2d",
                    "--case", "custom", "--boundary", d/"linear-trace.boundaries", "--output", d/"coupled",
                    "--gas-r", 1, "--gamma", 1.4, "--density", 1, "--pressure", 2.1, "--u", 0,
                    "--conductivity", .37, "--viscosity", .02, "--wall-gradient", "linear",
                    "--flux", "hllc", "--order", 2, "--integrator", "sdirk2", "--max-step", .02,
                    "--end-time", 4, "--max-seconds", 300, "--checkpoint-every", 20])
    finally:
        for d in root.glob("L*"):
            if d.is_dir(): compress_histories(d)
        report["filesSha256"] = {str(f.relative_to(root)): hashlib.sha256(f.read_bytes()).hexdigest()
                                for f in root.rglob("*") if f.is_file() and "foam" not in f.parts and f.name != "run.json"}
        save()
    return bool(report["failures"])


if __name__ == "__main__":
    raise SystemExit(main())
