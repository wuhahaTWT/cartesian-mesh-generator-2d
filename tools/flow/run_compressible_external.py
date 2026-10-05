#!/usr/bin/env python3
"""Generate a native conformal circular-body grid and run a finite-time external flow."""
import argparse
import json
import math
from pathlib import Path
import subprocess
import time
from run_compressible_channel import compress_histories, digest


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output", type=Path, required=True)
    p.add_argument("--far-spans", type=float, default=2.)
    p.add_argument("--wall-cells", type=float, default=16.)
    p.add_argument("--mach", type=float, default=.2)
    p.add_argument("--max-step", type=float, default=1e-8)
    p.add_argument("--end-time", type=float, default=2e-7)
    p.add_argument("--max-seconds", type=float, default=180.)
    args = p.parse_args()
    root = args.output.resolve();root.mkdir(parents=True, exist_ok=False)
    radius = 5e-5
    body = root/"body.xy"
    body.write_text("".join(f"{radius*math.cos(2*math.pi*i/48):.17g} {radius*math.sin(2*math.pi*i/48):.17g}\n" for i in range(48)))
    cli, flow = (Path("build")/x for x in ("cartmesh2d_cli", "cartmesh2d_euler_cli"))
    binaries = {str(x): digest(x) for x in (cli, flow)}
    report = {"scope": "finite-time laminar external start; no steady drag or external physical validation", "commands": [], "inputs": vars(args).copy(), "binarySha256": binaries}
    report["inputs"]["output"] = str(root)
    def save():
        (root/"run.json").write_text(json.dumps(report, indent=2)+"\n")
    def run(name, cmd):
        cmd = list(map(str, cmd));start = time.perf_counter()
        with (root/(name+".log")).open("w") as log:
            result = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT)
        report["commands"].append({"name": name, "argv": cmd, "returncode": result.returncode, "seconds": time.perf_counter()-start});save()
        if result.returncode:
            raise RuntimeError(f"{name} failed with {result.returncode}; retained logs/state in {root}")
    try:
        run("mesh", [cli, body, root/"mesh", 7, .25, .1, "exterior", root/"foam", 0, 0,
                     "--size-field", "--far-field-spans", args.far_spans, "--wall-cells-per-span", args.wall_cells,
                     "--cells-per-level", 3, "--far-level", 0, "--max-safe-wall-level", 11])
        rho = 101325/(287.05*300);speed = args.mach*math.sqrt(1.4*287.05*300)
        run("flow", [flow, "--mesh", root/"mesh.solver.cm2d", "--output", root/"flow", "--case", "external",
                     "--density", rho, "--pressure", 101325, "--u", speed, "--viscosity", 1.846e-5,
                     "--conductivity", .025759, "--wall-model", "no-slip", "--wall-thermal", "temperature", "--wall-value", 300,
                     "--flux", "hllc", "--order", 2, "--integrator", "sdirk2", "--max-step", args.max_step,
                     "--end-time", args.end_time, "--max-seconds", args.max_seconds, "--checkpoint-every", 20])
        report["native"] = json.loads((root/"flow.json").read_text());report["completed"] = True
    except Exception as error:
        report["failure"] = str(error);raise
    finally:
        report["histories"] = compress_histories(root);save()


if __name__ == "__main__":
    main()
