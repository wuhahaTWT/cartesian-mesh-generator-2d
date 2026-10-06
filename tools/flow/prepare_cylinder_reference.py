#!/usr/bin/env python3
"""Prepare a native Re-controlled cylinder for a low-Mach wake-frequency comparison.

Runs only the native mesher and boundary exporter. The saved flow command uses
the existing Euler research driver; no flow solution or qualification is claimed.
"""
import argparse
import json
import math
from pathlib import Path
import subprocess
import time

from run_compressible_channel import digest

REPORT = "https://ntrs.nasa.gov/api/citations/19930092207/downloads/19930092207.pdf"
NASA_CASE = "https://www.grc.nasa.gov/www/wind/valid/lamcyl/Study1_files/Study1.html"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--build", type=Path, default=Path("build"))
    parser.add_argument("--reynolds", type=float, default=100.)
    parser.add_argument("--mach", type=float, default=.1)
    parser.add_argument("--diameter", type=float, default=1e-4)
    parser.add_argument("--segments", type=int, default=192)
    parser.add_argument("--phase", type=float, default=0., help="fraction of one contour segment; .17 currently has a retained native quality failure")
    parser.add_argument("--wall-cells", type=float, default=48.)
    parser.add_argument("--far-spans", type=float, default=15.5)
    parser.add_argument("--end-transits", type=float, default=200.)
    parser.add_argument("--step-transits", type=float, default=.02)
    parser.add_argument("--wall-budget", type=float, default=1800.)
    args = parser.parse_args()
    positive = (args.reynolds, args.mach, args.diameter, args.wall_cells,
                args.far_spans, args.end_transits, args.step_transits, args.wall_budget)
    if not all(math.isfinite(x) and x > 0 for x in positive):
        parser.error("all physical, resolution and runtime inputs must be finite and positive")
    if not 50 < args.reynolds < 150 or not 0 < args.mach < 1:
        parser.error("this comparison uses the experimental 50 < Re < 150 fit and subsonic inflow")
    if not 16 <= args.segments <= 4096 or not math.isfinite(args.phase) or not 0 <= args.phase < 1:
        parser.error("segments must be in [16, 4096] and phase in [0, 1)")
    if args.step_transits > args.end_transits:
        parser.error("step-transits cannot exceed end-transits")
    build = args.build.resolve()
    binaries = {name: build/name for name in
                ("cartmesh2d_cli", "cartmesh2d_euler_cli", "cartmesh2d_euler_diffusion_benchmark")}
    for binary in binaries.values():
        if not binary.is_file():
            parser.error(f"build the required native executable first: {binary}")
    root = args.output.resolve()
    root.mkdir(parents=True, exist_ok=False)
    # Fixed constants match the existing coupled-diffusion native driver.
    gamma, gas_r, temperature = 1.4, 287.05, 300.
    viscosity, conductivity = 1.846e-5, .025758750694444447
    speed = args.mach * math.sqrt(gamma * gas_r * temperature)
    density = args.reynolds * viscosity / (speed * args.diameter)
    pressure = density * gas_r * temperature
    transit = args.diameter / speed
    body = root/"body.xy"
    body.write_text("".join(
        f"{.5*args.diameter*math.cos(2*math.pi*(i+args.phase)/args.segments):.17g} "
        f"{.5*args.diameter*math.sin(2*math.pi*(i+args.phase)/args.segments):.17g}\n"
        for i in range(args.segments)))
    mesh, boundary = root/"mesh.solver.cm2d", root/"flow.boundaries"
    flow = [str(binaries["cartmesh2d_euler_diffusion_benchmark"]), str(mesh), str(boundary),
            str(root/"flow"), "corrected", repr(args.end_transits*transit),
            repr(args.step_transits*transit), repr(args.wall_budget)]
    report = {
        "format": "cartmesh2d-cylinder-reference-input-v1", "prepared": False,
        "flowExecuted": False, "physicalQualification": False,
        "scope": "Constant-property subsonic continuum-model wake; low-speed experimental fit is a comparison target, not an exact solution or automatic pass gate.",
        "requested": {**vars(args), "output": str(root), "build": str(build)},
        "physics": {"gamma": gamma, "gasConstant": gas_r, "temperatureK": temperature,
                    "densityKgM3": density, "pressurePa": pressure, "velocityMS": speed,
                    "dynamicViscosityPaS": viscosity, "thermalConductivityWmK": conductivity,
                    "diameterM": args.diameter, "reynolds": args.reynolds, "mach": args.mach,
                    "diameterTransitSeconds": transit, "wall": "static no-slip adiabatic",
                    "prandtl": viscosity*gamma*gas_r/((gamma-1)*conductivity)},
        "reference": {
            "source": REPORT, "reportPages": [8, 11], "equation": "2a",
            "strouhalDefinition": "f*D/U", "empiricalStrouhal": .212*(1-21.2/args.reynolds),
            "validReynoldsIntervalExclusive": [50, 150],
            "reportedPrecision": "Author estimates best-fit line accuracy at 1%; not a confidence interval for every measurement. End effects were not corrected.",
            "limits": "Low-speed wind-tunnel correlation; finite Mach, wall model, contour, domain, space/time resolution and sampled cycles must be assessed separately.",
            "fourPercentIsNotLowReScatter": "Report p11 refers to extrapolating equation 2b to Re=10000, not uncertainty of equation 2a.",
            "independentNumericalCase": NASA_CASE,
            "nasaCaseIsDifferent": "NASA case uses Mach .2, Re150 and a far boundary 200D from the centre; its drag table is not adopted as a target for this case."
        },
        "recommendedMeasurement": {
            "signal": "Sum native wall momentum-y flux by accepted physical time; use lift or off-axis wake velocity, not drag's doubled frequency.",
            "normalization": "Cd,Cl = wall force per unit depth / (0.5*rho*U^2*D)",
            "requirements": "Retain startup; select a settled periodic interval explicitly; compare successive cycle amplitudes and periods, space/time/domain changes and full process cost.",
            "noPeriodicStateYet": True
        },
        "nativeFlowCommand": flow,
        "nativeFlowEnvironment": {"CARTMESH_RESEARCH_PRECONDITIONER": "diagonal"},
        "restart": "Use a new output prefix and append the previous native .checkpoint path; retain each original history and checkpoint.",
        "commands": [], "binarySha256": {k: digest(v) for k, v in binaries.items()},
        "preparerSha256": digest(Path(__file__).resolve())
    }
    def save():
        (root/"case.json").write_text(json.dumps(report, ensure_ascii=False, indent=2)+"\n")
    def run(name, command):
        command = list(map(str, command))
        start = time.perf_counter()
        with (root/(name+".log")).open("w") as log:
            result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT)
        report["commands"].append({"name": name, "argv": command, "returncode": result.returncode,
                                   "seconds": time.perf_counter()-start})
        save()
        if result.returncode:
            raise RuntimeError(f"native {name} failed; all inputs and diagnostics retained in {root}")
    save()
    try:
        run("mesh", [binaries["cartmesh2d_cli"], body, root/"mesh", 7, .25, .1, "exterior", root/"foam", 0, 0,
                     "--size-field", "--far-field-spans", args.far_spans, "--wall-cells-per-span", args.wall_cells,
                     "--cells-per-level", 3, "--far-level", 0, "--max-safe-wall-level", 11])
        run("boundaries", [binaries["cartmesh2d_euler_cli"], "--mesh", mesh, "--case", "external",
                          "--density", density, "--u", speed, "--v", 0, "--pressure", pressure,
                          "--gamma", gamma, "--gas-r", gas_r, "--viscosity", viscosity,
                          "--conductivity", conductivity, "--wall-model", "no-slip",
                          "--wall-thermal", "insulated", "--export-boundaries", boundary])
        report["inputSha256"] = {x.name: digest(x) for x in (body, mesh, boundary)}
        report["meshResolution"] = json.loads((root/"mesh.resolution.json").read_text())
        report["meshQuality"] = json.loads((root/"mesh.construction-quality.json").read_text())
        report["prepared"] = True
    except Exception as error:
        report["failure"] = str(error)
        raise
    finally:
        save()
    print(json.dumps({"case": str(root/"case.json"), "prepared": True, "flowExecuted": False,
                      "reynolds": args.reynolds, "empiricalStrouhal": report["reference"]["empiricalStrouhal"]}))


if __name__ == "__main__":
    main()
