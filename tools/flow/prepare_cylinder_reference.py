#!/usr/bin/env python3
"""Prepare a native cylinder for a selected low-Mach reference comparison.

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
FORNBERG = "https://www.colorado.edu/amath/sites/default/files/attached-files/jcp_80_fl_p_cyl.pdf"


def steady_commands(binaries, root, physics, end_time, maximum_step, wall_budget):
    """Describe existing native CLI/diagnostic operations without running a flow."""
    mesh, boundary = root/"mesh.solver.cm2d", root/"flow.boundaries"
    prefix = root/"steady"
    command = [str(binaries["cartmesh2d_euler_cli"]), "--mesh", str(mesh),
               "--output", str(prefix), "--case", "external"]
    for flag, key in (("--density", "densityKgM3"), ("--u", "velocityMS"),
                      ("--pressure", "pressurePa"), ("--gamma", "gamma"),
                      ("--gas-r", "gasConstant"), ("--viscosity", "dynamicViscosityPaS"),
                      ("--conductivity", "thermalConductivityWmK")):
        command += [flag, repr(physics[key])]
    command += ["--v", "0", "--wall-model", "no-slip", "--wall-thermal", "insulated",
                "--flux", "hllc", "--order", "2", "--integrator", "sdirk2",
                "--mode", "steady", "--steady-scale", repr(physics["diameterTransitSeconds"]),
                "--steady-tolerance", "1e-5", "--end-time", repr(end_time),
                "--max-step", repr(maximum_step), "--max-seconds", repr(wall_budget)]
    snapshot = [str(binaries["cartmesh2d_euler_diffusion_benchmark"]), str(mesh), str(boundary),
                str(root/"steady-instantaneous"), "corrected", repr(end_time),
                repr(maximum_step), repr(wall_budget), str(prefix)+".checkpoint"]
    return {
        "nativeSteadyCommand": command,
        "nativeSteadyRequirements": "Uses the existing three scaled stopping criteria at 1e-5; this is a numerical stopping request, not a drag-accuracy gate. Endpoint or budget exhaustion preserves the last accepted state and returns failure when the criteria are unmet.",
        "nativeInstantaneousCommand": snapshot,
        "nativeInstantaneousEnvironment": {
            "CARTMESH_RESEARCH_CHECKPOINT_CASE": "external",
            "CARTMESH_RESEARCH_SNAPSHOT": "1",
            "CARTMESH_RESEARCH_RECONSTRUCTION": "1"
        },
        "nativeInstantaneousRequirements": "Read only the saved CLI physical state, export native instantaneous fluxes to a distinct prefix. SNAPSHOT skips integration; RECONSTRUCTION enables flux export. This does not establish steady convergence."
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--build", type=Path, default=Path("build"))
    parser.add_argument("--reference", choices=("roshko-frequency", "fornberg-steady-drag"),
                        default="roshko-frequency")
    parser.add_argument("--reynolds", type=float, help="defaults to 100 for Roshko, 20 for Fornberg")
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
    steady_reference = args.reference == "fornberg-steady-drag"
    if args.reynolds is None:
        args.reynolds = 20. if steady_reference else 100.
    positive = (args.reynolds, args.mach, args.diameter, args.wall_cells,
                args.far_spans, args.end_transits, args.step_transits, args.wall_budget)
    if not all(math.isfinite(x) and x > 0 for x in positive):
        parser.error("all physical, resolution and runtime inputs must be finite and positive")
    if not 0 < args.mach < 1:
        parser.error("this preparer requires subsonic inflow")
    if steady_reference and args.reynolds != 20.:
        parser.error("the steady drag comparison currently supports only the stable Re=20 reference")
    if not steady_reference and not 50 < args.reynolds < 150:
        parser.error("the Roshko frequency comparison requires 50 < Re < 150")
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
        "scope": "Constant-property subsonic continuum-model cylinder; the selected reference is a comparison target, not an exact solution or automatic pass gate.",
        "requested": {**vars(args), "output": str(root), "build": str(build)},
        "physics": {"gamma": gamma, "gasConstant": gas_r, "temperatureK": temperature,
                    "densityKgM3": density, "pressurePa": pressure, "velocityMS": speed,
                    "dynamicViscosityPaS": viscosity, "thermalConductivityWmK": conductivity,
                    "diameterM": args.diameter, "reynolds": args.reynolds, "mach": args.mach,
                    "diameterTransitSeconds": transit, "wall": "static no-slip adiabatic",
                    "prandtl": viscosity*gamma*gas_r/((gamma-1)*conductivity)},
        "reference": {
            "kind": args.reference,
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
    if steady_reference:
        report["reference"] = {
            "kind": args.reference, "source": FORNBERG,
            "title": "A numerical study of steady viscous flow past a circular cylinder",
            "author": "Bengt Fornberg", "year": 1980,
            "doi": "10.1017/S0022112080000419", "journalPages": [844, 845, 846],
            "type": "incompressible numerical reference", "reynolds": 20.,
            "dragCoefficient": 2.0001, "authorEstimatedNumericalError": .0002,
            "dragDefinition": "force per unit depth / (rho*U^2*radius) = force per unit depth / (0.5*rho*U^2*D)",
            "reportedPrecision": "Author's estimate is not a complete uncertainty bound: it excludes possible effects of a restrictive upstream computational region (p846).",
            "limits": "Incompressible steady result with specialized far-field treatment; finite Mach, polygon geometry, finite domain/open boundaries and space/time errors must be assessed separately. The reported estimate is not our acceptance threshold.",
            "noFrequencyReference": True
        }
        report["recommendedMeasurement"] = {
            "signal": "Sum native instantaneous wall momentum-x flux per unit depth, retaining pressure and viscous contributions separately.",
            "normalization": "Cd = wall force per unit depth / (0.5*rho*U^2*D)",
            "requirements": "Establish steady residual, field-change and force histories; a physical-time endpoint or exhausted budget is not convergence. Compare mesh, domain and Mach sensitivity before attributing a difference from the incompressible reference to spatial error.",
            "noSteadyStateYet": True
        }
        report.update(steady_commands(binaries, root, report["physics"],
                                      args.end_transits*transit, args.step_transits*transit,
                                      args.wall_budget))
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
    summary = {"case": str(root/"case.json"), "prepared": True, "flowExecuted": False,
               "reynolds": args.reynolds, "reference": args.reference}
    if steady_reference:
        summary["referenceDragCoefficient"] = report["reference"]["dragCoefficient"]
    else:
        summary["empiricalStrouhal"] = report["reference"]["empiricalStrouhal"]
    print(json.dumps(summary))


if __name__ == "__main__":
    main()
