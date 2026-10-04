#!/usr/bin/env python3
"""Prepare and run a small real Cut-cell laminar channel through the native CLIs.

This is input/command orchestration and output readback, not another discretized
equation verifier. Existing native quality and solver failures remain fatal.
"""
import argparse
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path
import subprocess
import shutil
import time


def stream_digest(stream):
    value = hashlib.sha256()
    for chunk in iter(lambda: stream.read(1024*1024), b""):
        value.update(chunk)
    return value.hexdigest()


def digest(path):
    with path.open("rb") as stream:
        return stream_digest(stream)


def compress_histories(root):
    """Retain every step losslessly without leaving large raw CSV histories."""
    records = {}
    for path in sorted(root.glob("*.history.csv")):
        raw_hash, raw_bytes = digest(path), path.stat().st_size
        compressed = path.with_suffix(path.suffix+".gz")
        with path.open("rb") as source, compressed.open("xb") as target:
            with gzip.GzipFile(filename="", mode="wb", fileobj=target, mtime=0) as archive:
                shutil.copyfileobj(source, archive)
        with gzip.open(compressed, "rb") as restored:
            if stream_digest(restored) != raw_hash:
                raise IOError(f"compressed history readback failed: {path}")
        records[path.name] = {"file": compressed.name, "uncompressedSha256": raw_hash, "uncompressedBytes": raw_bytes, "compressedBytes": compressed.stat().st_size}
        path.unlink()
    return records


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build", type=Path, default=Path("build"))
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--mach", type=float, choices=(0.2, 0.5), default=0.2)
    parser.add_argument("--level", type=int, choices=(5, 6), default=5)
    parser.add_argument("--wall-temperature", type=float, choices=(300., 330.), default=300.)
    parser.add_argument("--shape", choices=("straight", "smooth"), default="straight")
    parser.add_argument("--flow-times", type=float, default=4.)
    args = parser.parse_args()
    if not math.isfinite(args.flow_times) or args.flow_times <= 0:
        parser.error("--flow-times must be finite and positive")
    root = args.output.resolve()
    root.mkdir(parents=True, exist_ok=False)
    build = args.build.resolve()
    cli = next((p for p in (build/"cartmesh2d_cli", build/"Release/cartmesh2d_cli.exe", build/"cartmesh2d_cli.exe") if p.is_file()), None)
    flow = next((p for p in (build/"cartmesh2d_euler_cli", build/"Release/cartmesh2d_euler_cli.exe", build/"cartmesh2d_euler_cli.exe") if p.is_file()), None)
    if cli is None or flow is None:
        raise FileNotFoundError("build cartmesh2d_cli and cartmesh2d_euler_cli first")
    report = {"format": "cartmesh2d-compressible-channel-v1", "qualification": "finite-time native workflow; not steady or mesh-independent qualification", "commands": [], "completed": False}
    def save():
        (root/"run.json").write_text(json.dumps(report, ensure_ascii=False, indent=2)+"\n")
    def run(name, command, expected=0):
        command = [str(x) for x in command]
        started = time.perf_counter()
        with (root/(name+".log")).open("w") as log:
            result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, check=False)
        report["commands"].append({"name": name, "argv": command, "returncode": result.returncode, "expectedReturncode": expected, "seconds": time.perf_counter()-started})
        save()
        if result.returncode != expected:
            raise RuntimeError(f"{name} returned {result.returncode}; see {root/(name+'.log')}")
    def summary(name):
        return json.loads((root/(name+".json")).read_text())

    height, length, temperature, pressure, gas_r, gamma = .0001, .0004, 300., 101325., 287.05, 1.4
    density = pressure/(gas_r*temperature)
    speed = args.mach*math.sqrt(gamma*gas_r*temperature)
    viscosity = 1.846e-5
    conductivity = viscosity*(gamma*gas_r/(gamma-1))/.72
    end_time = args.flow_times*length/speed
    report["inputs"] = {"heightM": height, "lengthM": length, "mach": args.mach, "reynoldsHeight": density*speed*height/viscosity, "wallTemperatureK": args.wall_temperature, "muPaS": viscosity, "kWmK": conductivity, "gamma": gamma, "gasR": gas_r, "referenceDensity": density, "referencePressurePa": pressure, "referenceTemperatureK": temperature, "referenceSpeedMs": speed, "endTimeS": end_time, "flowThroughTimes": args.flow_times, "shape": args.shape, "level": args.level, "pressureOutletBackflow": "explicit failure", "openDiffusion": "zero traction and Fourier heat flux"}
    # Explicit interior fluid region. The upper wall contracts smoothly by 8%;
    # its sampled polyline is passed unchanged to the native Cut-cell mesher.
    points = [(0., 0.), (length, 0.)]
    for i in range(16, -1, -1):
        x = length*i/16
        y = height*(1-(.08*math.sin(math.pi*i/16)**2 if args.shape == "smooth" else 0.))
        points.append((x, y))
    geometry = root/"channel.xy"
    geometry.write_text("".join(f"{x:.17g} {y:.17g}\n" for x, y in points))
    save()
    try:
        run("mesh", [cli, geometry, root/"mesh", args.level, .125, .1, "interior", root/"openfoam", args.level-1, 0])
        common = [flow, "--mesh", root/"mesh.solver.cm2d", "--case", "channel", "--outlet-pressure", pressure,
                  "--gamma", gamma, "--gas-r", gas_r, "--density", density, "--pressure", pressure, "--u", speed,
                  "--viscosity", viscosity, "--conductivity", conductivity, "--wall-model", "no-slip",
                  "--wall-thermal", "temperature", "--wall-value", args.wall_temperature, "--wall-gradient", "quadratic",
                  "--flux", "hllc", "--order", 2, "--cfl", .4, "--end-time", end_time,
                  "--max-steps", 100000, "--max-seconds", 240, "--checkpoint-every", 1000]
        run("full", [*common, "--output", root/"full"])
        full = summary("full")
        if not full["targetReached"] or full["time"] != end_time:
            raise RuntimeError("full run did not reach the physical endpoint")
        run("limited", [*common, "--output", root/"limited", "--max-steps", 8], expected=2)
        limited = summary("limited")
        if limited["targetReached"] or limited["acceptedSteps"] != 8 or "accepted-step budget" not in limited["failure"]:
            raise RuntimeError("iteration budget was not reported as a retained-state failure")
        run("resumed", [*common, "--output", root/"resumed", "--restart", root/"limited.checkpoint"])
        identical = (root/"full.checkpoint").read_bytes() == (root/"resumed.checkpoint").read_bytes()
        if not identical:
            raise RuntimeError("resumed final checkpoint differs from uninterrupted trajectory")
        # Also exercise custom boundary file parsing at the same eight-step
        # accepted endpoint without duplicating the long physical trajectory.
        custom = common.copy()
        custom[custom.index("channel")] = "custom"
        i = custom.index("--outlet-pressure")
        del custom[i:i+2]
        for flag in ("--wall-model", "--wall-thermal", "--wall-value"):
            i = custom.index(flag)
            del custom[i:i+2]
        run("custom", [*custom, "--boundary", root/"full.boundaries", "--output", root/"custom", "--max-steps", 8], expected=2)
        if (root/"limited.cells.csv").read_bytes() != (root/"custom.cells.csv").read_bytes():
            raise RuntimeError("custom pressure outlet roundtrip changed the accepted field")
        with (root/"full.cells.csv").open() as field:
            rows = list(csv.DictReader(field))
        report["fieldRanges"] = {key: [min(float(row[key]) for row in rows), max(float(row[key]) for row in rows)] for key in ("rho", "p", "temperature", "mach")}
        if any(not math.isfinite(v) for extent in report["fieldRanges"].values() for v in extent):
            raise RuntimeError("nonfinite exported field")
        report.update(completed=True, cells=full["cells"], faces=full["faces"], acceptedSteps=full["acceptedSteps"], rejectedCandidates=full["rejectedCandidates"], elapsedSeconds=full["elapsedSeconds"], finalCheckpointIdentical=identical, customBoundaryFieldIdentical=True)
    except Exception as error:
        report["failure"] = str(error)
        raise
    finally:
        report["compressedHistory"] = compress_histories(root)
        report["sha256"] = {p.name: digest(p) for p in sorted(root.iterdir()) if p.is_file() and p.name != "run.json"}
        report["binarySha256"] = {p.name: digest(p) for p in (cli, flow)}
        save()
    print(json.dumps({key: value for key, value in report.items() if key not in ("commands", "sha256", "binarySha256")}, indent=2))


if __name__ == "__main__":
    main()
