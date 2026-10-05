#!/usr/bin/env python3
"""Native cloud campaign: input preparation, process costs and real output readback.

No Python discretization of the governing equations. Every solve uses the C++
conservative operator and its unmodified topology/quality/balance checks.
"""
import argparse
import csv
import json
import math
from pathlib import Path
import subprocess
import sys
import time


def read_field(path):
    with path.open() as stream:
        return list(csv.DictReader(stream))


def compare(a, b):
    left, right = read_field(a), read_field(b)
    if len(left) != len(right):
        raise ValueError("different grids in temporal comparison")
    result = {}
    for key in ("rho", "u", "v", "p", "temperature"):
        weighted = area = maximum = 0.
        for x, y in zip(left, right):
            if (x["x"], x["y"], x["area"]) != (y["x"], y["y"], y["area"]):
                raise ValueError("different geometry in temporal comparison")
            delta = float(x[key])-float(y[key])
            volume = float(x["area"])
            area += volume
            weighted += volume*delta*delta
            maximum = max(maximum, abs(delta))
        result[key] = {"areaRms": math.sqrt(weighted/area), "maximum": maximum}
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    root = args.output.resolve()
    root.mkdir(parents=True, exist_ok=False)
    cases = [
        ("smooth-100ns", ["--shape", "smooth", "--max-step", "1e-7"]),
        ("smooth-50ns", ["--shape", "smooth", "--max-step", "5e-8"]),
        ("smooth-25ns", ["--shape", "smooth", "--max-step", "2.5e-8"]),
        ("smooth-reference", ["--shape", "smooth", "--integrator", "explicit", "--time-step-control", "stage-guarded", "--cfl", ".05"]),
        ("hot-20ns", ["--mach", ".5", "--wall-temperature", "330", "--max-step", "2e-8"]),
        ("hot-10ns", ["--mach", ".5", "--wall-temperature", "330", "--max-step", "1e-8"]),
        ("hot-5ns", ["--mach", ".5", "--wall-temperature", "330", "--max-step", "5e-9"]),
        ("hot-reference", ["--mach", ".5", "--wall-temperature", "330", "--integrator", "explicit", "--time-step-control", "stage-guarded", "--cfl", ".025"]),
        ("low-cold-total", ["--mach", ".05", "--wall-temperature", "270", "--inlet-model", "total", "--max-step", "2e-8"]),
        ("high-hot-total", ["--mach", ".7", "--wall-temperature", "330", "--inlet-model", "total", "--max-step", "1e-8"]),
        ("perturbed-smooth", ["--shape", "smooth", "--shape-phase", ".5", "--initial-pressure-perturbation", ".002", "--max-step", "2e-8"]),
        ("dense-cold", ["--shape", "smooth", "--shape-phase", "-.5", "--level", "6", "--wall-temperature", "280", "--max-step", "1e-8"]),
        ("long-outlet", ["--length", ".0006", "--shape", "smooth", "--outlet-ratio", ".995", "--max-step", "2e-8"]),
        ("short-outlet", ["--length", ".0003", "--outlet-ratio", "1.002", "--max-step", "2e-8"]),
        ("backflow-rejected", ["--outlet-ratio", "1.5", "--max-step", "1e-8"]),
    ]
    report = {"scope": "native Linux finite-time campaign; failures retained; no App qualification", "cases": {}, "comparisons": {}}
    def save():
        (root/"campaign.json").write_text(json.dumps(report, indent=2)+"\n")
    for name, options in cases:
        command = [sys.executable, "tools/flow/run_compressible_channel.py", "--output", str(root/name),
                   "--integrator", "sdirk2", "--flow-times", ".125", "--single-run", "--max-seconds", "180", *options]
        started = time.perf_counter()
        with (root/(name+".log")).open("w") as log:
            completed = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT)
        entry = {"argv": command, "returncode": completed.returncode, "completeProcessSeconds": time.perf_counter()-started}
        summary = root/name/"full.json"
        if summary.exists():
            entry["native"] = json.loads(summary.read_text())
            rows = read_field(root/name/"full.cells.csv")
            entry["ranges"] = {key: [min(float(r[key]) for r in rows), max(float(r[key]) for r in rows)] for key in ("rho", "p", "temperature", "mach", "area")}
        report["cases"][name] = entry
        save()
        print(name, completed.returncode, entry["completeProcessSeconds"], flush=True)
    for family, names in (("smooth", ("100ns", "50ns", "25ns")), ("hot", ("20ns", "10ns", "5ns"))):
        ref = root/(family+"-reference")/"full.cells.csv"
        if ref.exists() and report["cases"][family+"-reference"].get("native", {}).get("targetReached"):
            for suffix in names:
                path = root/(family+"-"+suffix)/"full.cells.csv"
                if path.exists() and report["cases"][family+"-"+suffix].get("native", {}).get("targetReached"):
                    report["comparisons"][family+"-"+suffix] = compare(path, ref)
    report["unexpectedFailures"] = [name for name, entry in report["cases"].items()
                                    if name != "backflow-rejected" and entry["returncode"] != 0]
    rejected = report["cases"]["backflow-rejected"].get("native", {})
    report["backflowRejectedWithInitialState"] = rejected.get("steps") == 0 and "backflow" in rejected.get("failure", "")
    save()
    if report["unexpectedFailures"] or not report["backflowRejectedWithInitialState"]:
        raise RuntimeError("campaign retained unexpected failures; inspect campaign.json")


if __name__ == "__main__":
    main()
