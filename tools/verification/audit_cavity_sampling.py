#!/usr/bin/env python3
"""Recheck archived cavity fields with current and legacy sampling, without solving.

Exit 1 retains a failed original-threshold sequence. Never changes saved fields,
source summaries, reference values, or numerical tolerances to obtain a pass.
"""
import argparse
import copy
import hashlib
import json
import math
from pathlib import Path
import shlex
import sys

import verify_native_flow as verifier


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def calibration():
    f = lambda x, y: .7 + 1.3*x - .8*y
    errors, old_errors, curvature = [], [], []
    for n in (14, 30, 62):
        for scale in (1e-6, 1., 1e6):
            rows = [dict(x=scale*(i+.5)/n, y=scale*(j+.5)/n,
                         u=f((i+.5)/n, (j+.5)/n)) for j in range(n) for i in range(n)]
            walls = [dict(x=.5*scale, y=y*scale, u=f(.5, y)) for y in (0., 1.)]
            for y, _ in verifier.GHIA_U[1:-1]:
                errors.append(verifier.affine_sample(rows, .5*scale, y*scale, "u", boundary=walls)-f(.5, y))
                old_errors.append(verifier.idw(rows, .5*scale, y*scale, "u", boundary=walls)-f(.5, y))
        rows = [dict(x=(i+.5)/n, y=(j+.5)/n, u=((i+.5)/n)**2+((j+.5)/n)**2)
                for j in range(n) for i in range(n)]
        walls = [dict(x=.5, y=y, u=.25+y*y) for y in (0., 1.)]
        diffs = [verifier.affine_sample(rows, .5, y, "u", boundary=walls)-.25-y*y
                 for y, _ in verifier.GHIA_U[1:-1]]
        curvature.append(math.sqrt(math.fsum(e*e for e in diffs)/len(diffs)))
    return dict(affineMaxAbs=max(map(abs, errors)), legacyAffineMaxAbs=max(map(abs, old_errors)),
                quadraticGridSizes=[14, 30, 62], quadraticRms=curvature,
                quadraticOrders=[math.log(curvature[i]/curvature[i+1])/math.log(b/a)
                                 for i, (a, b) in enumerate(((14, 30), (30, 62)))])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", type=Path, action="append", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--max-iterations", type=int, default=4000,
                        help="original producer iteration budget, not a new solve")
    args = parser.parse_args()
    controls = verifier.argument_parser().parse_args([])
    controls.max_iterations = args.max_iterations
    report = dict(format="cartmesh2d-cavity-sampling-audit-v1", valid=False,
                  execution="independent readback only; no solver run", cases=[], issues=[],
                  command=shlex.join(sys.argv), sources=[], sequenceChecks={},
                  verifierSourceSha256=sha256(Path(verifier.__file__)),
                  producerSourceSha256=sha256(Path(__file__)),
                  thresholds={k: v for k, v in vars(controls).items()
                              if isinstance(v, (float, int))},
                  limitations=["Finite-grid Ghia table is retained; not an analytic exact solution",
                               "Affine fit reproduces linear fields; does not guarantee bounds for arbitrary data",
                               "This audit does not prove CFD engineering qualification"])
    try:
        if args.max_iterations <= 0:
            raise ValueError("max-iterations must be positive")
        seen = set()
        for path in args.summary:
            source = json.loads(path.read_text())
            report["sources"].append(dict(path=str(path.resolve()), sha256=sha256(path)))
            for old in source["cases"]:
                if old.get("case") != "cavity":
                    continue
                mesh, prefix = Path(old["mesh"]), Path(old["prefix"])
                if str(prefix.resolve()) in seen:
                    raise ValueError("Duplicate cavity artifacts")
                seen.add(str(prefix.resolve()))
                if sha256(mesh) != old["meshSha256"]:
                    raise ValueError(f"Mesh hash mismatch: {mesh}")
                declared = {str(Path(k).resolve()): v for k, v in old["artifactSha256"].items()}
                hashes = {"mesh": sha256(mesh)}
                for suffix in ("cells.csv", "faces.csv", "json", "residuals.csv"):
                    artifact = Path(str(prefix) + "." + suffix)
                    hashes[suffix] = sha256(artifact)
                    if hashes[suffix] != declared.get(str(artifact.resolve())):
                        raise ValueError(f"Artifact hash missing or changed: {artifact}")
                result = verifier.verify_case(mesh, prefix, "cavity", old["nu"], old["speed"], controls)
                result["label"] = old["label"]
                benchmark = result["benchmark"]
                report["cases"].append(dict(label=old["label"], scheme=result["native"]["convection"],
                    mesh=str(mesh), prefix=str(prefix), hashes=hashes, verification=result,
                    newRmse=benchmark["centrelineRmse"], legacyRmse=benchmark["legacyIdw"]["centrelineRmse"]))
        if not report["cases"]:
            raise ValueError("No cavity fields in source summaries")
        for scheme in sorted({c["scheme"] for c in report["cases"]}):
            cases = [c["verification"] for c in report["cases"] if c["scheme"] == scheme]
            legacy = copy.deepcopy(cases)
            for case in legacy:
                case["benchmark"] = case["benchmark"]["legacyIdw"]
            report["sequenceChecks"][scheme] = dict(current=verifier.sequence_checks(cases),
                                                    legacy=verifier.sequence_checks(legacy))
        report["calibration"] = calibration()
        report["valid"] = (all(c["verification"]["valid"] for c in report["cases"])
                           and all(s["current"]["valid"] for s in report["sequenceChecks"].values()))
        if not report["valid"]:
            report["issues"].append("Existing cavity case or refinement gate remains failed")
    except (OSError, ValueError, KeyError) as exc:
        report["issues"].append(str(exc))
    verifier.write_json(args.output, report)
    print(json.dumps({k: report[k] for k in ("valid", "issues", "sequenceChecks")}, indent=2))
    return 0 if report["valid"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
