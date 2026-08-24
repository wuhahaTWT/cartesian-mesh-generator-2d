#!/usr/bin/env python3
"""Configure and evaluate deterministic harmonic manufactured solutions.

The utility is deliberately independent of the cartmesh2d library.  It reads
OpenFOAM's own ``writeCellCentres``/``writeCellVolumes`` products, prescribes
Dirichlet data at the actual boundary-face centres, and evaluates the solved
cell field with volume-weighted norms.
"""

import argparse
import json
import math
import re
from pathlib import Path


def foam_payload(path: Path) -> str:
    text = path.read_text(encoding="utf-8")
    end = text.find("}", text.find("FoamFile"))
    if end < 0:
        raise ValueError(f"{path}: missing FoamFile header")
    return text[end + 1 :]


def scalar_internal(path: Path) -> list[float]:
    text = foam_payload(path)
    uniform = re.search(r"internalField\s+uniform\s+([^;\s]+)\s*;", text)
    if uniform:
        return [float(uniform.group(1))]
    match = re.search(
        r"internalField\s+nonuniform\s+List<scalar>\s+(\d+)\s*\((.*?)\)\s*;",
        text, re.S)
    if not match:
        raise ValueError(f"{path}: unsupported internalField")
    values = [float(value) for value in match.group(2).split()]
    if len(values) != int(match.group(1)):
        raise ValueError(f"{path}: internalField count mismatch")
    return values


def boundary_blocks(path: Path) -> dict[str, tuple[str, list[float]]]:
    text = foam_payload(path)
    start = text.find("boundaryField")
    if start < 0:
        raise ValueError(f"{path}: missing boundaryField")
    body = text[start:]
    blocks = {}
    for name, block in re.findall(r"\n\s*([A-Za-z_][A-Za-z0-9_]*)\s*\{(.*?)\n\s*\}", body, re.S):
        kind = re.search(r"\btype\s+([A-Za-z_][A-Za-z0-9_]*)\s*;", block)
        if not kind:
            continue
        if kind.group(1) == "empty":
            blocks[name] = ("empty", [])
            continue
        uniform = re.search(r"\bvalue\s+uniform\s+([^;\s]+)\s*;", block)
        if uniform:
            blocks[name] = ("uniform", [float(uniform.group(1))])
            continue
        values = re.search(
            r"\bvalue\s+nonuniform\s+List<scalar>\s+(\d+)\s*\((.*?)\)\s*;",
            block, re.S)
        if not values:
            raise ValueError(f"{path}: unsupported boundary values for {name}")
        parsed = [float(value) for value in values.group(2).split()]
        if len(parsed) != int(values.group(1)):
            raise ValueError(f"{path}: boundary count mismatch for {name}")
        blocks[name] = ("nonuniform", parsed)
    return blocks


def expanded(kind: str, values: list[float], count: int) -> list[float]:
    if kind == "uniform":
        return values * count
    if kind == "nonuniform" and len(values) == count:
        return values
    raise ValueError("coordinate patch representations do not agree")


def exact_value(x: float, y: float, field: str) -> float:
    if field == "constant":
        return 1.25
    if field == "linear":
        return 0.7 * x - 0.3 * y + 0.2
    if field == "quadratic":
        return x * x - y * y + 0.15 * x - 0.05 * y
    raise ValueError(f"unknown field {field}")


def header(object_name: str, location: str, field_class: str = "dictionary") -> str:
    return f"""FoamFile
{{
    version 2.0;
    format ascii;
    class {field_class};
    location \"{location}\";
    object {object_name};
}}
"""


def write_scalar_field(case: Path, field: str, initialise_exact: bool) -> None:
    cx = scalar_internal(case / "0" / "Cx")
    cy = scalar_internal(case / "0" / "Cy")
    if len(cx) != len(cy):
        raise ValueError("Cx/Cy cell counts differ")
    internal = ([exact_value(x, y, field) for x, y in zip(cx, cy)]
                if initialise_exact else [0.0] * len(cx))
    bx = boundary_blocks(case / "0" / "Cx")
    by = boundary_blocks(case / "0" / "Cy")
    if bx.keys() != by.keys():
        raise ValueError("Cx/Cy patch sets differ")
    lines = [header("T", "0", "volScalarField"),
             "dimensions [0 0 0 0 0 0 0];\n",
             f"internalField nonuniform List<scalar>\n{len(internal)}\n(\n"]
    lines.extend(f"{value:.17g}\n" for value in internal)
    lines.append(")\n;\n\nboundaryField\n{\n")
    for name in bx:
        xkind, xvalues = bx[name]
        ykind, yvalues = by[name]
        lines.append(f"    {name}\n    {{\n")
        if xkind == "empty" or ykind == "empty":
            if xkind != ykind:
                raise ValueError(f"empty patch mismatch for {name}")
            lines.append("        type empty;\n")
        else:
            count = max(len(xvalues), len(yvalues))
            xs = expanded(xkind, xvalues, count)
            ys = expanded(ykind, yvalues, count)
            values = [exact_value(x, y, field) for x, y in zip(xs, ys)]
            lines.append("        type fixedValue;\n")
            lines.append(f"        value nonuniform List<scalar>\n        {count}\n        (\n")
            lines.extend(f"        {value:.17g}\n" for value in values)
            lines.append("        )\n        ;\n")
        lines.append("    }\n")
    lines.append("}\n")
    (case / "0" / "T").write_text("".join(lines), encoding="utf-8")


def configure(case: Path, field: str, initialise_exact: bool) -> None:
    write_scalar_field(case, field, initialise_exact)
    (case / "constant" / "transportProperties").write_text(
        header("transportProperties", "constant") +
        "DT [0 2 -1 0 0 0 0] 1;\n", encoding="utf-8")
    (case / "system" / "controlDict").write_text(
        header("controlDict", "system") +
        "application laplacianFoam;\nstartFrom startTime;\nstartTime 0;\n"
        "stopAt endTime;\nendTime 30;\ndeltaT 1;\nwriteControl timeStep;\n"
        "writeInterval 30;\nwriteFormat ascii;\nwritePrecision 17;\n"
        "writeCompression off;\ntimeFormat general;\ntimePrecision 6;\n"
        "runTimeModifiable false;\n", encoding="utf-8")
    (case / "system" / "fvSchemes").write_text(
        header("fvSchemes", "system") +
        "ddtSchemes { default steadyState; }\n"
        "gradSchemes { default Gauss linear; }\n"
        "divSchemes { default none; }\n"
        "laplacianSchemes { default Gauss linear corrected; }\n"
        "interpolationSchemes { default linear; }\n"
        "snGradSchemes { default corrected; }\n"
        "fluxRequired { default no; T; }\n", encoding="utf-8")
    (case / "system" / "fvSolution").write_text(
        header("fvSolution", "system") +
        "solvers\n{\n    T\n    {\n        solver PCG;\n        preconditioner DIC;\n"
        "        tolerance 1e-13;\n        relTol 0;\n        maxIter 20000;\n    }\n}\n",
        encoding="utf-8")


def evaluate(case: Path, time_name: str, field: str, log: Path | None) -> dict:
    cx = scalar_internal(case / "0" / "Cx")
    cy = scalar_internal(case / "0" / "Cy")
    volumes = scalar_internal(case / "0" / "V")
    numerical = scalar_internal(case / time_name / "T")
    if not (len(cx) == len(cy) == len(volumes) == len(numerical)):
        raise ValueError("cell field counts differ")
    exact = [exact_value(x, y, field) for x, y in zip(cx, cy)]
    errors = [value - reference for value, reference in zip(numerical, exact)]
    total_volume = sum(volumes)
    l1 = sum(v * abs(e) for v, e in zip(volumes, errors)) / total_volume
    l2 = math.sqrt(sum(v * e * e for v, e in zip(volumes, errors)) / total_volume)
    linf = max(map(abs, errors), default=0.0)
    mean_error = sum(v * e for v, e in zip(volumes, errors)) / total_volume
    report = {"valid": all(math.isfinite(value) for value in numerical),
            "field": field, "time": time_name, "cell_count": len(cx),
            "total_volume": total_volume, "l1": l1, "l2": l2,
            "linf": linf, "volume_weighted_mean_error": mean_error,
            "min_numerical": min(numerical), "max_numerical": max(numerical)}
    if log:
        matches = re.findall(
            r"Solving for T, Initial residual = ([^,]+), Final residual = ([^,]+), No Iterations (\d+)",
            log.read_text(encoding="utf-8"))
        if not matches:
            raise ValueError(f"{log}: no T solver records")
        initial, final, iterations = matches[-1]
        report["linear_system_initial_residual"] = float(initial)
        report["linear_system_final_residual"] = float(final)
        report["linear_system_iterations"] = int(iterations)
        report["valid"] = report["valid"] and float(final) <= 1.0e-12
    return report


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=("configure", "evaluate"))
    parser.add_argument("case", type=Path)
    parser.add_argument("--field", choices=("constant", "linear", "quadratic"), required=True)
    parser.add_argument("--time", default="30")
    parser.add_argument("--initialise-exact", action="store_true")
    parser.add_argument("--log", type=Path)
    parser.add_argument("--report", type=Path)
    args = parser.parse_args()
    if args.mode == "configure":
        configure(args.case, args.field, args.initialise_exact)
        report = {"valid": True, "field": args.field, "configured": str(args.case)}
    else:
        report = evaluate(args.case, args.time, args.field, args.log)
    encoded = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(encoded, encoding="utf-8")
    print(encoded, end="")


if __name__ == "__main__":
    main()
