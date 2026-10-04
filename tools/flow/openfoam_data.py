"""OpenFOAM field and boundary readers for visualization."""
from pathlib import Path
import re
import math

def payload(path: Path) -> str:
    text = path.read_text(encoding="utf-8")
    at = text.find("FoamFile")
    end = text.find("}", at)
    if at < 0 or end < 0:
        raise ValueError(f"{path}: missing FoamFile header")
    return text[end + 1:]


def counted(path: Path) -> tuple[int, str]:
    match = re.search(r"\b(\d+)\s*\n\s*\((.*)\)\s*$", payload(path), re.S)
    if not match:
        raise ValueError(f"{path}: malformed counted list")
    return int(match.group(1)), match.group(2)


def read_points(path: Path):
    count, body = counted(path)
    rows = [tuple(map(float, row.split())) for row in re.findall(r"\(([^()]*)\)", body)]
    if len(rows) != count or any(len(row) != 3 for row in rows):
        raise ValueError(f"{path}: point count/vector mismatch")
    return rows


def read_faces(path: Path):
    count, body = counted(path)
    faces = []
    for declared, row in re.findall(r"(\d+)\s*\(([^()]*)\)", body):
        face = [int(value) for value in row.split()]
        if len(face) != int(declared):
            raise ValueError(f"{path}: face arity mismatch")
        faces.append(face)
    if len(faces) != count:
        raise ValueError(f"{path}: face count mismatch")
    return faces


def read_labels(path: Path):
    count, body = counted(path)
    values = [int(value) for value in re.findall(r"\b-?\d+\b", body)]
    if len(values) != count:
        raise ValueError(f"{path}: label count mismatch")
    return values


def read_boundary(path: Path):
    count, body = counted(path)
    patches = []
    for name, fields in re.findall(r"([A-Za-z_][A-Za-z0-9_]*)\s*\{(.*?)\}", body, re.S):
        values = dict(re.findall(r"(type|nFaces|startFace)\s+([^;]+);", fields))
        if set(values) != {"type", "nFaces", "startFace"}:
            raise ValueError(f"{path}: incomplete patch {name}")
        patches.append({"name": name, "type": values["type"].strip(),
                        "nFaces": int(values["nFaces"]),
                        "startFace": int(values["startFace"])})
    if len(patches) != count:
        raise ValueError(f"{path}: patch count mismatch")
    return patches


def header(cls: str, location: str, obj: str) -> str:
    return ("FoamFile\n{\n    version 2.0;\n    format ascii;\n"
            f"    class {cls};\n    location \"{location}\";\n    object {obj};\n}}\n\n")


def scalar_field(path: Path):
    text = payload(path)
    uniform = re.search(r"internalField\s+uniform\s+([^;\s]+)", text)
    if uniform:
        value = float(uniform.group(1))
        if not math.isfinite(value): raise ValueError(f"{path}: non-finite scalar")
        return [value]
    match = re.search(r"internalField\s+nonuniform\s+List<scalar>\s+(\d+)\s*\((.*?)\)\s*;", text, re.S)
    if not match:
        raise ValueError(f"{path}: unsupported scalar field")
    values = [float(v) for v in match.group(2).split()]
    if len(values) != int(match.group(1)) or any(not math.isfinite(v) for v in values):
        raise ValueError(f"{path}: scalar count or finiteness mismatch")
    return values


def vector_field(path: Path):
    text = payload(path)
    uniform = re.search(r"internalField\s+uniform\s*\(([^()]*)\)", text)
    if uniform:
        row = tuple(float(v) for v in uniform.group(1).split())
        if len(row) != 3 or any(not math.isfinite(v) for v in row): raise ValueError(f"{path}: invalid internal vector")
        return [row]
    match = re.search(r"internalField\s+nonuniform\s+List<vector>\s+(\d+)\s*\((.*?)\)\s*;", text, re.S)
    if not match: raise ValueError(f"{path}: unsupported vector field")
    rows = [tuple(float(v) for v in row.split()) for row in re.findall(r"\(([^()]*)\)", match.group(2))]
    if len(rows) != int(match.group(1)) or any(len(row) != 3 or any(not math.isfinite(v) for v in row) for row in rows):
        raise ValueError(f"{path}: vector count or finiteness mismatch")
    return rows


def boundary_type(path: Path, name: str):
    text = payload(path); at = text.find("boundaryField")
    match = re.search(rf"\b{re.escape(name)}\s*\{{(.*?)\n\s*\}}", text[at:], re.S)
    if not match: raise ValueError(f"{path}: missing boundary {name}")
    kind = re.search(r"\btype\s+([^;\s]+)", match.group(1))
    return kind.group(1) if kind else None


def boundary_values(path: Path, vector: bool = False):
    """Read solved boundary values, preserving uniform/nonuniform cardinality."""
    text = payload(path); at = text.find("boundaryField")
    if at < 0:
        raise ValueError(f"{path}: missing boundaryField")
    values = {}
    for name, block in re.findall(r"\n\s*([A-Za-z_][A-Za-z0-9_]*)\s*\{(.*?)\n\s*\}", text[at:], re.S):
        if vector:
            uni = re.search(r"\bvalue\s+uniform\s*\(([^()]*)\)", block)
            if uni:
                row = tuple(float(v) for v in uni.group(1).split())
                if len(row) != 3 or any(not math.isfinite(v) for v in row): raise ValueError(f"{path}: invalid vector value for {name}")
                values[name] = [row]
                continue
            non = re.search(r"\bvalue\s+nonuniform\s+List<vector>\s+(\d+)\s*\((.*?)\)\s*;", block, re.S)
            if non:
                rows = [tuple(float(v) for v in row.split()) for row in re.findall(r"\(([^()]*)\)", non.group(2))]
                if len(rows) != int(non.group(1)) or any(len(row) != 3 or any(not math.isfinite(v) for v in row) for row in rows):
                    raise ValueError(f"{path}: vector value count mismatch for {name}")
                values[name] = rows
        else:
            uni = re.search(r"\bvalue\s+uniform\s+([^;\s]+)", block)
            if uni:
                value = float(uni.group(1))
                if not math.isfinite(value): raise ValueError(f"{path}: non-finite scalar value for {name}")
                values[name] = [value]
                continue
            non = re.search(r"\bvalue\s+nonuniform\s+List<scalar>\s+(\d+)\s*\((.*?)\)\s*;", block, re.S)
            if non:
                rows = [float(v) for v in non.group(2).split()]
                if len(rows) != int(non.group(1)) or any(not math.isfinite(v) for v in rows): raise ValueError(f"{path}: scalar value count mismatch for {name}")
                values[name] = rows
    return values


def face_area(face, points):
    origin = points[face[0]]; area = [0.0, 0.0, 0.0]
    for i in range(1, len(face) - 1):
        a = points[face[i]]; b = points[face[i + 1]]
        u = [a[j] - origin[j] for j in range(3)]; v = [b[j] - origin[j] for j in range(3)]
        cross = [u[1] * v[2] - u[2] * v[1], u[2] * v[0] - u[0] * v[2], u[0] * v[1] - u[1] * v[0]]
        for j in range(3): area[j] += 0.5 * cross[j]
    return area

