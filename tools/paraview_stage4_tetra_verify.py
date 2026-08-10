#!/usr/bin/env pvpython
"""VTK/ParaView 独立验证阶段四的四面体显式分解。"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import paraview
from vtkmodules.vtkCommonCore import vtkVersion
from vtkmodules.vtkFiltersGeneral import vtkCellValidator
from vtkmodules.vtkIOXML import vtkXMLUnstructuredGridReader


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tetrahedra", required=True, type=Path)
    parser.add_argument("--project-report", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    arguments = parser.parse_args()

    project = json.loads(arguments.project_report.read_text(encoding="utf-8"))
    reader = vtkXMLUnstructuredGridReader()
    reader.SetFileName(str(arguments.tetrahedra.resolve()))
    reader.Update()
    grid = reader.GetOutput()
    if reader.GetErrorCode() != 0 or grid.GetNumberOfCells() == 0:
        raise RuntimeError(f"tetra VTU 读取失败，错误码={reader.GetErrorCode()}")
    non_tetra = sum(grid.GetCellType(cell) != 10
                    for cell in range(grid.GetNumberOfCells()))
    volume = grid.GetCellData().GetArray("tetra_volume")
    source_piece = grid.GetCellData().GetArray("source_piece_id")
    region = grid.GetCellData().GetArray("global_region_id")
    if volume is None or source_piece is None or region is None:
        raise RuntimeError("tetra VTU 缺少体积、原始片或 region 字段")
    volume_sum = sum(volume.GetTuple1(cell)
                     for cell in range(grid.GetNumberOfCells()))
    source_piece_count = len({int(source_piece.GetTuple1(cell))
                              for cell in range(grid.GetNumberOfCells())})
    validator = vtkCellValidator()
    validator.SetInputConnection(reader.GetOutputPort())
    validator.Update()
    states = validator.GetOutput().GetCellData().GetArray("ValidityState")
    invalid = sum(states.GetTuple1(cell) != 0
                  for cell in range(states.GetNumberOfTuples()))
    expected_volume = float(project["explicitFluidPieceVolume"])
    volume_error = abs(volume_sum - expected_volume)
    volume_tolerance = 2.0e-11 * max(1.0, abs(expected_volume))
    status = (non_tetra == 0 and invalid == 0 and
              source_piece_count == int(project["explicitFluidPieceCount"]) and
              volume_error <= volume_tolerance)
    result = {
        "schema": "cartmesh-stage4-paraview-tetra-verification-v1",
        "status": "pass" if status else "fail",
        "validatedScope": "cut_cell_fluid_polyhedron_pieces_only",
        "completeSolverVolumeMeshValidated": False,
        "paraviewVersion": paraview.__version__,
        "vtkVersion": vtkVersion.GetVTKVersion(),
        "readerErrorCode": reader.GetErrorCode(),
        "tetrahedronCount": grid.GetNumberOfCells(),
        "nonTetraCellCount": non_tetra,
        "vtkCellValidatorInvalidCount": invalid,
        "sourcePieceCount": source_piece_count,
        "expectedSourcePieceCount": int(project["explicitFluidPieceCount"]),
        "tetrahedronVolumeSum": volume_sum,
        "expectedFluidPieceVolume": expected_volume,
        "absoluteVolumeDifference": volume_error,
        "globalRegionIds": sorted({int(region.GetTuple1(cell))
                                   for cell in range(grid.GetNumberOfCells())}),
        "input": str(arguments.tetrahedra),
        "inputBytes": arguments.tetrahedra.stat().st_size,
        "inputSha256": sha256(arguments.tetrahedra),
    }
    arguments.report.parent.mkdir(parents=True, exist_ok=True)
    arguments.report.write_text(json.dumps(result, indent=2) + "\n",
                                encoding="utf-8")
    print(json.dumps(result, ensure_ascii=False, indent=2))
    if not status:
        raise RuntimeError("阶段四 tetra 外部验证失败")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
