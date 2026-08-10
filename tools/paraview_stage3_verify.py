#!/usr/bin/env pvpython
"""ParaView/VTK 独立读取并渲染阶段三 Cut-cell 与嵌入边界。"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import paraview
from vtkmodules.vtkCommonCore import vtkLookupTable, vtkVersion
from vtkmodules.vtkCommonDataModel import vtkDataObject, vtkPlane
from vtkmodules.vtkFiltersCore import vtkCutter, vtkThreshold
from vtkmodules.vtkFiltersGeometry import vtkGeometryFilter
from vtkmodules.vtkFiltersGeneral import vtkCellValidator
from vtkmodules.vtkIOImage import vtkPNGWriter
from vtkmodules.vtkIOXML import vtkXMLPolyDataReader, vtkXMLUnstructuredGridReader
from vtkmodules.vtkRenderingCore import (
    vtkActor,
    vtkDataSetMapper,
    vtkPolyDataMapper,
    vtkRenderWindow,
    vtkRenderer,
    vtkWindowToImageFilter,
)
import vtkmodules.vtkRenderingOpenGL2  # noqa: F401


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def screenshot(window: vtkRenderWindow, path: Path) -> dict[str, object]:
    path.parent.mkdir(parents=True, exist_ok=True)
    window.Render()
    capture = vtkWindowToImageFilter()
    capture.SetInput(window)
    capture.SetInputBufferTypeToRGBA()
    capture.ReadFrontBufferOff()
    capture.Update()
    writer = vtkPNGWriter()
    writer.SetFileName(str(path.resolve()))
    writer.SetInputConnection(capture.GetOutputPort())
    writer.Write()
    if not path.exists() or path.stat().st_size == 0:
        raise RuntimeError(f"未生成截图：{path}")
    return {"path": str(path), "bytes": path.stat().st_size, "sha256": sha256(path)}


def window(renderer: vtkRenderer) -> vtkRenderWindow:
    renderer.SetBackground(0.96, 0.97, 0.99)
    result = vtkRenderWindow()
    result.SetOffScreenRendering(1)
    result.SetSize(1100, 850)
    result.AddRenderer(renderer)
    return result


def boundary_actor(reader: vtkXMLPolyDataReader) -> vtkActor:
    mapper = vtkPolyDataMapper()
    mapper.SetInputConnection(reader.GetOutputPort())
    actor = vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetColor(0.90, 0.20, 0.12)
    actor.GetProperty().SetOpacity(0.88)
    actor.GetProperty().EdgeVisibilityOn()
    actor.GetProperty().SetEdgeColor(0.30, 0.04, 0.02)
    return actor


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mesh", required=True, type=Path)
    parser.add_argument("--boundary", required=True, type=Path)
    parser.add_argument("--polyhedra", required=True, type=Path)
    parser.add_argument("--overview", required=True, type=Path)
    parser.add_argument("--slice", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    parser.add_argument("--slice-origin", required=True, nargs=3, type=float)
    parser.add_argument("--slice-normal", nargs=3, type=float, default=(0.0, 0.0, 1.0))
    arguments = parser.parse_args()

    mesh_reader = vtkXMLUnstructuredGridReader()
    mesh_reader.SetFileName(str(arguments.mesh.resolve()))
    mesh_reader.Update()
    grid = mesh_reader.GetOutput()
    if mesh_reader.GetErrorCode() != 0 or grid.GetNumberOfCells() == 0:
        raise RuntimeError(f"阶段三 VTU 读取失败，错误码={mesh_reader.GetErrorCode()}")
    if any(grid.GetCellType(index) != 12 for index in range(grid.GetNumberOfCells())):
        raise RuntimeError("阶段三背景 VTU 包含非 VTK_HEXAHEDRON 单元")
    fraction = grid.GetCellData().GetArray("fluid_volume_fraction")
    cut_field = grid.GetCellData().GetArray("cut_cell")
    if fraction is None or cut_field is None:
        raise RuntimeError("VTU 缺少 fluid_volume_fraction 或 cut_cell")
    cut_count = sum(cut_field.GetTuple1(index) == 1.0 for index in range(grid.GetNumberOfCells()))
    fraction_range = fraction.GetRange()

    boundary_reader = vtkXMLPolyDataReader()
    boundary_reader.SetFileName(str(arguments.boundary.resolve()))
    boundary_reader.Update()
    boundary = boundary_reader.GetOutput()
    if boundary_reader.GetErrorCode() != 0 or boundary.GetNumberOfCells() == 0:
        raise RuntimeError(f"嵌入边界 VTP 读取失败，错误码={boundary_reader.GetErrorCode()}")
    area_field = boundary.GetCellData().GetArray("area")
    normal_field = boundary.GetCellData().GetArray("fluid_outward_normal")
    if area_field is None or normal_field is None:
        raise RuntimeError("嵌入边界 VTP 缺少 area 或 fluid_outward_normal")
    embedded_area = sum(area_field.GetTuple1(index) for index in range(area_field.GetNumberOfTuples()))

    polyhedron_reader = vtkXMLUnstructuredGridReader()
    polyhedron_reader.SetFileName(str(arguments.polyhedra.resolve()))
    polyhedron_reader.Update()
    polyhedra = polyhedron_reader.GetOutput()
    if polyhedron_reader.GetErrorCode() != 0 or polyhedra.GetNumberOfCells() == 0:
        raise RuntimeError(
            f"显式流体 polyhedron VTU 读取失败，错误码={polyhedron_reader.GetErrorCode()}"
        )
    if any(polyhedra.GetCellType(index) != 42 for index in range(polyhedra.GetNumberOfCells())):
        raise RuntimeError("显式流体分解包含非 VTK_POLYHEDRON 单元")
    piece_volume_field = polyhedra.GetCellData().GetArray("piece_volume")
    if piece_volume_field is None:
        raise RuntimeError("显式流体 polyhedron 缺少 piece_volume")
    polyhedron_volume = sum(
        piece_volume_field.GetTuple1(index)
        for index in range(piece_volume_field.GetNumberOfTuples())
    )
    validator = vtkCellValidator()
    validator.SetInputConnection(polyhedron_reader.GetOutputPort())
    validator.Update()
    validity = validator.GetOutput().GetCellData().GetArray("ValidityState")
    if validity is None:
        raise RuntimeError("vtkCellValidator 未生成 ValidityState")
    invalid_polyhedra = sum(
        validity.GetTuple1(index) != 0.0 for index in range(validity.GetNumberOfTuples())
    )
    if invalid_polyhedra:
        raise RuntimeError(f"vtkCellValidator 报告 {invalid_polyhedra} 个无效 polyhedron")

    threshold = vtkThreshold()
    threshold.SetInputConnection(mesh_reader.GetOutputPort())
    threshold.SetInputArrayToProcess(
        0, 0, 0, vtkDataObject.FIELD_ASSOCIATION_CELLS, "cut_cell"
    )
    threshold.SetLowerThreshold(1.0)
    threshold.SetUpperThreshold(1.0)
    threshold.SetThresholdFunction(vtkThreshold.THRESHOLD_BETWEEN)
    threshold.Update()
    if threshold.GetOutput().GetNumberOfCells() != cut_count:
        raise RuntimeError("ParaView Cut-cell 阈值数量与字段计数不一致")
    cut_surface = vtkGeometryFilter()
    cut_surface.SetInputConnection(threshold.GetOutputPort())
    cut_mapper = vtkPolyDataMapper()
    cut_mapper.SetInputConnection(cut_surface.GetOutputPort())
    cut_actor = vtkActor()
    cut_actor.SetMapper(cut_mapper)
    cut_actor.GetProperty().SetColor(0.10, 0.55, 0.92)
    cut_actor.GetProperty().SetOpacity(0.32)
    cut_actor.GetProperty().EdgeVisibilityOn()
    cut_actor.GetProperty().SetEdgeColor(0.03, 0.15, 0.28)

    overview_renderer = vtkRenderer()
    overview_renderer.AddActor(cut_actor)
    overview_renderer.AddActor(boundary_actor(boundary_reader))
    polyhedron_surface = vtkGeometryFilter()
    polyhedron_surface.SetInputConnection(polyhedron_reader.GetOutputPort())
    polyhedron_mapper = vtkPolyDataMapper()
    polyhedron_mapper.SetInputConnection(polyhedron_surface.GetOutputPort())
    polyhedron_actor = vtkActor()
    polyhedron_actor.SetMapper(polyhedron_mapper)
    polyhedron_actor.GetProperty().SetColor(0.18, 0.72, 0.42)
    polyhedron_actor.GetProperty().SetOpacity(0.24)
    overview_renderer.AddActor(polyhedron_actor)
    overview_window = window(overview_renderer)
    overview_renderer.ResetCamera()
    overview_renderer.GetActiveCamera().Azimuth(35)
    overview_renderer.GetActiveCamera().Elevation(24)
    overview_renderer.GetActiveCamera().Zoom(1.12)
    overview = screenshot(overview_window, arguments.overview)

    plane = vtkPlane()
    plane.SetOrigin(*arguments.slice_origin)
    plane.SetNormal(*arguments.slice_normal)
    cutter = vtkCutter()
    cutter.SetCutFunction(plane)
    cutter.SetInputConnection(mesh_reader.GetOutputPort())
    cutter.GenerateTrianglesOff()
    cutter.Update()
    slice_data = cutter.GetOutput()
    if slice_data.GetNumberOfCells() == 0:
        raise RuntimeError("指定切片平面未穿过背景网格")
    lookup = vtkLookupTable()
    lookup.SetNumberOfTableValues(5)
    lookup.SetTableValue(0, 0.25, 0.25, 0.28, 1.0)
    lookup.SetTableValue(1, 0.16, 0.44, 0.82, 1.0)
    lookup.SetTableValue(2, 0.10, 0.66, 0.78, 1.0)
    lookup.SetTableValue(3, 0.20, 0.75, 0.45, 1.0)
    lookup.SetTableValue(4, 0.92, 0.86, 0.20, 1.0)
    lookup.SetTableRange(0.0, 1.0)
    lookup.Build()
    slice_mapper = vtkDataSetMapper()
    slice_mapper.SetInputConnection(cutter.GetOutputPort())
    slice_mapper.SetScalarModeToUseCellFieldData()
    slice_mapper.SelectColorArray("fluid_volume_fraction")
    slice_mapper.SetLookupTable(lookup)
    slice_mapper.SetScalarRange(0.0, 1.0)
    slice_actor = vtkActor()
    slice_actor.SetMapper(slice_mapper)
    slice_actor.GetProperty().EdgeVisibilityOn()
    slice_actor.GetProperty().SetEdgeColor(0.10, 0.12, 0.16)
    slice_renderer = vtkRenderer()
    slice_renderer.AddActor(slice_actor)
    slice_renderer.AddActor(boundary_actor(boundary_reader))
    slice_window = window(slice_renderer)
    slice_renderer.ResetCamera()
    slice = screenshot(slice_window, arguments.slice)

    result = {
        "schema": "cartmesh-stage3-paraview-verification-v1",
        "status": "pass",
        "geometryTopologyValidated": True,
        "solverReadyCutCellMesh": False,
        "externalCfdCheckerAccepted": False,
        "paraviewVersion": paraview.__version__,
        "vtkVersion": vtkVersion.GetVTKVersion(),
        "meshReaderErrorCode": mesh_reader.GetErrorCode(),
        "boundaryReaderErrorCode": boundary_reader.GetErrorCode(),
        "polyhedronReaderErrorCode": polyhedron_reader.GetErrorCode(),
        "backgroundPointCount": grid.GetNumberOfPoints(),
        "backgroundCellCount": grid.GetNumberOfCells(),
        "cutCellCount": cut_count,
        "thresholdSelectedCellCount": threshold.GetOutput().GetNumberOfCells(),
        "fluidVolumeFractionRange": list(fraction_range),
        "embeddedPolygonCount": boundary.GetNumberOfCells(),
        "embeddedPointCount": boundary.GetNumberOfPoints(),
        "embeddedBoundaryArea": embedded_area,
        "fluidPolyhedronPieceCount": polyhedra.GetNumberOfCells(),
        "fluidPolyhedronPieceVolume": polyhedron_volume,
        "vtkCellValidatorInvalidCount": invalid_polyhedra,
        "sliceCellCount": slice_data.GetNumberOfCells(),
        "sliceOrigin": arguments.slice_origin,
        "sliceNormal": arguments.slice_normal,
        "overview": overview,
        "slice": slice,
    }
    serialized = json.dumps(result, ensure_ascii=False, indent=2, sort_keys=True)
    arguments.report.parent.mkdir(parents=True, exist_ok=True)
    arguments.report.write_text(serialized + "\n", encoding="utf-8")
    print(serialized)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
