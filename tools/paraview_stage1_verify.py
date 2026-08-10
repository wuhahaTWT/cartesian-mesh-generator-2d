#!/usr/bin/env pvpython
"""用 ParaView/VTK 独立读取 STL/VTU，验证相交单元并生成叠加图和分类切片。"""

from __future__ import annotations

import argparse
import hashlib
import json
import pathlib

import paraview
from vtkmodules.vtkCommonCore import vtkLookupTable, vtkVersion
from vtkmodules.vtkCommonDataModel import vtkDataObject, vtkPlane
from vtkmodules.vtkFiltersCore import vtkCutter, vtkThreshold
from vtkmodules.vtkFiltersGeometry import vtkGeometryFilter
from vtkmodules.vtkIOGeometry import vtkSTLReader
from vtkmodules.vtkIOImage import vtkPNGWriter
from vtkmodules.vtkIOXML import vtkXMLUnstructuredGridReader
from vtkmodules.vtkRenderingCore import (
    vtkActor,
    vtkDataSetMapper,
    vtkPolyDataMapper,
    vtkRenderWindow,
    vtkRenderer,
    vtkWindowToImageFilter,
)
import vtkmodules.vtkRenderingOpenGL2  # noqa: F401 - 注册 OpenGL 后端


def sha256(path: pathlib.Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_screenshot(window: vtkRenderWindow, path: pathlib.Path) -> dict[str, object]:
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
        raise RuntimeError(f"未能生成截图：{path}")
    return {"path": str(path), "bytes": path.stat().st_size, "sha256": sha256(path)}


def surface_actor(reader: vtkSTLReader) -> vtkActor:
    mapper = vtkPolyDataMapper()
    mapper.SetInputConnection(reader.GetOutputPort())
    actor = vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetColor(0.82, 0.84, 0.88)
    actor.GetProperty().SetOpacity(0.38)
    actor.GetProperty().EdgeVisibilityOn()
    actor.GetProperty().SetEdgeColor(0.18, 0.20, 0.24)
    return actor


def configure_window(renderer: vtkRenderer) -> vtkRenderWindow:
    renderer.SetBackground(0.96, 0.97, 0.99)
    window = vtkRenderWindow()
    window.SetOffScreenRendering(1)
    window.SetSize(1100, 850)
    window.AddRenderer(renderer)
    return window


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--surface", required=True, type=pathlib.Path)
    parser.add_argument("--mesh", required=True, type=pathlib.Path)
    parser.add_argument("--overview", required=True, type=pathlib.Path)
    parser.add_argument("--slice", required=True, type=pathlib.Path)
    parser.add_argument("--report", required=True, type=pathlib.Path)
    parser.add_argument("--slice-origin", type=float, nargs=3, required=True)
    parser.add_argument("--slice-normal", type=float, nargs=3, default=(0.0, 0.0, 1.0))
    arguments = parser.parse_args()

    mesh_reader = vtkXMLUnstructuredGridReader()
    mesh_reader.SetFileName(str(arguments.mesh.resolve()))
    mesh_reader.Update()
    grid = mesh_reader.GetOutput()
    if mesh_reader.GetErrorCode() != 0 or grid.GetNumberOfCells() == 0:
        raise RuntimeError(f"VTU 读取失败，错误码={mesh_reader.GetErrorCode()}")
    if any(grid.GetCellType(index) != 12 for index in range(grid.GetNumberOfCells())):
        raise RuntimeError("VTU 包含非 VTK_HEXAHEDRON 单元")

    field_name = "stl_cell_classification"
    field = grid.GetCellData().GetArray(field_name)
    if field is None or field.GetNumberOfTuples() != grid.GetNumberOfCells():
        raise RuntimeError(f"缺少完整的 {field_name} 字段")
    counts = [0, 0, 0, 0]
    for index in range(field.GetNumberOfTuples()):
        value = field.GetTuple1(index)
        rounded = int(round(value))
        if value != rounded or rounded < 0 or rounded > 3:
            raise RuntimeError(f"分类字段在单元 {index} 含非法值 {value}")
        counts[rounded] += 1

    stl_reader = vtkSTLReader()
    stl_reader.SetFileName(str(arguments.surface.resolve()))
    stl_reader.Update()
    surface = stl_reader.GetOutput()
    if stl_reader.GetErrorCode() != 0 or surface.GetNumberOfCells() == 0:
        raise RuntimeError(f"STL 读取失败，错误码={stl_reader.GetErrorCode()}")

    threshold = vtkThreshold()
    threshold.SetInputConnection(mesh_reader.GetOutputPort())
    threshold.SetInputArrayToProcess(
        0, 0, 0, vtkDataObject.FIELD_ASSOCIATION_CELLS, field_name
    )
    threshold.SetLowerThreshold(2.0)
    threshold.SetUpperThreshold(2.0)
    threshold.SetThresholdFunction(vtkThreshold.THRESHOLD_BETWEEN)
    threshold.Update()
    intersected_count = threshold.GetOutput().GetNumberOfCells()
    if intersected_count != counts[2]:
        raise RuntimeError("ParaView 精确阈值选中数与分类字段计数不一致")

    intersected_surface = vtkGeometryFilter()
    intersected_surface.SetInputConnection(threshold.GetOutputPort())
    intersected_mapper = vtkPolyDataMapper()
    intersected_mapper.SetInputConnection(intersected_surface.GetOutputPort())
    intersected_actor = vtkActor()
    intersected_actor.SetMapper(intersected_mapper)
    intersected_actor.GetProperty().SetColor(0.08, 0.58, 0.90)
    intersected_actor.GetProperty().SetOpacity(0.70)
    intersected_actor.GetProperty().EdgeVisibilityOn()
    intersected_actor.GetProperty().SetEdgeColor(0.03, 0.14, 0.24)

    overview_renderer = vtkRenderer()
    overview_renderer.AddActor(intersected_actor)
    overview_renderer.AddActor(surface_actor(stl_reader))
    overview_window = configure_window(overview_renderer)
    overview_renderer.ResetCamera()
    overview_renderer.GetActiveCamera().Azimuth(35)
    overview_renderer.GetActiveCamera().Elevation(25)
    overview_renderer.GetActiveCamera().Zoom(1.12)
    overview_result = write_screenshot(overview_window, arguments.overview)

    plane = vtkPlane()
    plane.SetOrigin(*arguments.slice_origin)
    plane.SetNormal(*arguments.slice_normal)
    cutter = vtkCutter()
    cutter.SetCutFunction(plane)
    cutter.SetInputConnection(mesh_reader.GetOutputPort())
    cutter.GenerateTrianglesOff()
    cutter.Update()
    slice_output = cutter.GetOutput()
    if slice_output.GetNumberOfCells() == 0:
        raise RuntimeError("指定平面没有产生分类切片")

    lookup = vtkLookupTable()
    lookup.SetNumberOfTableValues(4)
    lookup.SetTableValue(0, 0.88, 0.91, 0.95, 1.0)
    lookup.SetTableValue(1, 0.18, 0.72, 0.42, 1.0)
    lookup.SetTableValue(2, 0.08, 0.58, 0.90, 1.0)
    lookup.SetTableValue(3, 0.88, 0.16, 0.18, 1.0)
    lookup.SetTableRange(0.0, 3.0)
    lookup.Build()
    slice_mapper = vtkDataSetMapper()
    slice_mapper.SetInputConnection(cutter.GetOutputPort())
    slice_mapper.SetScalarModeToUseCellFieldData()
    slice_mapper.SelectColorArray(field_name)
    slice_mapper.SetLookupTable(lookup)
    slice_mapper.SetScalarRange(0.0, 3.0)
    slice_actor = vtkActor()
    slice_actor.SetMapper(slice_mapper)
    slice_actor.GetProperty().EdgeVisibilityOn()
    slice_actor.GetProperty().SetEdgeColor(0.12, 0.14, 0.18)

    slice_renderer = vtkRenderer()
    slice_renderer.AddActor(slice_actor)
    slice_renderer.AddActor(surface_actor(stl_reader))
    slice_window = configure_window(slice_renderer)
    slice_renderer.ResetCamera()
    slice_result = write_screenshot(slice_window, arguments.slice)

    result = {
        "reader": "ParaView/VTK STL and XML readers",
        "paraviewVersion": paraview.__version__,
        "vtkVersion": vtkVersion.GetVTKVersion(),
        "surfaceTriangleCount": surface.GetNumberOfCells(),
        "pointCount": grid.GetNumberOfPoints(),
        "cellCount": grid.GetNumberOfCells(),
        "cellType": "VTK_HEXAHEDRON",
        "field": field_name,
        "classificationCounts": {
            "outside": counts[0],
            "inside": counts[1],
            "intersected": counts[2],
            "conflict": counts[3],
        },
        "thresholdValue": 2,
        "thresholdSelectedCellCount": intersected_count,
        "sliceOrigin": arguments.slice_origin,
        "sliceNormal": arguments.slice_normal,
        "sliceCellCount": slice_output.GetNumberOfCells(),
        "sliceCellTypes": sorted(
            {slice_output.GetCellType(index) for index in range(slice_output.GetNumberOfCells())}
        ),
        "overview": overview_result,
        "slice": slice_result,
        "status": "pass",
    }
    serialized = json.dumps(result, ensure_ascii=False, indent=2, sort_keys=True)
    arguments.report.parent.mkdir(parents=True, exist_ok=True)
    arguments.report.write_text(serialized + "\n", encoding="utf-8")
    print(serialized)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
