#!/usr/bin/env pvpython
"""通过 ParaView 的 VTK 运行时读取 VTU 并生成离屏 PNG。"""

from __future__ import annotations

import argparse
import json
import pathlib

import paraview
from vtkmodules.vtkCommonCore import vtkVersion
from vtkmodules.vtkCommonDataModel import vtkDataObject
from vtkmodules.vtkFiltersCore import vtkThreshold
from vtkmodules.vtkFiltersGeometry import vtkGeometryFilter
from vtkmodules.vtkIOXML import vtkXMLUnstructuredGridReader
from vtkmodules.vtkRenderingCore import (
    vtkActor,
    vtkPolyDataMapper,
    vtkRenderWindow,
    vtkRenderer,
    vtkWindowToImageFilter,
)
from vtkmodules.vtkIOImage import vtkPNGWriter
import vtkmodules.vtkRenderingOpenGL2  # noqa: F401 - 注册 OpenGL 后端


def main() -> int:
    parser = argparse.ArgumentParser(description="使用 ParaView 读取并离屏渲染阶段 0 VTU")
    parser.add_argument("mesh", type=pathlib.Path)
    parser.add_argument("screenshot", type=pathlib.Path)
    parser.add_argument("--field", help="渲染前把指定单元字段筛选到 [0.5, 1.0]")
    parser.add_argument("--report", type=pathlib.Path, help="把 JSON 结果写入指定路径")
    arguments = parser.parse_args()

    reader = vtkXMLUnstructuredGridReader()
    reader.SetFileName(str(arguments.mesh.resolve()))
    reader.Update()
    grid = reader.GetOutput()
    if reader.GetErrorCode() != 0:
        raise RuntimeError(f"ParaView/VTK 读取器错误码：{reader.GetErrorCode()}")
    if grid.GetNumberOfPoints() == 0 or grid.GetNumberOfCells() == 0:
        raise RuntimeError("ParaView/VTK 读取器返回了空网格")
    if any(grid.GetCellType(cell) != 12 for cell in range(grid.GetNumberOfCells())):
        raise RuntimeError("ParaView/VTK 读取器发现了非六面体的阶段 0 单元")

    pipeline = reader
    selected_cell_count = grid.GetNumberOfCells()
    if arguments.field:
        if grid.GetCellData().GetArray(arguments.field) is None:
            raise RuntimeError(f"缺少请求的单元字段：{arguments.field}")
        threshold = vtkThreshold()
        threshold.SetInputConnection(reader.GetOutputPort())
        threshold.SetInputArrayToProcess(
            0, 0, 0, vtkDataObject.FIELD_ASSOCIATION_CELLS, arguments.field
        )
        threshold.SetLowerThreshold(0.5)
        threshold.SetUpperThreshold(1.0)
        threshold.SetThresholdFunction(vtkThreshold.THRESHOLD_BETWEEN)
        threshold.Update()
        selected_cell_count = threshold.GetOutput().GetNumberOfCells()
        if selected_cell_count == 0:
            raise RuntimeError("指定阈值没有选中任何单元")
        pipeline = threshold

    surface = vtkGeometryFilter()
    surface.SetInputConnection(pipeline.GetOutputPort())
    mapper = vtkPolyDataMapper()
    mapper.SetInputConnection(surface.GetOutputPort())
    actor = vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetColor(0.2, 0.65, 0.95)
    actor.GetProperty().EdgeVisibilityOn()
    actor.GetProperty().SetEdgeColor(0.08, 0.12, 0.18)

    renderer = vtkRenderer()
    renderer.SetBackground(0.96, 0.97, 0.99)
    renderer.AddActor(actor)
    window = vtkRenderWindow()
    window.SetOffScreenRendering(1)
    window.SetSize(1000, 800)
    window.AddRenderer(renderer)
    renderer.ResetCamera()
    renderer.GetActiveCamera().Azimuth(35)
    renderer.GetActiveCamera().Elevation(25)
    renderer.GetActiveCamera().Zoom(1.15)
    window.Render()

    arguments.screenshot.parent.mkdir(parents=True, exist_ok=True)
    capture = vtkWindowToImageFilter()
    capture.SetInput(window)
    capture.SetInputBufferTypeToRGBA()
    capture.ReadFrontBufferOff()
    capture.Update()
    writer = vtkPNGWriter()
    writer.SetFileName(str(arguments.screenshot.resolve()))
    writer.SetInputConnection(capture.GetOutputPort())
    writer.Write()
    if not arguments.screenshot.exists() or arguments.screenshot.stat().st_size == 0:
        raise RuntimeError("ParaView/VTK 未能生成截图")

    result = {
        "reader": "ParaView vtkXMLUnstructuredGridReader",
        "paraviewVersion": paraview.__version__,
        "vtkVersion": vtkVersion.GetVTKVersion(),
        "pointCount": grid.GetNumberOfPoints(),
        "cellCount": grid.GetNumberOfCells(),
        "selectedCellCount": selected_cell_count,
        "cellType": "VTK_HEXAHEDRON",
        "cellDataFields": [
            grid.GetCellData().GetArrayName(index)
            for index in range(grid.GetCellData().GetNumberOfArrays())
        ],
        "screenshot": str(arguments.screenshot),
        "status": "pass",
    }
    serialized = json.dumps(result, sort_keys=True)
    if arguments.report:
        arguments.report.parent.mkdir(parents=True, exist_ok=True)
        arguments.report.write_text(serialized + "\n", encoding="utf-8")
    print(serialized)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
