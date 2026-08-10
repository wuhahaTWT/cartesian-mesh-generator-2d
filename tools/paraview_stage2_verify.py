#!/usr/bin/env pvpython
"""Read adaptive VTU/STL with ParaView, verify fields, and render level views."""

from __future__ import annotations

import argparse
import hashlib
import json
import pathlib

import paraview
from vtkmodules.vtkCommonCore import vtkLookupTable, vtkVersion
from vtkmodules.vtkCommonDataModel import vtkDataObject, vtkPlane
from vtkmodules.vtkFiltersCore import vtkCutter, vtkThreshold
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
import vtkmodules.vtkRenderingOpenGL2  # noqa: F401


def sha256(path: pathlib.Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def screenshot(renderer: vtkRenderer, path: pathlib.Path, size: tuple[int, int]) -> dict[str, object]:
    window = vtkRenderWindow()
    window.SetOffScreenRendering(1)
    window.SetSize(*size)
    window.AddRenderer(renderer)
    renderer.SetBackground(0.96, 0.97, 0.99)
    window.Render()
    capture = vtkWindowToImageFilter()
    capture.SetInput(window)
    capture.SetInputBufferTypeToRGBA()
    capture.ReadFrontBufferOff()
    capture.Update()
    path.parent.mkdir(parents=True, exist_ok=True)
    writer = vtkPNGWriter()
    writer.SetFileName(str(path.resolve()))
    writer.SetInputConnection(capture.GetOutputPort())
    writer.Write()
    if not path.exists() or path.stat().st_size == 0:
        raise RuntimeError(f"failed to render {path}")
    return {"path": str(path), "bytes": path.stat().st_size, "sha256": sha256(path)}


def surface_actor(reader: vtkSTLReader) -> vtkActor:
    mapper = vtkPolyDataMapper()
    mapper.SetInputConnection(reader.GetOutputPort())
    actor = vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetColor(0.15, 0.17, 0.20)
    actor.GetProperty().SetOpacity(0.28)
    actor.GetProperty().EdgeVisibilityOn()
    actor.GetProperty().SetEdgeColor(0.05, 0.05, 0.05)
    return actor


def level_lookup(minimum: int, maximum: int) -> vtkLookupTable:
    lookup = vtkLookupTable()
    lookup.SetNumberOfTableValues(maximum - minimum + 1)
    palette = [
        (0.19, 0.40, 0.72, 1.0),
        (0.12, 0.66, 0.78, 1.0),
        (0.25, 0.72, 0.45, 1.0),
        (0.94, 0.73, 0.20, 1.0),
        (0.91, 0.34, 0.20, 1.0),
        (0.64, 0.20, 0.58, 1.0),
    ]
    for index in range(maximum - minimum + 1):
        lookup.SetTableValue(index, *palette[index % len(palette)])
    lookup.SetTableRange(float(minimum), float(maximum))
    lookup.Build()
    return lookup


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--surface", required=True, type=pathlib.Path)
    parser.add_argument("--mesh", required=True, type=pathlib.Path)
    parser.add_argument("--overview", required=True, type=pathlib.Path)
    parser.add_argument("--slice", required=True, type=pathlib.Path)
    parser.add_argument("--report", required=True, type=pathlib.Path)
    parser.add_argument("--slice-origin", required=True, type=float, nargs=3)
    parser.add_argument("--slice-normal", type=float, nargs=3, default=(0.0, 0.0, 1.0))
    arguments = parser.parse_args()

    mesh_reader = vtkXMLUnstructuredGridReader()
    mesh_reader.SetFileName(str(arguments.mesh.resolve()))
    mesh_reader.Update()
    grid = mesh_reader.GetOutput()
    if mesh_reader.GetErrorCode() or grid.GetNumberOfCells() == 0:
        raise RuntimeError(f"VTU reader error {mesh_reader.GetErrorCode()}")
    if any(grid.GetCellType(index) != 12 for index in range(grid.GetNumberOfCells())):
        raise RuntimeError("adaptive VTU contains non-hexahedron cells")
    level_name = "octree_level"
    classification_name = "stl_cell_classification"
    level_array = grid.GetCellData().GetArray(level_name)
    classification_array = grid.GetCellData().GetArray(classification_name)
    if level_array is None or classification_array is None:
        raise RuntimeError("adaptive VTU lacks level or classification cell field")
    levels: dict[int, int] = {}
    classifications = [0, 0, 0, 0]
    for cell_id in range(grid.GetNumberOfCells()):
        level_value = level_array.GetTuple1(cell_id)
        level = int(round(level_value))
        if level_value != level or level < 0 or level > 21:
            raise RuntimeError(f"invalid octree level at cell {cell_id}")
        levels[level] = levels.get(level, 0) + 1
        classification_value = classification_array.GetTuple1(cell_id)
        classification = int(round(classification_value))
        if classification_value != classification or classification not in range(4):
            raise RuntimeError(f"invalid classification at cell {cell_id}")
        classifications[classification] += 1

    threshold = vtkThreshold()
    threshold.SetInputConnection(mesh_reader.GetOutputPort())
    threshold.SetInputArrayToProcess(
        0, 0, 0, vtkDataObject.FIELD_ASSOCIATION_CELLS, classification_name
    )
    threshold.SetLowerThreshold(2.0)
    threshold.SetUpperThreshold(2.0)
    threshold.SetThresholdFunction(vtkThreshold.THRESHOLD_BETWEEN)
    threshold.Update()
    if threshold.GetOutput().GetNumberOfCells() != classifications[2]:
        raise RuntimeError("exact intersected threshold count mismatch")

    stl_reader = vtkSTLReader()
    stl_reader.SetFileName(str(arguments.surface.resolve()))
    stl_reader.Update()
    if stl_reader.GetErrorCode() or stl_reader.GetOutput().GetNumberOfCells() == 0:
        raise RuntimeError(f"STL reader error {stl_reader.GetErrorCode()}")

    minimum_level = min(levels)
    maximum_level = max(levels)
    lookup = level_lookup(minimum_level, maximum_level)
    overview_mapper = vtkDataSetMapper()
    overview_mapper.SetInputConnection(mesh_reader.GetOutputPort())
    overview_mapper.SetScalarModeToUseCellFieldData()
    overview_mapper.SelectColorArray(level_name)
    overview_mapper.SetLookupTable(lookup)
    overview_mapper.SetScalarRange(minimum_level, maximum_level)
    overview_actor = vtkActor()
    overview_actor.SetMapper(overview_mapper)
    overview_actor.GetProperty().SetOpacity(0.72)
    overview_actor.GetProperty().EdgeVisibilityOn()
    overview_actor.GetProperty().SetEdgeColor(0.10, 0.11, 0.14)
    overview_renderer = vtkRenderer()
    overview_renderer.AddActor(overview_actor)
    overview_renderer.AddActor(surface_actor(stl_reader))
    overview_renderer.ResetCamera()
    overview_renderer.GetActiveCamera().Azimuth(32)
    overview_renderer.GetActiveCamera().Elevation(22)
    overview_result = screenshot(overview_renderer, arguments.overview, (1150, 850))

    plane = vtkPlane()
    plane.SetOrigin(*arguments.slice_origin)
    plane.SetNormal(*arguments.slice_normal)
    cutter = vtkCutter()
    cutter.SetCutFunction(plane)
    cutter.SetInputConnection(mesh_reader.GetOutputPort())
    cutter.GenerateTrianglesOff()
    cutter.Update()
    if cutter.GetOutput().GetNumberOfCells() == 0:
        raise RuntimeError("requested adaptive level slice is empty")
    slice_mapper = vtkDataSetMapper()
    slice_mapper.SetInputConnection(cutter.GetOutputPort())
    slice_mapper.SetScalarModeToUseCellFieldData()
    slice_mapper.SelectColorArray(level_name)
    slice_mapper.SetLookupTable(lookup)
    slice_mapper.SetScalarRange(minimum_level, maximum_level)
    slice_actor = vtkActor()
    slice_actor.SetMapper(slice_mapper)
    slice_actor.GetProperty().EdgeVisibilityOn()
    slice_actor.GetProperty().SetEdgeColor(0.08, 0.09, 0.12)
    slice_renderer = vtkRenderer()
    slice_renderer.AddActor(slice_actor)
    slice_renderer.AddActor(surface_actor(stl_reader))
    slice_renderer.ResetCamera()
    slice_result = screenshot(slice_renderer, arguments.slice, (1150, 850))

    result = {
        "reader": "ParaView/VTK STL and XML readers",
        "paraviewVersion": paraview.__version__,
        "vtkVersion": vtkVersion.GetVTKVersion(),
        "surfaceTriangleCount": stl_reader.GetOutput().GetNumberOfCells(),
        "pointCount": grid.GetNumberOfPoints(),
        "leafCount": grid.GetNumberOfCells(),
        "levelField": level_name,
        "leafCountByLevel": {str(level): count for level, count in sorted(levels.items())},
        "classificationCounts": {
            "outside": classifications[0],
            "inside": classifications[1],
            "intersected": classifications[2],
            "conflict": classifications[3],
        },
        "intersectedThresholdSelectedCellCount": threshold.GetOutput().GetNumberOfCells(),
        "sliceOrigin": arguments.slice_origin,
        "sliceNormal": arguments.slice_normal,
        "sliceCellCount": cutter.GetOutput().GetNumberOfCells(),
        "overview": overview_result,
        "slice": slice_result,
        "status": "pass",
    }
    arguments.report.parent.mkdir(parents=True, exist_ok=True)
    arguments.report.write_text(
        json.dumps(result, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(result, ensure_ascii=False, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
