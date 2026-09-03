#!/bin/sh
set -eu

SCRIPT_DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
DESKTOP_DIR=$(CDPATH= cd -- "$SCRIPT_DIR/.." && pwd)
PROJECT_DIR=$(CDPATH= cd -- "$DESKTOP_DIR/.." && pwd)
BUILD_DIR=${CARTMESH2D_BUILD_DIR:-"$PROJECT_DIR/build-desktop"}

cmake -S "$PROJECT_DIR" -B "$BUILD_DIR" -DCMAKE_BUILD_TYPE=Release
cmake --build "$BUILD_DIR" --target cartmesh2d_cli cartmesh2d_dxf_cli cartmesh2d_hybrid_cli -j 8

mkdir -p "$DESKTOP_DIR/runtime/bin" "$DESKTOP_DIR/runtime/examples"
cp "$BUILD_DIR/cartmesh2d_cli" "$DESKTOP_DIR/runtime/bin/cartmesh2d_cli"
cp "$BUILD_DIR/cartmesh2d_dxf_cli" "$DESKTOP_DIR/runtime/bin/cartmesh2d_dxf_cli"
cp "$BUILD_DIR/cartmesh2d_hybrid_cli" "$DESKTOP_DIR/runtime/bin/cartmesh2d_hybrid_cli"
cp "$PROJECT_DIR/examples/dxf/spline_circle_mm.dxf" "$DESKTOP_DIR/runtime/examples/spline_circle_mm.dxf"

cd "$DESKTOP_DIR"
npm run pack:mac
