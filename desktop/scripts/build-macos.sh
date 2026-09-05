#!/bin/sh
set -eu

SCRIPT_DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
DESKTOP_DIR=$(CDPATH= cd -- "$SCRIPT_DIR/.." && pwd)
PROJECT_DIR=$(CDPATH= cd -- "$DESKTOP_DIR/.." && pwd)
BUILD_DIR=${CARTMESH2D_BUILD_DIR:-"$PROJECT_DIR/build-desktop"}

# Pin the system toolchain.  mesasdk puts its own c++ on PATH, and cmake picks that up
# by default; the resulting binary loads @rpath/lib/libstdc++.6.dylib, which only
# resolves with DYLD_LIBRARY_PATH set, so the packaged .app dies on launch.
CXX_COMPILER=${CARTMESH2D_CXX:-/usr/bin/clang++}
cmake -S "$PROJECT_DIR" -B "$BUILD_DIR" -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_CXX_COMPILER="$CXX_COMPILER"
cmake --build "$BUILD_DIR" --target cartmesh2d_cli cartmesh2d_hybrid_cli cartmesh2d_dxf_cli -j 8

rm -rf "$DESKTOP_DIR/runtime"
mkdir -p "$DESKTOP_DIR/runtime/bin" "$DESKTOP_DIR/runtime/samples"
for tool in cartmesh2d_cli cartmesh2d_hybrid_cli cartmesh2d_dxf_cli; do
  cp "$BUILD_DIR/$tool" "$DESKTOP_DIR/runtime/bin/$tool"
  # A bundled tool must run with no environment help at all.
  if otool -L "$DESKTOP_DIR/runtime/bin/$tool" | grep -q '@rpath/lib/libstdc++'; then
    echo "$tool links a non-system libstdc++; it would not run inside the .app" >&2
    exit 1
  fi
done

# Sample files named by desktop/src/core/samples.js.  Copied flat because the app
# resolves them as resources/samples/<file>.
for sample in \
  acceptance/circle.xy \
  complex/naca2412_dense.xy \
  complex/thick_cambered_airfoil.xy \
  complex/two_obstacles.xy \
  complex/superellipse_24.xy \
  complex/serpentine_body.xy \
  complex/gear_star.xy \
  complex/annulus.xy \
  complex/nozzle_profile.xy \
  h4_3/narrow_gap.xy \
  h4_3/sharp_trailing_edge.xy
do
  cp "$PROJECT_DIR/examples/$sample" "$DESKTOP_DIR/runtime/samples/$(basename "$sample")"
done

# Fail loudly here rather than shipping a sample the picker cannot open.
node -e '
const { SAMPLES } = require("./src/core/samples");
const fs = require("node:fs");
const missing = SAMPLES.filter(s => !fs.existsSync(`runtime/samples/${s.file}`));
if (missing.length) {
  console.error("missing sample files: " + missing.map(s => s.file).join(", "));
  process.exit(1);
}
console.log(`bundled ${SAMPLES.length} samples`);
'

cd "$DESKTOP_DIR"
npm test
npm run pack:mac
