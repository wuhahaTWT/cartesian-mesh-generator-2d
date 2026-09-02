#!/bin/bash
# R2 tidy-up byte-identity gate.
#
# Runs a fixed set of old command lines against two builds and compares every
# mesh product byte for byte. This is the only evidence that a simplification
# removed nothing the generator actually used.
#
# usage: byte_identity.sh <baseline-repo> <candidate-repo>
set -uo pipefail
export DYLD_LIBRARY_PATH=/Applications/mesasdk/lib

# Resolve both repos to absolute paths before any cd: run_all cd's into each one,
# so a relative "." would silently keep pointing at the first repo and compare a
# build against itself.
BASE="$(cd "${1:?baseline repo}" && pwd)"
CAND="$(cd "${2:?candidate repo}" && pwd)"
OUT=/tmp/tidy_cmp
rm -rf "$OUT"; mkdir -p "$OUT/base" "$OUT/cand"

run_all() {
  local repo="$1" out="$2"
  cd "$repo" || return 1
  echo "running $repo" >&2
  ./build/cartmesh2d_cli examples/acceptance/circle.xy "$out/cut_l8" \
    8 0.25 0.10 exterior "$out/cut_l8-case" 6 > "$out/cut_l8.log" 2>&1
  ./build/cartmesh2d_cli examples/acceptance/airfoil_like.xy "$out/cut_airfoil" \
    7 0.30 0.02 exterior "$out/cut_airfoil-case" 4 > "$out/cut_airfoil.log" 2>&1
  ./build/cartmesh2d_cli examples/complex/two_obstacles.xy "$out/cut_two" \
    6 0.25 0.10 exterior "$out/cut_two-case" 4 > "$out/cut_two.log" 2>&1
  ./build/cartmesh2d_hybrid_cli examples/acceptance/circle.xy "$out/hy_circle" \
    6 3 6 4 0.02 1.2 1.0 "$out/hy_circle-case" 0.01 > "$out/hy_circle.log" 2>&1
  ./build/cartmesh2d_hybrid_cli examples/complex/superellipse_24.xy "$out/hy_super" \
    6 3 6 3 0.015 1.15 1.0 "$out/hy_super-case" 0.01 > "$out/hy_super.log" 2>&1
  ./build/cartmesh2d_hybrid_cli examples/h4_3/concave_l.xy "$out/hy_concave" \
    8 3 8 4 0.012 1.15 1.0 "$out/hy_concave-case" 0.01 > "$out/hy_concave.log" 2>&1
  ./build/cartmesh2d_hybrid_cli examples/h4_3/narrow_gap.xy "$out/hy_narrow" \
    8 3 8 4 0.012 1.15 1.0 "$out/hy_narrow-case" 0.01 > "$out/hy_narrow.log" 2>&1
  ./build/cartmesh2d_hybrid_cli examples/h4_3/sharp_trailing_edge.xy "$out/hy_sharp" \
    8 3 8 4 0.012 1.15 1.0 "$out/hy_sharp-case" 0.01 > "$out/hy_sharp.log" 2>&1
  ./build/cartmesh2d_hybrid_cli examples/h4_3/narrow_gap.xy "$out/hy_q33" \
    8 3 8 4 0.012 1.15 1.0 "$out/hy_q33-case" 0.01 \
    --q3-termination-grouped > "$out/hy_q33.log" 2>&1
  ./build/cartmesh2d_hybrid_cli examples/h4_3/narrow_gap.xy "$out/hy_q41" \
    8 3 8 4 0.012 1.15 1.0 "$out/hy_q41-case" 0.01 \
    --q4-termination-construction --q3-termination-grouped > "$out/hy_q41.log" 2>&1
  ./build/cartmesh2d_dxf_cli examples/dxf/airfoil_like.dxf \
    "$out/dxf.xy" 0.001 "$out/dxf.json" > "$out/dxf.log" 2>&1
  ./build/cartmesh2d_boundary_layer_cli examples/acceptance/circle.xy \
    "$out/bl_circle" 4 first 0.02 1.2 exterior > "$out/bl_circle.log" 2>&1
}

run_all "$BASE" "$OUT/base"
run_all "$CAND" "$OUT/cand"

# Timings are wall time and legitimately differ; the profile file is the only
# product allowed to. Everything else must match byte for byte.
IDENTICAL=0; DIFFERENT=0; MISSING=0
while IFS= read -r rel; do
  case "$rel" in
    *.log|*.profile.json) continue ;;
  esac
  a="$OUT/base/$rel"; b="$OUT/cand/$rel"
  if [ ! -f "$b" ]; then echo "MISSING   $rel"; MISSING=$((MISSING+1)); continue; fi
  if cmp -s "$a" "$b"; then IDENTICAL=$((IDENTICAL+1)); else echo "DIFFERS   $rel"; DIFFERENT=$((DIFFERENT+1)); fi
done < <(cd "$OUT/base" && find . -type f | sed 's|^\./||' | sort)

echo
echo "identical=$IDENTICAL differs=$DIFFERENT missing=$MISSING"
[ "$DIFFERENT" -eq 0 ] && [ "$MISSING" -eq 0 ]
