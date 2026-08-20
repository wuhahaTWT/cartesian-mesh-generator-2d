# Complex geometry stress check — 2026-08-20

This is an exploratory robustness check beyond the formal Stage 2D-0~V acceptance set. It does not redefine the already-closed acceptance gate.

GitHub Actions run #24 (`32336795285`) kept the formal Stage 0~6 suite green (15/15) and then exercised harder single-loop geometries.

## Results

- `nozzle_profile.xy` — PASS at maxLevel 7. End-to-end CLI, stabilization, topology audit and SVG rendering all passed.
  - leaf_count = 1480
  - source_cells = 904
  - small_cells = 4
  - stabilized_cells = 900
  - vertices = 1273
  - edges = 2172
  - min_cell_area = 3.17186e-05
  - min_edge_length = 1.31417e-05
  - max_edge_aspect_ratio = 4993.64
  - max_centroid_skewness = 0.4222
  - all topology audit counters = 0
- `gear_star.xy` — FAIL. 220 unsupported Cut-cell leaf cases at maxLevel 8. Increasing refinement from the earlier level-6 run did not remove the issue.
- `serpentine_body.xy` — FAIL. 1 unsupported Cut-cell leaf case at maxLevel 7.
- `naca2412_dense.xy` — FAIL. Cut-cell construction reached topology assembly, but source global topology audit failed at maxLevel 8.
- `superellipse_24.xy` — later smoke attempt also reached Cut-cell construction but failed source global topology audit at maxLevel 7.

## Interpretation

The first native-2D product pipeline is complete for its formal acceptance envelope, but it is not yet robust for arbitrary complex single-loop geometry. The main next robustness work should target:

1. concave polygon / AABB clipping that can return multiple connected components without fabricating bridges;
2. explicit multi-component per-background-cell representation or deterministic refinement-until-single-component policy;
3. topology tolerance/canonicalization around curved or nearly tangent boundary intersections;
4. sharp/thin trailing-edge handling;
5. a permanent complex-geometry stress suite kept separate from the core acceptance suite.

Do not describe the current implementation as supporting arbitrary complex 2D geometry until these stress cases are closed.
