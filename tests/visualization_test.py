#!/usr/bin/env python3
from __future__ import annotations

import json
import subprocess
import sys
import tempfile
from pathlib import Path


def main() -> int:
    repo_root = Path(__file__).resolve().parents[2]
    renderer = repo_root / "cartmesh2d" / "tools" / "visualization" / "render_cm2d.py"
    with tempfile.TemporaryDirectory(prefix="cartmesh2d_v_") as tmp:
        root = Path(tmp)
        cm2d = root / "fixture.cm2d"
        quality = root / "fixture.quality.json"
        viz = root / "fixture.viz.json"
        svg = root / "fixture.svg"
        cm2d.write_text(
            """CM2D 1
VERTICES 6
0 0 0
1 1 0
2 2 0
3 2 1
4 1 1
5 0 1
EDGES 7
0 0 1 0 -1 1
1 1 4 0 1 0
2 4 5 0 -1 1
3 5 0 0 -1 1
4 1 2 1 -1 1
5 2 3 1 -1 1
6 3 4 1 -1 1
CELLS 2
0 10 3 1 4 0 1 4 5 4 0 1 2 3
1 20 4 1 4 1 2 3 4 4 4 5 6 1
AUDIT 0 0 0 0 0 0 0
END
""",
            encoding="utf-8",
        )
        quality.write_text(
            json.dumps(
                {
                    "valid": True,
                    "counts": {"vertices": 6, "edges": 7, "cells": 2,
                               "source_cut_cells": 1, "source_small_cells": 1},
                    "quality": {"min_cell_area": 1.0, "min_edge_length": 1.0,
                                "max_cell_edge_length_ratio": 1.0,
                                "max_centroid_vertex_mean_offset_normalized": 0.0},
                    "topology_audit": {
                        "duplicate_vertices": 0, "duplicate_edges": 0,
                        "orphan_internal_edges": 0, "non_manifold_edges": 0,
                        "unclassified_boundary_edges": 0, "open_cell_loops": 0,
                        "area_mismatches": 0,
                    },
                },
                indent=2,
            ),
            encoding="utf-8",
        )
        viz.write_text(
            json.dumps(
                {
                    "format": "cartmesh2d-viz-v1",
                    "small_alpha_threshold": 0.2,
                    "source_small_cell_count": 1,
                    "merged_small_cell_count": 1,
                    "source_cells": [
                        {"source_id": 10, "source_key": 3, "level": 3, "kind": "cut",
                         "area_fraction": 0.1, "small": True, "small_status": "candidate",
                         "source_topology_cell_id": 0, "target_topology_cell_id": 1,
                         "background_bounds": [0.0, 0.0, 1.0, 1.0], "centroid": [0.5, 0.5]},
                        {"source_id": 20, "source_key": 4, "level": 4, "kind": "full",
                         "area_fraction": 1.0, "small": False, "small_status": "stable",
                         "source_topology_cell_id": 1, "target_topology_cell_id": None,
                         "background_bounds": [1.0, 0.0, 2.0, 1.0], "centroid": [1.5, 0.5]},
                    ],
                },
                indent=2,
            ),
            encoding="utf-8",
        )
        subprocess.run(
            [sys.executable, str(renderer), str(cm2d), str(svg),
             "--quality", str(quality), "--viz", str(viz), "--labels"],
            check=True,
        )
        text = svg.read_text(encoding="utf-8")
        required = [
            "<svg", 'id="source-background"', 'id="cells"', 'id="edges"',
            'class="edge-internal"', 'class="edge-embedded"', 'class="source-small"',
            'class="small-dot"', 'data-cell="0"', 'data-cell="1"', 'data-kind="cut"',
            "Quality summary", "duplicate_edges: 0", "viz small: 1",
            "cartmesh2d — stabilized solver topology",
        ]
        missing = [token for token in required if token not in text]
        if missing:
            raise RuntimeError(f"visualization SVG missing expected content: {missing}")

        bad = root / "bad.cm2d"
        bad_svg = root / "bad.svg"
        bad.write_text(cm2d.read_text(encoding="utf-8").replace(
            "AUDIT 0 0 0 0 0 0 0", "AUDIT 0 1 0 0 0 0 0"), encoding="utf-8")
        subprocess.run([sys.executable, str(renderer), str(bad), str(bad_svg)], check=True)
        if "INVALID TOPOLOGY / QUALITY" not in bad_svg.read_text(encoding="utf-8"):
            raise RuntimeError("renderer must visibly flag non-zero topology audit")

    print("cartmesh2d 2D-V visualization tests passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
