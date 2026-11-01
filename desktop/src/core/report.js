'use strict';

// Turn the two CLIs' very different output into one shape the renderer can draw.
//
// The pure CLI prints one `key=value` per line; the hybrid CLI prints many pairs per
// line.  Neither is the authoritative source: the JSON reports are, and stdout is
// used for the live log and for status words only.

function parseKeyValues(text) {
  const values = {};
  for (const line of text.split(/\r?\n/)) {
    const whole = line.match(/^\s*([a-z0-9_]+)=(.*)$/i);
    // A line holding a single pair may have a value containing spaces (a path), so
    // it is taken verbatim.  Only split lines that really do carry several pairs.
    if (whole && (line.match(/\b[a-z0-9_]+=/gi) || []).length === 1) {
      values[whole[1]] = whole[2].trim();
      continue;
    }
    for (const pair of line.matchAll(/([a-z0-9_]+)=(\S+)/gi)) values[pair[1]] = pair[2];
  }
  return values;
}

const SOLVER_METRIC_LABELS = [
  ['max_non_orthogonality_deg', '最大非正交角', '°', 'max'],
  ['max_internal_skewness', '最大内部歪斜', '', 'max'],
  ['max_boundary_skewness', '最大边界歪斜', '', 'max'],
  ['min_face_weight', '最小面权重', '', 'min'],
  ['min_volume_ratio', '最小邻接面积比', '', 'min'],
  ['max_cell_aspect', '最大原生形状比', '', 'max']
];

function summarizeSolverQuality(metrics, valid, report = {}) {
  if (!metrics) return null;
  const policy = report.policy || {};
  return {
    valid: Boolean(valid), advisory: true,
    nativeFvCompatible: report.native_fv_compatible ?? null,
    rows: SOLVER_METRIC_LABELS
      .filter(([key]) => metrics[key] !== null && Number.isFinite(Number(metrics[key])))
      .map(([key, label, unit, direction]) => {
        const limit = Number.isFinite(policy[key]) ? policy[key] : null;
        const value = Number(metrics[key]);
        return { key, label, unit, limit, direction, value,
          pass: limit === null ? null : direction === 'max' ? value <= limit : value >= limit };
      })
  };
}

// The topology and Solver gates are independent, so they are
// reported side by side rather than collapsed into a single verdict.
function normalizeResult({ method, stdout, reports, paths, mesh, incomplete = false }) {
  const values = parseKeyValues(stdout);
  const hybrid = reports.hybrid || null;
  const isHybrid = method === 'hybrid';
  const fellBack = isHybrid && values.mesh_mode === 'pure_cutcell_fallback';

  // The mesh file is the fallback rather than a nicety: a run that dies during the
  // OpenFOAM export never printed its summary block, but it did write the mesh.
  const counts = {
    cells: Number(values.solver_output_cells || values.solver_cells || values.stabilized_cells ||
      (hybrid && hybrid.cell_count) || (mesh && mesh.cells.length) || 0),
    vertices: Number(values.solver_output_vertices || values.vertices ||
      (mesh && mesh.vertices.length) || 0),
    faces: Number(values.solver_output_faces || values.edges || (mesh && mesh.edges.length) || 0),
    leaves: Number(values.leaf_count || 0),
    cutCells: Number(values.cut_cells || values.remainder_cut || 0),
    layerCells: Number(values.layer_cells || (hybrid && hybrid.boundary_layer_cell_count) || 0)
  };

  const solverMetrics = reports.solverQuality ? reports.solverQuality.metrics : null;
  const solverValid = reports.solverQuality ? reports.solverQuality.valid : undefined;

  const written = values.openfoam === 'written' || Boolean(values.openfoam_case);
  return {
    requestedMethod: method,
    actualMethod: fellBack ? 'cutcell-fallback' : method,
    counts,
    gates: {
      // Both CLIs return EXIT_FAILURE when the topology audit fails, so reaching a
      // zero exit *is* the internal audit passing. Partial runs remain unconfirmed.
      topology: { pass: incomplete ? null : true, label: '拓扑不变量' },
      solver: summarizeSolverQuality(solverMetrics, solverValid, reports.solverQuality || {})
    },
    externalCheck: reports.externalCheck || null,
    sizeField: reports.sizeField || null,
    resolution: reports.resolution || null,
    sizing: reports.sizing || null,
    openFoam: {
      written,
      path: paths.casePath,
      // The hybrid CLI reports `openfoam=written` without a count; the exported cell
      // count is the solver mesh's own.
      cells: Number(values.openfoam_cells || (written ? counts.cells : 0))
    },
    timings: Object.fromEntries(Object.entries(values)
      .filter(([key]) => key.startsWith('timing_'))
      .map(([key, value]) => [key.replace('timing_', ''), Number(value)])),
    raw: values
  };
}

module.exports = { parseKeyValues, summarizeSolverQuality, normalizeResult };
