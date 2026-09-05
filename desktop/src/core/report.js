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

const CELL_TYPE_LABELS = {
  cartesian: 'Cartesian',
  remainder_cut: 'Cut-cell 余域',
  transition: '过渡环',
  termination: '层终止',
  boundary_layer: '边界层'
};

function summarizeContract(contract) {
  if (!contract) return null;
  const byType = Object.entries(contract.by_cell_type || {}).map(([type, row]) => ({
    type,
    label: CELL_TYPE_LABELS[type] || type,
    status: row.rated ? row.status : 'OBSERVED',
    rated: Boolean(row.rated),
    cellCount: row.cell_count || 0,
    hard: row.hard_issue_count || 0,
    preferred: row.preferred_issue_count || 0
  }));
  const metrics = {};
  for (const issue of contract.issues || []) {
    const key = `${issue.level}:${issue.metric}`;
    metrics[key] = (metrics[key] || 0) + 1;
  }
  return { status: contract.status, byType, metricCounts: metrics, limits: contract.contract };
}

const SOLVER_METRIC_LABELS = [
  ['max_non_orthogonality_deg', '最大非正交角', '°', 70, 'max'],
  ['max_internal_skewness', '最大内部歪斜', '', 4, 'max'],
  ['max_boundary_skewness', '最大边界歪斜', '', 4, 'max'],
  ['min_face_weight', '最小面权重', '', 0.05, 'min'],
  ['min_volume_ratio', '最小体积比', '', 0.01, 'min'],
  ['max_cell_aspect', '最大单元长宽比', '', 1000, 'max']
];

function summarizeSolverQuality(metrics, valid) {
  if (!metrics) return null;
  return {
    valid: Boolean(valid),
    rows: SOLVER_METRIC_LABELS
      .filter(([key]) => Number.isFinite(Number(metrics[key])))
      .map(([key, label, unit, limit, direction]) => ({
        key, label, unit, limit, direction,
        value: Number(metrics[key]),
        pass: direction === 'max' ? Number(metrics[key]) <= limit : Number(metrics[key]) >= limit
      }))
  };
}

// The four gates are independent; passing one says nothing about another, so they are
// reported side by side rather than collapsed into a single verdict.
function normalizeResult({ method, stdout, reports, paths, mesh }) {
  const values = parseKeyValues(stdout);
  const contract = summarizeContract(reports.contract);
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

  const solverMetrics = isHybrid
    ? (reports.solverQuality ? reports.solverQuality.metrics : null)
    : (reports.contract ? reports.contract.legacy_hard_safety : null);
  const solverValid = isHybrid
    ? (reports.solverQuality ? reports.solverQuality.valid : undefined)
    : (reports.contract ? reports.contract.legacy_hard_safety?.valid : undefined);

  const written = values.openfoam === 'written' || Boolean(values.openfoam_case);
  return {
    requestedMethod: method,
    actualMethod: fellBack ? 'cutcell-fallback' : method,
    counts,
    gates: {
      // Both CLIs return EXIT_FAILURE when the topology audit fails, so reaching a
      // zero exit *is* the audit passing.  normalizeResult is only called on success.
      topology: { pass: true, label: '拓扑不变量' },
      solver: summarizeSolverQuality(solverMetrics, solverValid),
      contract
    },
    sizeField: reports.sizeField || null,
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

module.exports = { parseKeyValues, summarizeContract, summarizeSolverQuality, normalizeResult, CELL_TYPE_LABELS };
