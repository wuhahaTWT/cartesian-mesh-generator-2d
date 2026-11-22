'use strict';

// Result cards, quality advisories, legend, size histogram, and SVG time-series
// monitors for the desktop renderer. Kept independent of app.js state so layout
// and chart rendering can evolve without bloating the main controller.
(function () {
  const $ = id => document.getElementById(id);
  const fmt = value => Number(value || 0).toLocaleString('en-US');

  const gateRow = (label, verdict, detail) =>
    `<div class="gate"><b>${label}</b><span class="verdict ${verdict}">${verdict}</span>` +
    `<details class="gate-detail"><summary>详情</summary><span class="detail">${detail}</span></details></div>`;

  function renderCounters(result) {
    if (result.background) {
      $('counters').innerHTML = [
        ['完整单元', result.counts.cells],
        ['外部', result.counts.classification[0]],
        ['内部', result.counts.classification[1]],
        ['相交（未裁切）', result.counts.classification[2]]
      ].map(([label, n]) => `<div><span>${label}</span><b>${fmt(n)}</b></div>`).join('');
      return;
    }
    const count = value => (value ? fmt(value) : '—');
    const rows = [
      ['solver cells', count(result.counts.cells)],
      ['vertices', count(result.counts.vertices)],
      ['faces', count(result.counts.faces)],
      ['quadtree leaves', count(result.counts.leaves)],
      ['cut cells', count(result.counts.cutCells)],
      ['openfoam cells', result.openFoam.written ? fmt(result.openFoam.cells) : '未导出']
    ];
    if (result.counts.layerCells) rows.splice(4, 0, ['layer cells', fmt(result.counts.layerCells)]);
    $('counters').innerHTML = rows
      .map(([label, value]) => `<div><span>${label}</span><b>${value}</b></div>`).join('');
  }

  function renderGates(result) {
    if (result.background) {
      $('gates').innerHTML = gateRow('背景网格', '已生成',
        '完整域覆盖与记录一致性已核对。物体内部保留；几何轮廓不是数值壁面，未生成流体求解拓扑或OpenFOAM算例。');
      return;
    }
    const parts = [];
    parts.push(gateRow('内部拓扑检查',
      result.gates.topology.pass === null ? '未确认' : result.gates.topology.pass ? 'PASS' : 'FAIL',
      result.gates.topology.pass ? '生成器内部检查通过；外部 checkMesh 需另外执行。' : '本次未完整成功，不能据此确认通过。'));
    const solver = result.gates.solver;
    if (solver) {
      const worst = solver.rows.filter(row => row.pass === false);
      parts.push(gateRow('网格形状建议', worst.length ? 'WARN' : 'INFO',
        (worst.length ? worst.map(row => `${row.label} ${row.value.toPrecision(4)}（参考 ${row.limit}）`).join('　') :
          solver.rows.map(row => `${row.label} ${row.value.toPrecision(4)}`).join('　')) +
        '。仅供优化参考，不阻止导出或计算；不是外部 checkMesh 结果。'));
      if (solver.nativeFvCompatible === false) {
        parts.push(gateRow('当前原生离散', 'WARN',
          '存在当前离散无法使用的面几何；网格仍可导出，由目标 CFD 软件检查和处理。'));
      }
    }
    const directional = result.resolution?.directional_connectivity;
    if (directional) {
      const minimum = Number.isFinite(directional.minimum) ? directional.minimum.toPrecision(5) : '未测得';
      parts.push(gateRow('挤出方向诊断建议', directional.valid ? 'INFO' : 'WARN',
        `最小 determinant ${minimum}，参考值 ${directional.threshold}；需关注单元 ${directional.failed_cell_ids.length}。` +
        '仅为本项目诊断，不阻止导出；外部 checkMesh 结果单独显示。'));
    }
    if (result.actualMethod === 'cutcell-fallback') {
      parts.push(gateRow('方法', 'WARN', '请求了贴体边界层，实际落到纯 Cut-cell fallback'));
    }
    if (!result.openFoam.written) {
      parts.push(gateRow('OpenFOAM 导出', 'WARN', '未写出。请查看实际拓扑、几何或文件导出错误；经验质量建议不阻止导出。'));
    }
    const external = result.externalCheck;
    parts.push(gateRow('OpenFOAM checkMesh',
      external ? (external.checks.some(check => check.outcome !== 'passed') ? 'WARN' : 'INFO') : '未运行',
      external ? external.checks.map(check =>
        `${check.name === 'standard' ? '标准检查' : '扩展检查'}：${({ passed: '工具未报告失败项', reported_issues: '工具报告问题', tool_error: '工具未完成' })[check.outcome] || check.outcome}，退出码 ${check.exit_code}`
      ).join('；') :
        '尚无该导出网格的外部检查结果；本项目建议不代表成熟工具检查。'));
    $('gates').innerHTML = parts.join('');
  }

  function renderHistogram(histogram, mesh, basis) {
    const { levelColour } = window.MeshView;
    const container = $('histogram');
    container.replaceChildren();
    if (!histogram.length) return;
    const title = document.createElement('div');
    title.className = 'title';
    title.textContent = basis === 'size' ? '每尺寸档单元数（粗→细）' : '每层级单元数';
    container.appendChild(title);
    const peak = Math.max(...histogram.map(row => row.count));
    for (const row of histogram) {
      const bar = document.createElement('div');
      bar.className = 'bar';
      const label = document.createElement('em');
      label.textContent = `${basis === 'size' ? '档' : 'L'}${row.level}`;
      const fill = document.createElement('i');
      fill.style.width = `${Math.max(1, Math.round((row.count / peak) * 108))}px`;
      fill.style.background = levelColour(row.level, mesh.minLevel, mesh.maxLevel);
      const count = document.createElement('b');
      count.textContent = fmt(row.count);
      bar.append(label, fill, count);
      container.appendChild(bar);
    }
  }

  function renderLegend(mesh, basis, view) {
    const { RAMP } = window.MeshView;
    const container = $('legend');
    container.replaceChildren();
    if (mesh.background && view.mode === 'classification') {
      for (const [i, label] of ['外部', '内部', '相交'].entries()) {
        const entry = document.createElement('span');
        entry.textContent = `${label} ${fmt(mesh.classificationCounts[i])}　`;
        entry.style.color = ['#79b8d4', '#9da9bc', '#efac56'][i];
        container.appendChild(entry);
      }
      container.hidden = false;
      return;
    }
    const fieldMode = ['speed', 'pressure', 'temperature', 'adaptive'].includes(view.mode);
    const coloured = view.mode === 'level' || fieldMode;
    const ramp = document.createElement('div');
    ramp.className = 'ramp';
    const palette = fieldMode ? view.fieldPalette() : RAMP;
    for (const colour of palette) {
      const swatch = document.createElement('span');
      swatch.style.background = colour;
      ramp.appendChild(swatch);
    }
    const ends = document.createElement('div');
    ends.className = 'ends';
    const coarse = document.createElement('span');
    const fine = document.createElement('span');
    if (fieldMode) {
      const range = view.fieldRange;
      const unit = view.mode === 'adaptive' ? '无量纲指示子' : view.mode === 'temperature' ? 'K' : view.mode === 'speed' ? 'm/s' : 'm²/s²';
      coarse.textContent = range ? `${range.min.toPrecision(4)} ${unit}` : '最小';
      fine.textContent = range ? `${range.max.toPrecision(4)} ${unit}` : '最大';
    } else {
      const prefix = basis === 'size' ? '档' : 'L';
      coarse.textContent = `${prefix}${mesh.minLevel} 粗`;
      fine.textContent = `细 ${prefix}${mesh.maxLevel}`;
    }
    ends.append(coarse, fine);

    const keys = document.createElement('div');
    keys.className = 'keys';
    const theme = view.theme();
    for (const [colour, text] of [[theme.wall, '物面（嵌入边界）'], [theme.domain, '计算域边界']]) {
      const key = document.createElement('span');
      const dash = document.createElement('i');
      dash.style.background = colour;
      key.append(dash, document.createTextNode(text));
      keys.appendChild(key);
    }
    if (coloured) container.append(ramp, ends);
    container.appendChild(keys);
    container.hidden = false;
  }

  function renderResolution(resolution) {
    const el = $('resolutionResult');
    el.hidden = !resolution;
    if (!resolution) return;
    const actual = resolution.actual;
    const exceed = resolution.wall_owner_tangential_exceedance_length_fraction;
    el.textContent = `最终求解网格实测：Lref=${resolution.reference_length.toPrecision(5)} m；` +
      `单元等效尺寸 h/Lref 中位 ${actual.sqrt_area_over_reference?.p50.toPrecision(4) ?? '—'}；` +
      `壁面单元切向跨度 P95 ${actual.wall_owner_tangential_extent_over_reference?.p95.toPrecision(4) ?? '—'}；` +
      `切向跨度超过壁面目标的壁长占比 ${exceed == null ? '未设置目标' : (100 * exceed).toFixed(2) + '%'}。` +
      '切向跨度与法向高度分开测量；完整尺寸报告随结果包保存。';
    const layers = resolution.boundary_layers;
    if (layers && layers.status !== 'not_requested') {
      const percent = value => value == null ? '未测得' : (100 * value).toFixed(2) + '%';
      el.textContent += layers.status === 'no_layers_retained'
        ? ' 本次为纯网格回退，未保留边界层。'
        : ` 边界层最终保留 ${layers.retained_cells}/${layers.requested_cells} 个请求单元；` +
          `首层覆盖壁长 ${percent(layers.first_layer_wall_length_fraction)}，` +
          `完整请求层数覆盖壁长 ${percent(layers.full_requested_layers_wall_length_fraction)}。` +
          `首层法向高度超过目标的壁长占比 ${percent(layers.first_layer_height_exceedance_wall_length_fraction)}。` +
          (layers.status === 'incomplete' ? ' 层身份或连通性核对不完整，详见报告。' : '');
    }
  }

  function flowConvectionLabel(summary, catalog) {
    const scheme = catalog?.flowConvectionSchemes?.[summary.convection];
    const label = scheme?.label || (summary.convection === 'limited-linear'
      ? '线性迎风（限制重构）' : '一阶迎风');
    return `${label}${summary.convectionInferred ? '（旧摘要缺字段，按一阶迎风推断）' : ''}`;
  }

  function flowPressurePreconditionerLabel(summary, catalog) {
    if (summary.pressurePreconditionerInferred || summary.pressurePreconditioner === 'legacy-unspecified') {
      return '旧结果未记录';
    }
    const scheme = catalog?.flowPressurePreconditioners?.[summary.pressurePreconditioner];
    return scheme?.label || (summary.pressurePreconditioner === 'aggregation'
      ? '多重网格（试验）' : '标准（IC0）');
  }

  function flowOutletBackflowLabel(summary) {
    return summary.outletBackflow === 'normal-inlet' ? '允许法向回流（试验）' : '检测到回流时停止';
  }

  function renderFlowResult(summary, catalog) {
    const container = $('flowResult');
    container.replaceChildren();
    const convectionText = flowConvectionLabel(summary, catalog);
    const pressurePreconditionerText = flowPressurePreconditionerLabel(summary, catalog);
    const pressureText = summary.pressureDiscretization === 'shared-face-gauss'
      ? '共享面压力' : '旧结果未记录';
    const stateLine = document.createElement('div');
    stateLine.className = 'flow-state';
    stateLine.textContent = summary.converged
      ? `已收敛 · ${summary.iterations} 次迭代 · 对流：${convectionText}`
      : `到达 ${summary.iterations} 次迭代上限，结果有效但未收敛 · 对流：${convectionText}`;
    if (summary.temporalDiscretization) {
      stateLine.textContent = `已接受 t=${summary.acceptedTime.toPrecision(6)} s · 本次 ${summary.completedSteps} 步 · 最后一步内迭代 ${summary.iterations} 次`;
    }
    container.appendChild(stateLine);
    const rows = [
      ...(summary.temporalDiscretization ? [['最后一步步长（s）', summary.dt], ['最后一步最大 CFL', summary.maxCourant]] : []),
      ...(summary.timeStepControl === 'adaptive-cfl-retry' ? [['时间步控制', '自动 CFL / 重试'], ['拒绝的试算次数', String(summary.rejectedSteps)], ['CFL 上限', summary.targetCourant]] : []),
      ['对流格式', convectionText],
      ['压力求解', pressurePreconditionerText],
      ['压力—速度耦合', (summary.coupling ?? 'simple').toUpperCase()],
      ['线性迭代精度', summary.linearPolicy === 'adaptive' ? (summary.converged ? '自适应 · 严格复核通过' : '自适应 · 尚未收敛') : '固定精度'],
      ['速度松弛系数', summary.velocityRelaxation ?? 0.6],
      ['压力松弛系数', summary.pressureRelaxation ?? 0.25],
      ...(summary.coupling === 'pimple' ? [['动量校正次数', summary.momentumCorrectors]] : []),
      ['压力校正次数', summary.pressureCorrectionPasses ?? 4],
      ['稳态加速', summary.steadyAcceleration === 'anderson'
        ? `历史迭代 · 接受 ${summary.accelerationAccepted} / 舍弃 ${summary.accelerationRejected}` : '关闭'],
      ['出口回流', flowOutletBackflowLabel(summary)],
      ['压力离散', pressureText],
      ['局部连续性（无量纲）', summary.continuity],
      ['全局不平衡（m²/s）', summary.globalImbalance],
      ['全局相对不平衡', summary.globalRelativeImbalance],
      ['速度变化', summary.velocityChange],
      ['流动停止容差', summary.tolerance],
      ['压力变化', summary.pressureChange],
      ['动量残差', summary.momentumResidual]
    ];
    if (summary.outletBackflowFaces !== undefined) rows.push(['回流出口面数', String(summary.outletBackflowFaces)]);
    if (summary.outletInflow !== undefined) rows.push(['出口流入量（m²/s）', summary.outletInflow]);
    if (summary.namedBoundaryFluxes) {
      rows.push(['边界体积流量', '每单位厚度；净流出为正，净流入为负']);
      for (const flux of summary.namedBoundaryFluxes) {
        rows.push([`${flux.name} · 净流量（m²/s）`, flux.net],
          [`${flux.name} · 流入／流出（m²/s）`, `${flux.inflow.toExponential(3)} / ${flux.outflow.toExponential(3)}`],
          [`${flux.name} · 平均法向速度（m/s）`, flux.normalMeanVelocity]);
      }
    }
    if (summary.namedWallLoads) {
      rows.push(['壁面载荷', '流体作用于壁面；力÷ρ÷单位厚度'],
        ['力矩基准', '原点 (0,0)，逆时针为正']);
      for (const load of summary.namedWallLoads) {
        rows.push([`${load.name} · Fx（m³/s²）`, load.forceX],
          [`${load.name} · Fy（m³/s²）`, load.forceY],
          [`${load.name} · 力矩（m⁴/s²）`, load.torque]);
      }
    }
    for (const [label, value] of rows) {
      const item = document.createElement('div');
      const caption = document.createElement('span');
      caption.textContent = label;
      const number = document.createElement('b');
      number.textContent = typeof value === 'number' ? value.toExponential(3) : value;
      item.append(caption, number);
      container.appendChild(item);
    }
    container.hidden = false;
  }

  function renderThermalResult(summary, last) {
    const container = $('thermalResult');
    container.replaceChildren();
    const heading = document.createElement('div');
    heading.className = 'flow-state';
    heading.textContent = `温度与流动已接受 t=${Number(summary.acceptedTime ?? summary.time).toPrecision(6)} s · 本次 ${summary.steps} 步 · 单向恒物性输运`;
    container.appendChild(heading);
    for (const [label, value] of [
      ['时间步（s）', summary.dt], ['最低温度（K）', summary.minValue], ['最高温度（K）', summary.maxValue],
      ['热扩散率（m²/s）', summary.diffusivity], ['温度积分（K·m²）', last?.heatContent],
      ['全域平衡（K·m²/s）', last?.globalBalance], ['温度残差', last?.scalarResidual], ['最后一步最大 CFL', last?.maxCourant]
    ]) {
      const item = document.createElement('div');
      const caption = document.createElement('span');
      const number = document.createElement('b');
      caption.textContent = label;
      number.textContent = Number.isFinite(value) ? value.toExponential(4) : '未记录';
      item.append(caption, number);
      container.appendChild(item);
    }
    container.hidden = false;
  }

  function renderTimeSeriesMonitor(rows, { timelineId, svgId, metricId, captionId, prefix }) {
    const items = rows || [];
    $(timelineId).hidden = items.length === 0;
    const svg = $(svgId);
    svg.replaceChildren();
    if (!items.length) return;
    const metric = $(metricId).value;
    const data = items.filter(row => Number.isFinite(row[metric]));
    if (!data.length) return;
    const stride = Math.max(1, Math.ceil(data.length / 700));
    const sampled = data.filter((_r, i) => i % stride === 0);
    if (sampled.at(-1) !== data.at(-1)) sampled.push(data.at(-1));
    const t0 = data[0].time, t1 = data.at(-1).time;
    let lo = Infinity, hi = -Infinity;
    for (const r of data) { lo = Math.min(lo, r[metric]); hi = Math.max(hi, r[metric]); }
    const margin = hi === lo ? Math.max(1e-12, Math.abs(hi) * 0.02) : 0.05 * (hi - lo);
    lo -= margin; hi += margin;
    const x = t => 94 + 592 * (t - t0) / (t1 - t0 || 1);
    const y = v => 100 - 82 * (v - lo) / (hi - lo);
    const add = (tag, attrs, text) => {
      const el = document.createElementNS('http://www.w3.org/2000/svg', tag);
      for (const [key, value] of Object.entries(attrs)) el.setAttribute(key, String(value));
      if (text !== undefined) el.textContent = text;
      svg.appendChild(el);
    };
    add('path', { d: 'M94 12 V100 H690', fill: 'none', stroke: 'currentColor', opacity: 0.3 });
    add('path', { d: sampled.map((r, i) => `${i ? 'L' : 'M'}${x(r.time)},${y(r[metric])}`).join(' '), fill: 'none', stroke: 'currentColor', 'stroke-width': 1.8 });
    const last = sampled.at(-1);
    add('circle', { cx: x(last.time), cy: y(last[metric]), r: 2.5, fill: 'currentColor' });
    const digits = Math.min(12, Math.max(3, 2 + Math.ceil(Math.log10(Math.max(Math.abs(lo), Math.abs(hi)) / (hi - lo) || 1))));
    for (const [tx, ty, label] of [
      [90, 20, hi.toPrecision(digits)], [90, 101, lo.toPrecision(digits)],
      [135, 122, `${t0.toPrecision(4)} s`], [675, 122, `${t1.toPrecision(4)} s`]
    ]) {
      add('text', { x: tx, y: ty, fill: 'currentColor', 'text-anchor': 'end', 'font-size': 11 }, label);
    }
    $(captionId).textContent = `${prefix} · 最新 ${last[metric].toPrecision(Math.max(5, digits))} · 曲线按需抽样显示，导出 CSV 保留本次全部步数。`;
  }

  window.CartMeshResultsView = Object.freeze({
    renderCounters,
    renderGates,
    renderHistogram,
    renderLegend,
    renderResolution,
    renderFlowResult,
    renderThermalResult,
    renderTimeSeriesMonitor
  });
})();
