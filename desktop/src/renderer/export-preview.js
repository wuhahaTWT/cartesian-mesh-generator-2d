'use strict';

// Draw the exported solver mesh independently of the interactive camera/theme.
window.CartMeshExport = {
  render(mesh, result) {
    if (!mesh?.cells?.length) throw new Error('没有可绘制的最终网格。');
    const canvas = document.createElement('canvas');
    canvas.width = 2000; canvas.height = 1080;
    const ctx = canvas.getContext('2d');
    ctx.fillStyle = '#fff'; ctx.fillRect(0, 0, canvas.width, canvas.height);
    const text = (value, x, y, size = 24, align = 'left', color = '#252d36') => {
      ctx.fillStyle = color; ctx.font = `${size}px -apple-system, sans-serif`;
      ctx.textAlign = align; ctx.fillText(value, x, y);
    };
    const n = mesh.cells.length.toLocaleString('en-US');
    text(`网格预览  ·  ${n} 个单元`, 1000, 52, 32, 'center');
    text('全景 · 最终求解网格', 510, 98, 24, 'center');
    text('壁面局部 · 真实单元', 1500, 98, 24, 'center');

    const b = mesh.bounds;
    const span = Math.max(b.maxX - b.minX, b.maxY - b.minY);
    const wall = mesh.edges.filter(edge => edge.patch === 1);
    const edge = wall[Math.floor(wall.length / 3)] || mesh.edges[0];
    const a = mesh.vertices[edge.a], z = mesh.vertices[edge.b];
    const cx = (a[0] + z[0]) / 2, cy = (a[1] + z[1]) / 2;
    const zoomSpan = Math.max(span * 0.075, Math.hypot(a[0] - z[0], a[1] - z[1]) * 12);
    const detail = { minX: cx - zoomSpan / 2, maxX: cx + zoomSpan / 2,
      minY: cy - zoomSpan / 2, maxY: cy + zoomSpan / 2 };
    const layers = new Map();
    for (const column of result.resolution?.boundary_layers?.columns || []) {
      for (const [index, id] of (column.solver_cell_ids || []).entries()) {
        if (Number.isInteger(id) && id >= 0 && id < mesh.cells.length) layers.set(id, index);
      }
    }
    const colors = ['#287ba7', '#58b5b9', '#acd5c3', '#f2ebbc'];
    function panel(frame, rect, markDetail) {
      const pad = Math.max(frame.maxX - frame.minX, frame.maxY - frame.minY) * 0.025;
      const scale = Math.min(rect.w / (frame.maxX - frame.minX + pad * 2), rect.h / (frame.maxY - frame.minY + pad * 2));
      const centerX = (frame.minX + frame.maxX) / 2, centerY = (frame.minY + frame.maxY) / 2;
      const X = x => rect.x + rect.w / 2 + (x - centerX) * scale;
      const Y = y => rect.y + rect.h / 2 - (y - centerY) * scale;
      const visible = (p, q) => !(Math.max(p[0], q[0]) < centerX - rect.w / scale / 2 ||
        Math.min(p[0], q[0]) > centerX + rect.w / scale / 2 ||
        Math.max(p[1], q[1]) < centerY - rect.h / scale / 2 ||
        Math.min(p[1], q[1]) > centerY + rect.h / scale / 2);
      ctx.save(); ctx.beginPath(); ctx.rect(rect.x, rect.y, rect.w, rect.h); ctx.clip();
      // Layer membership comes from the final mesh's lineage report.
      for (const [id, layer] of layers) {
        const ids = mesh.cells[id].vertices;
        ctx.beginPath();
        ids.forEach((v, i) => { const p = mesh.vertices[v]; i ? ctx.lineTo(X(p[0]), Y(p[1])) : ctx.moveTo(X(p[0]), Y(p[1])); });
        ctx.closePath(); ctx.fillStyle = colors[Math.min(layer, colors.length - 1)]; ctx.fill();
      }
      for (const boundary of [false, true]) {
        ctx.beginPath();
        for (const e of mesh.edges) {
          if ((e.patch !== 0) !== boundary) continue;
          const p = mesh.vertices[e.a], q = mesh.vertices[e.b];
          if (!visible(p, q)) continue;
          ctx.moveTo(X(p[0]), Y(p[1])); ctx.lineTo(X(q[0]), Y(q[1]));
        }
        ctx.strokeStyle = boundary ? '#344c60' : '#94a2b0';
        ctx.lineWidth = boundary ? 1.25 : 0.42; ctx.stroke();
      }
      if (markDetail) {
        ctx.strokeStyle = '#277eac'; ctx.lineWidth = 2; ctx.setLineDash([7, 5]);
        ctx.strokeRect(X(detail.minX), Y(detail.maxY), zoomSpan * scale, zoomSpan * scale);
      }
      ctx.restore();
      ctx.strokeStyle = '#55616e'; ctx.lineWidth = 1; ctx.strokeRect(rect.x, rect.y, rect.w, rect.h);
      const format = value => Number(value.toPrecision(4)).toString();
      for (let i = 0; i <= 4; i++) {
        const x = centerX + (i / 4 - .5) * rect.w / scale;
        const y = centerY + (.5 - i / 4) * rect.h / scale;
        text(format(x), rect.x + rect.w * i / 4, rect.y + rect.h + 29, 18, 'center');
        text(format(y), rect.x - 12, rect.y + rect.h * i / 4 + 6, 18, 'right');
      }
      text('x (m)', rect.x + rect.w / 2, rect.y + rect.h + 58, 20, 'center');
      text('y (m)', rect.x, rect.y - 14, 20);
    }
    panel(b, { x: 105, y: 145, w: 840, h: 790 }, true);
    panel(detail, { x: 1100, y: 145, w: 840, h: 790 }, false);
    const gate = value => (value?.pass ?? value?.valid) === true ? '通过' : (value?.pass ?? value?.valid) === false ? '未通过' : '未检查';
    text(`内部拓扑：${gate(result.gates?.topology)}   |   内部 Solver：${gate(result.gates?.solver)}   |   外部 CFD 检查：需另行运行`, 1000, 1042, 22, 'center');
    return canvas.toDataURL('image/png');
  }
};
