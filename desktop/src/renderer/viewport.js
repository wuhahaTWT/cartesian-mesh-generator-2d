'use strict';

// Mesh viewport: fit, wheel zoom at the cursor, drag to pan, cells coloured by
// refinement level.
//
// Zoom is not a nicety here.  A 10-body-span far field makes the body one
// twenty-first of the domain width, so a fit-to-domain view shows the wall as a few
// pixels and the whole point of the size field becomes invisible.

// These files are classic scripts sharing one global lexical scope, so the body is
// wrapped and only window.MeshView escapes.  Without it a top-level `levelColour`
// here collides with app.js destructuring the same name back off the namespace.
(function () {

// Sequential ramp, coarse -> fine.  One ramp is used by both the mesh and the level
// histogram so the two read as the same scale.
const RAMP = ['#1b3a4b', '#1d4f5e', '#216b66', '#3d8560', '#7d9a4e', '#b8a344', '#dd9b3c', '#f07f3c'];
function levelColour(level, minLevel, maxLevel) {
  if (!(maxLevel > minLevel)) return RAMP[RAMP.length - 1];
  const t = (level - minLevel) / (maxLevel - minLevel);
  const scaled = t * (RAMP.length - 1);
  const index = Math.min(RAMP.length - 1, Math.max(0, Math.round(scaled)));
  return RAMP[index];
}

// Three display styles.  The level colouring reads the size field at a glance, but a
// mesh figure for a paper wants plain lines, so `fill: null` means draw edges only.
const THEMES = {
  level: { fill: 'level', cellEdge: 'rgba(12,20,26,0.55)', domain: 'rgba(150,178,196,0.55)',
           wall: '#ff5a1f', unclassified: '#ff3b6b', region: '#7fd7ff', vertex: '#ffd8c4' },
  light: { fill: null, cellEdge: 'rgba(40,54,64,0.55)', domain: '#5a6b76',
           wall: '#c62828', unclassified: '#ff3b6b', region: '#1d78a8', vertex: '#c62828' },
  dark: { fill: null, cellEdge: 'rgba(198,214,226,0.42)', domain: 'rgba(198,214,226,0.75)',
          wall: '#ff5a1f', unclassified: '#ff3b6b', region: '#7fd7ff', vertex: '#ffd8c4' }
};
const CELLS_PER_PATH = 64;

class Viewport {
  constructor(canvas) {
    this.canvas = canvas;
    this.mesh = null;
    this.outline = null;
    this.regions = [];
    this.frame = null;
    this.mode = 'level';
    this.showGrid = true;
    this.showRegions = true;
    this.scale = 1;
    this.offset = { x: 0, y: 0 };
    this.dragging = null;
    this.meshCache = null;
    this.pendingDraw = null;
    this.attachInput();
    this.resizeObserver = new ResizeObserver(() => this.requestDraw());
    this.resizeObserver.observe(canvas);
  }

  theme() { return THEMES[this.mode] || THEMES.level; }

  attachInput() {
    this.canvas.addEventListener('wheel', event => {
      if (!this.mesh && !this.outline) return;
      event.preventDefault();
      const rect = this.canvas.getBoundingClientRect();
      const px = event.clientX - rect.left;
      const py = event.clientY - rect.top;
      // Keep the world point under the cursor fixed: the only zoom that feels like
      // direct manipulation.
      const before = this.toWorld(px, py);
      const factor = Math.exp(-event.deltaY * 0.0016);
      this.scale = Math.min(1e9, Math.max(1e-9, this.scale * factor));
      const after = this.toWorld(px, py);
      this.offset.x += before.x - after.x;
      this.offset.y += before.y - after.y;
      this.requestDraw();
    }, { passive: false });

    this.canvas.addEventListener('pointerdown', event => {
      this.dragging = { x: event.clientX, y: event.clientY };
      this.canvas.setPointerCapture(event.pointerId);
      this.canvas.classList.add('grabbing');
    });
    this.canvas.addEventListener('pointermove', event => {
      if (!this.dragging) return;
      // scale is CSS pixels per world unit, so a screen delta divides straight by it.
      // y is negated once because the world axis points up.
      this.offset.x -= (event.clientX - this.dragging.x) / this.scale;
      this.offset.y += (event.clientY - this.dragging.y) / this.scale;
      this.dragging = { x: event.clientX, y: event.clientY };
      this.requestDraw();
    });
    const release = () => {
      this.dragging = null;
      this.canvas.classList.remove('grabbing');
    };
    this.canvas.addEventListener('pointerup', release);
    this.canvas.addEventListener('pointercancel', release);
  }

  size() {
    const rect = this.canvas.getBoundingClientRect();
    return { width: Math.max(1, rect.width), height: Math.max(1, rect.height) };
  }

  toWorld(px, py) {
    const { height } = this.size();
    return { x: px / this.scale + this.offset.x, y: (height - py) / this.scale + this.offset.y };
  }

  fitTo(bounds, margin = 0.06) {
    if (!bounds) return;
    const { width, height } = this.size();
    this.lastSize = { width, height };
    const spanX = Math.max(1e-300, bounds.maxX - bounds.minX);
    const spanY = Math.max(1e-300, bounds.maxY - bounds.minY);
    this.scale = (1 - 2 * margin) * Math.min(width / spanX, height / spanY);
    this.offset.x = (bounds.minX + bounds.maxX) / 2 - width / (2 * this.scale);
    this.offset.y = (bounds.minY + bounds.maxY) / 2 - height / (2 * this.scale);
    this.draw();
  }

  setMesh(mesh) {
    this.mesh = mesh;
    this.outline = null;
    this.meshCache = null;
    this.fitTo(mesh.bounds);
  }

  // Before a mesh exists the chosen geometry is still worth drawing: it is how the
  // user confirms the importer read the file they meant.
  setOutline(loops) {
    this.mesh = null;
    this.meshCache = null;
    this.outline = loops;
    const points = loops.flat();
    if (!points.length) return;
    const xs = points.map(p => p[0]);
    const ys = points.map(p => p[1]);
    this.fitTo({ minX: Math.min(...xs), minY: Math.min(...ys),
                 maxX: Math.max(...xs), maxY: Math.max(...ys) }, 0.12);
  }

  clear() {
    this.mesh = null;
    this.outline = null;
    this.meshCache = null;
    this.draw();
  }

  // Wheel, pointer and resize events can arrive several times before Chromium paints
  // one frame.  Keep only the newest view transform instead of redrawing the same
  // large mesh for every event.  draw() itself remains synchronous for exports.
  requestDraw() {
    if (this.pendingDraw !== null) return;
    this.pendingDraw = window.requestAnimationFrame(() => {
      this.pendingDraw = null;
      this.draw();
    });
  }

  cancelPendingDraw() {
    if (this.pendingDraw === null) return;
    window.cancelAnimationFrame(this.pendingDraw);
    this.pendingDraw = null;
  }

  context() {
    const dpr = window.devicePixelRatio || 1;
    const { width, height } = this.size();
    if (this.lastSize && (width !== this.lastSize.width || height !== this.lastSize.height)) {
      const centreX = this.offset.x + this.lastSize.width / (2 * this.scale);
      const centreY = this.offset.y + this.lastSize.height / (2 * this.scale);
      // Use one reversible viewport measure. min(widthRatio,heightRatio)
      // shrinks on a sidebar round trip and never restores the original zoom.
      this.scale *= Math.min(width, height) / Math.min(this.lastSize.width, this.lastSize.height);
      this.offset.x = centreX - width / (2 * this.scale);
      this.offset.y = centreY - height / (2 * this.scale);
    }
    this.lastSize = { width, height };
    const pixelWidth = Math.floor(width * dpr);
    const pixelHeight = Math.floor(height * dpr);
    // Assigning either dimension clears the canvas and reallocates its backing store.
    // Do that on a real resize only, not for every pan or zoom frame.
    if (this.canvas.width !== pixelWidth) this.canvas.width = pixelWidth;
    if (this.canvas.height !== pixelHeight) this.canvas.height = pixelHeight;
    const ctx = this.canvas.getContext('2d');
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    ctx.clearRect(0, 0, width, height);
    return { ctx, width, height };
  }

  draw() {
    this.cancelPendingDraw();
    const { ctx, width, height } = this.context();
    const theme = this.theme();
    const project = ([x, y]) => [
      (x - this.offset.x) * this.scale,
      height - (y - this.offset.y) * this.scale
    ];
    if (this.outline) {
      this.drawOutline(ctx, project, theme);
      this.drawRegions(ctx, project, theme);
      return;
    }
    if (!this.mesh) return;
    const mesh = this.mesh;
    const cache = this.cachedMeshPaths(mesh);
    const overscan = 3 / this.scale;
    const visible = {
      minX: this.offset.x - overscan,
      maxX: this.offset.x + width / this.scale + overscan,
      minY: this.offset.y - overscan,
      maxY: this.offset.y + height / this.scale + overscan
    };
    ctx.save();
    // Paths stay in world coordinates.  Pan and zoom now change one canvas transform
    // rather than projecting every vertex in JavaScript on every frame.
    ctx.transform(this.scale, 0, 0, -this.scale,
                  -this.offset.x * this.scale,
                  height + this.offset.y * this.scale);
    for (const { level, chunks } of cache.levels) {
      const visibleChunks = chunks.filter(chunk => this.intersects(chunk.bounds, visible));
      if (theme.fill === 'level') {
        ctx.fillStyle = levelColour(level, mesh.minLevel, mesh.maxLevel);
        for (const chunk of visibleChunks) ctx.fill(chunk.path);
      }
      // Outlining cells narrower than ~3 px turns the mesh into a solid block and
      // costs the most time on the largest meshes, so it is skipped there.  Without a
      // fill there would be nothing left to see, so the floor drops to 1 px.
      const onScreen = (cache.domainSpan / Math.pow(2, level)) * this.scale;
      if (this.showGrid && onScreen >= (theme.fill ? 3 : 1)) {
        ctx.strokeStyle = theme.cellEdge;
        ctx.lineWidth = Math.min(theme.fill ? 1 : 0.7,
                                 Math.max(0.35, onScreen / 12)) / this.scale;
        for (const chunk of visibleChunks) ctx.stroke(chunk.path);
      }
    }
    this.drawBoundaries(ctx, cache, theme);
    ctx.restore();
    this.drawRegions(ctx, project, theme);
  }

  cachedMeshPaths(mesh) {
    if (this.meshCache && this.meshCache.mesh === mesh) return this.meshCache;
    const levels = new Map();
    for (const cell of mesh.cells) {
      if (!levels.has(cell.level)) levels.set(cell.level, []);
      const chunks = levels.get(cell.level);
      if (!chunks.length || chunks[chunks.length - 1].cells === CELLS_PER_PATH) {
        chunks.push({ path: new Path2D(), cells: 0,
          bounds: { minX: Infinity, minY: Infinity, maxX: -Infinity, maxY: -Infinity } });
      }
      const chunk = chunks[chunks.length - 1];
      const path = chunk.path;
      const first = mesh.vertices[cell.vertices[0]];
      path.moveTo(first[0], first[1]);
      this.extendBounds(chunk.bounds, first);
      for (let i = 1; i < cell.vertices.length; i++) {
        const point = mesh.vertices[cell.vertices[i]];
        path.lineTo(point[0], point[1]);
        this.extendBounds(chunk.bounds, point);
      }
      path.closePath();
      chunk.cells++;
    }
    const boundaries = new Map([[1, new Path2D()], [2, new Path2D()], [3, new Path2D()]]);
    for (const edge of mesh.edges) {
      const path = boundaries.get(edge.patch);
      if (!path) continue;
      const a = mesh.vertices[edge.a];
      const b = mesh.vertices[edge.b];
      path.moveTo(a[0], a[1]);
      path.lineTo(b[0], b[1]);
    }
    this.meshCache = {
      mesh,
      levels: [...levels].sort((a, b) => a[0] - b[0])
        .map(([level, chunks]) => ({ level, chunks })),
      boundaries,
      domainSpan: Math.max(mesh.bounds.maxX - mesh.bounds.minX,
                           mesh.bounds.maxY - mesh.bounds.minY)
    };
    return this.meshCache;
  }

  extendBounds(bounds, [x, y]) {
    if (x < bounds.minX) bounds.minX = x;
    if (x > bounds.maxX) bounds.maxX = x;
    if (y < bounds.minY) bounds.minY = y;
    if (y > bounds.maxY) bounds.maxY = y;
  }

  intersects(a, b) {
    return a.maxX >= b.minX && a.minX <= b.maxX &&
           a.maxY >= b.minY && a.minY <= b.maxY;
  }

  // Hand-placed regions are stated in body spans about the body centre, which is the
  // frame the panel edits in; drawing them in the same frame is what makes the numbers
  // checkable by eye.
  drawRegions(ctx, project, theme) {
    if (!this.showRegions || !this.frame || !this.regions.length) return;
    // refine() clips a region to the domain, so the overlay clips too: a rectangle
    // drawn outside the computational domain would promise refinement that never
    // happens.  The editor still shows the numbers as typed.
    const limit = this.mesh ? this.mesh.bounds : null;
    ctx.save();
    ctx.setLineDash([6, 4]);
    ctx.strokeStyle = theme.region;
    ctx.lineWidth = 1.4;
    ctx.font = '11px ui-monospace, Menlo, monospace';
    this.regions.forEach((box, index) => {
      let x0 = this.frame.centreX + box.xmin * this.frame.bodySpan;
      let x1 = this.frame.centreX + box.xmax * this.frame.bodySpan;
      let y0 = this.frame.centreY + box.ymin * this.frame.bodySpan;
      let y1 = this.frame.centreY + box.ymax * this.frame.bodySpan;
      if (limit) {
        x0 = Math.max(x0, limit.minX); x1 = Math.min(x1, limit.maxX);
        y0 = Math.max(y0, limit.minY); y1 = Math.min(y1, limit.maxY);
        if (!(x1 > x0) || !(y1 > y0)) return;
      }
      const a = project([x0, y0]);
      const b = project([x1, y1]);
      ctx.globalAlpha = box.active ? 1 : 0.6;
      ctx.strokeRect(Math.min(a[0], b[0]), Math.min(a[1], b[1]),
                     Math.abs(b[0] - a[0]), Math.abs(b[1] - a[1]));
      ctx.fillStyle = theme.region;
      ctx.fillText(`R${index + 1}  -${box.levelsBelowWall}`,
                   Math.min(a[0], b[0]) + 4, Math.min(a[1], b[1]) + 13);
    });
    ctx.restore();
  }

  drawBoundaries(ctx, cache, theme) {
    const stroke = (patch, colour, lineWidth) => {
      ctx.strokeStyle = colour;
      ctx.lineWidth = lineWidth / this.scale;
      ctx.stroke(cache.boundaries.get(patch));
    };
    stroke(2, theme.domain, 1);
    stroke(3, theme.unclassified, 2.4);
    stroke(1, theme.wall, 1.8);
  }

  drawOutline(ctx, project, theme) {
    ctx.lineJoin = 'round';
    for (const loop of this.outline) {
      ctx.beginPath();
      loop.forEach((point, index) => {
        const [x, y] = project(point);
        if (index === 0) ctx.moveTo(x, y); else ctx.lineTo(x, y);
      });
      ctx.closePath();
      ctx.fillStyle = this.mode === 'light' ? 'rgba(198,40,40,0.08)' : 'rgba(255,90,31,0.10)';
      ctx.fill();
      ctx.strokeStyle = theme.wall;
      ctx.lineWidth = 1.6;
      ctx.stroke();
      // Vertices are the wall's tangential resolution, and it being too coarse is the
      // known cause of the hybrid path's ceiling, so they are worth showing.
      ctx.fillStyle = theme.vertex;
      for (const point of loop) {
        const [x, y] = project(point);
        ctx.fillRect(x - 1.4, y - 1.4, 2.8, 2.8);
      }
    }
  }
}

// The renderer runs with contextIsolation on and cannot require(), so the one export
// is a namespace on window.
window.MeshView = { Viewport, levelColour, RAMP };

})();
