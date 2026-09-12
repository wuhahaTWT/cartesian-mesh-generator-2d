(function (root, factory) {
  'use strict';
  const api = factory();
  if (typeof module === 'object' && module.exports) module.exports = api;
  else root.CartMeshRaster = api;
}(typeof self !== 'undefined' ? self : (typeof window !== 'undefined' ? window : globalThis), function () {
  'use strict';

  const MAX_PIXELS = 2000000;
  const MAX_DIMENSION = 4096;
  const MAX_LOOPS = 64;
  const MAX_POINTS = 4096;
  const MAX_EXACT_POINTS = 65536;
  const MAX_RAW_EDGES = 262144;
  const MODES = new Set(['auto', 'dark', 'light', 'background', 'alpha']);

  function fail(message) { throw new Error(message); }

  function validateInput(image, options) {
    if (!image || typeof image !== 'object') fail('图片像素数据缺失。');
    const width = image.width;
    const height = image.height;
    if (!Number.isInteger(width) || !Number.isInteger(height) || width < 2 || height < 2)
      fail('图片宽高必须是至少 2 像素的整数。');
    if (width > MAX_DIMENSION || height > MAX_DIMENSION || width * height > MAX_PIXELS)
      fail('图片过大；请裁剪或缩小到 200 万像素以内（单边不超过 4096）。');
    if (!image.data || !ArrayBuffer.isView(image.data) || image.data.BYTES_PER_ELEMENT !== 1 ||
        image.data.length !== width * height * 4)
      fail('RGBA 像素字节长度与图片宽高不一致。');

    const mode = options.mode == null ? 'auto' : String(options.mode);
    if (!MODES.has(mode)) fail(`未知图片分割模式：${mode}`);
    const selection = options.selection == null ? 'largest' : String(options.selection);
    if (selection !== 'largest' && selection !== 'all') fail(`未知主体选择方式：${selection}`);
    let threshold = null;
    if (options.threshold !== undefined && options.threshold !== null && options.threshold !== '') {
      threshold = Number(options.threshold);
      if (!Number.isFinite(threshold) || threshold < 0 || threshold > 255)
        fail('阈值必须在 0 到 255 之间。');
      threshold = Math.round(threshold);
    }
    const defaultMinArea = Math.max(4, Math.ceil(width * height * 0.0002));
    const minAreaPx = options.minAreaPx == null ? defaultMinArea : Number(options.minAreaPx);
    if (!Number.isInteger(minAreaPx) || minAreaPx < 1 || minAreaPx > width * height)
      fail('最小主体面积必须是有效的正整数像素数。');
    const simplifyTolerance = options.simplifyTolerance == null ? 1 : Number(options.simplifyTolerance);
    if (!Number.isFinite(simplifyTolerance) || simplifyTolerance < 0 || simplifyTolerance > 8)
      fail('轮廓简化容差必须在 0 到 8 像素之间。');
    let seed = null;
    if (options.seed != null) {
      const x = Number(options.seed.x), y = Number(options.seed.y);
      if (!Number.isFinite(x) || !Number.isFinite(y)) fail('主体选点坐标无效。');
      seed = { x: Math.floor(x), y: Math.floor(y) };
      if (seed.x < 0 || seed.y < 0 || seed.x >= width || seed.y >= height)
        fail('主体选点落在图片范围之外。');
    }
    return { width, height, data: image.data, mode, selection, threshold, minAreaPx,
      simplifyTolerance, seed, fillHoles: options.fillHoles === true };
  }

  function borderIndices(width, height) {
    const out = [];
    for (let x = 0; x < width; x++) out.push(x, (height - 1) * width + x);
    for (let y = 1; y + 1 < height; y++) out.push(y * width, y * width + width - 1);
    return out;
  }

  function median(values) {
    values.sort((a, b) => a - b);
    const middle = values.length >> 1;
    return values.length % 2 ? values[middle] : (values[middle - 1] + values[middle]) / 2;
  }

  function otsu(values) {
    const histogram = new Uint32Array(256);
    for (let i = 0; i < values.length; i++) histogram[values[i]]++;
    let sum = 0;
    for (let i = 0; i < 256; i++) sum += i * histogram[i];
    let backgroundCount = 0, backgroundSum = 0, bestVariance = -1, best = 0;
    for (let i = 0; i < 255; i++) {
      backgroundCount += histogram[i];
      backgroundSum += i * histogram[i];
      if (!backgroundCount) continue;
      const foregroundCount = values.length - backgroundCount;
      if (!foregroundCount) break;
      const mean0 = backgroundSum / backgroundCount;
      const mean1 = (sum - backgroundSum) / foregroundCount;
      const variance = backgroundCount * foregroundCount * (mean0 - mean1) * (mean0 - mean1);
      if (variance > bestVariance) { bestVariance = variance; best = i; }
    }
    return best;
  }

  function makeMask(input, warnings) {
    const { width, height, data } = input;
    const total = width * height;
    const border = borderIndices(width, height);
    let mode = input.mode;
    let threshold = input.threshold;
    let values = null;
    let backgroundRgb = null;

    if (mode === 'auto') {
      let transparentBorder = 0, visiblePixels = 0;
      for (const index of border) if (data[index * 4 + 3] <= 16) transparentBorder++;
      for (let i = 0; i < total; i++) if (data[i * 4 + 3] > 32) visiblePixels++;
      if (transparentBorder / border.length >= 0.9 && visiblePixels > 0 && visiblePixels < total) {
        mode = 'alpha';
      } else {
        values = grayscale(data, total);
        const autoThreshold = threshold == null ? otsu(values) : threshold;
        const borderLuma = median(border.map(index => values[index]));
        mode = borderLuma > autoThreshold ? 'dark' : 'light';
      }
    }

    if (mode === 'alpha') {
      values = new Uint8Array(total);
      for (let i = 0; i < total; i++) values[i] = data[i * 4 + 3];
      if (threshold == null) threshold = otsu(values);
    } else if (mode === 'background') {
      const rs = [], gs = [], bs = [];
      for (const index of border) {
        const p = index * 4;
        rs.push(data[p]); gs.push(data[p + 1]); bs.push(data[p + 2]);
      }
      backgroundRgb = [median(rs), median(gs), median(bs)];
      values = new Uint8Array(total);
      for (let i = 0; i < total; i++) {
        const p = i * 4;
        const dr = data[p] - backgroundRgb[0];
        const dg = data[p + 1] - backgroundRgb[1];
        const db = data[p + 2] - backgroundRgb[2];
        values[i] = Math.round(Math.sqrt((dr * dr + dg * dg + db * db) / 3));
      }
      if (threshold == null) threshold = Math.max(12, otsu(values));
      warnings.push('背景色模式只按图片边框的中位颜色做色差分割，不是任意场景的语义分割。');
    } else {
      if (!values) values = grayscale(data, total);
      if (threshold == null) threshold = otsu(values);
    }

    const mask = new Uint8Array(total);
    let foregroundPixels = 0;
    for (let i = 0; i < total; i++) {
      const selected = mode === 'dark' ? values[i] <= threshold
        : mode === 'light' ? values[i] > threshold
          : values[i] > threshold;
      if (selected) { mask[i] = 1; foregroundPixels++; }
    }
    if (!foregroundPixels) fail('阈值结果为空；请调整模式、阈值或裁剪图片。');
    if (foregroundPixels === total) fail('阈值结果覆盖整张图片；请调整模式、阈值或裁剪图片。');
    return { mask, mode, threshold, foregroundPixels, backgroundRgb };
  }

  function grayscale(data, total) {
    const values = new Uint8Array(total);
    for (let i = 0; i < total; i++) {
      const p = i * 4;
      const alpha = data[p + 3] / 255;
      const visible = 0.2126 * data[p] + 0.7152 * data[p + 1] + 0.0722 * data[p + 2];
      values[i] = Math.round(visible * alpha + 255 * (1 - alpha));
    }
    return values;
  }

  function components4(mask, width, height) {
    const labels = new Int32Array(mask.length);
    const queue = new Int32Array(mask.length);
    const components = [];
    for (let start = 0; start < mask.length; start++) {
      if (!mask[start] || labels[start]) continue;
      const id = components.length + 1;
      let head = 0, tail = 0, size = 0, touchesBorder = false;
      labels[start] = id; queue[tail++] = start;
      while (head < tail) {
        const index = queue[head++];
        const x = index % width, y = (index - x) / width;
        size++;
        if (x === 0 || y === 0 || x + 1 === width || y + 1 === height) touchesBorder = true;
        if (x > 0) visit(index - 1);
        if (x + 1 < width) visit(index + 1);
        if (y > 0) visit(index - width);
        if (y + 1 < height) visit(index + width);
      }
      components.push({ id, size, touchesBorder, first: start });
      function visit(next) {
        if (mask[next] && !labels[next]) { labels[next] = id; queue[tail++] = next; }
      }
    }
    return { labels, components };
  }

  function selectComponents(segmentation, input, warnings) {
    const { labels, components } = components4(segmentation.mask, input.width, input.height);
    const kept = components.filter(component => component.size >= input.minAreaPx);
    const removed = components.filter(component => component.size < input.minAreaPx);
    const removedPixels = removed.reduce((sum, component) => sum + component.size, 0);
    if (removed.length)
      warnings.push(`已移除 ${removed.length} 个小噪点主体，共 ${removedPixels} 像素（面积门槛 ${input.minAreaPx} 像素）。`);
    if (!kept.length) fail(`没有主体达到最小面积 ${input.minAreaPx} 像素；请调低门槛或重新裁剪。`);

    let selected;
    if (input.seed) {
      const label = labels[input.seed.y * input.width + input.seed.x];
      if (!label) fail('主体选点落在背景上；请点选要保留的图形内部。');
      const component = components[label - 1];
      if (component.size < input.minAreaPx)
        fail(`主体选点对应的区域只有 ${component.size} 像素，低于最小面积 ${input.minAreaPx}。`);
      selected = [component];
    } else if (input.selection === 'all') {
      selected = kept;
    } else {
      selected = [kept.reduce((best, item) => item.size > best.size ? item : best)];
      if (kept.length > 1)
        warnings.push(`已按“最大主体”保留 1 个区域，另有 ${kept.length - 1} 个达到面积门槛的区域未选中。`);
    }
    if (selected.some(component => component.touchesBorder))
      fail('所选主体接触图片边缘，轮廓可能被裁断；请扩大画布边距或重新裁剪。');
    const selectedIds = new Set(selected.map(component => component.id));
    const mask = new Uint8Array(segmentation.mask.length);
    let selectedPixels = 0;
    for (let i = 0; i < labels.length; i++) {
      if (selectedIds.has(labels[i])) { mask[i] = 1; selectedPixels++; }
    }
    return { mask, labels, components, selected, removed, removedPixels, selectedPixels, keptCount: kept.length };
  }

  function fillEnclosedBackground(mask, width, height) {
    const outside = new Uint8Array(mask.length);
    const queue = new Int32Array(mask.length);
    let head = 0, tail = 0;
    const add = index => { if (!mask[index] && !outside[index]) { outside[index] = 1; queue[tail++] = index; } };
    for (let x = 0; x < width; x++) { add(x); add((height - 1) * width + x); }
    for (let y = 1; y + 1 < height; y++) { add(y * width); add(y * width + width - 1); }
    while (head < tail) {
      const index = queue[head++], x = index % width, y = (index - x) / width;
      if (x > 0) add(index - 1);
      if (x + 1 < width) add(index + 1);
      if (y > 0) add(index - width);
      if (y + 1 < height) add(index + width);
    }
    const seen = new Uint8Array(mask.length);
    let holes = 0, area = 0;
    for (let start = 0; start < mask.length; start++) {
      if (mask[start] || outside[start] || seen[start]) continue;
      holes++;
      head = 0; tail = 0; seen[start] = 1; queue[tail++] = start;
      while (head < tail) {
        const index = queue[head++], x = index % width, y = (index - x) / width;
        mask[index] = 1; area++;
        const visit = next => {
          if (!mask[next] && !outside[next] && !seen[next]) { seen[next] = 1; queue[tail++] = next; }
        };
        if (x > 0) visit(index - 1);
        if (x + 1 < width) visit(index + 1);
        if (y > 0) visit(index - width);
        if (y + 1 < height) visit(index + width);
      }
    }
    return { holes, area };
  }

  function tracePixelBoundaries(mask, width, height) {
    const from = [], to = [], direction = [];
    const stride = width + 1;
    const add = (ax, ay, bx, by, dir) => {
      if (from.length >= MAX_RAW_EDGES)
        fail('图片边界细节过多；请裁剪、去噪或提高阈值后重试。');
      from.push(ay * stride + ax); to.push(by * stride + bx); direction.push(dir);
    };
    for (let y = 0; y < height; y++) for (let x = 0; x < width; x++) {
      const index = y * width + x;
      if (!mask[index]) continue;
      if (y === 0 || !mask[index - width]) add(x, y, x + 1, y, 0);
      if (x + 1 === width || !mask[index + 1]) add(x + 1, y, x + 1, y + 1, 1);
      if (y + 1 === height || !mask[index + width]) add(x + 1, y + 1, x, y + 1, 2);
      if (x === 0 || !mask[index - 1]) add(x, y + 1, x, y, 3);
    }
    const outgoing = new Map();
    for (let i = 0; i < from.length; i++) {
      const list = outgoing.get(from[i]);
      if (list) list.push(i); else outgoing.set(from[i], [i]);
    }
    const used = new Uint8Array(from.length), loops = [];
    for (let first = 0; first < from.length; first++) {
      if (used[first]) continue;
      const loop = [];
      const startVertex = from[first];
      let edge = first;
      while (true) {
        if (used[edge]) fail('像素边界拓扑不闭合；请调整阈值或去噪。');
        used[edge] = 1;
        const vertex = from[edge], x = vertex % stride, y = (vertex - x) / stride;
        loop.push([x, y]);
        const end = to[edge];
        if (end === startVertex) break;
        const candidates = (outgoing.get(end) || []).filter(index => !used[index]);
        if (!candidates.length) fail('像素边界拓扑不闭合；请调整阈值或去噪。');
        const incoming = direction[edge];
        const priorities = [(incoming + 1) % 4, incoming, (incoming + 3) % 4, (incoming + 2) % 4];
        edge = candidates.slice().sort((a, b) => priorities.indexOf(direction[a]) - priorities.indexOf(direction[b]))[0];
      }
      loops.push(loop);
      if (loops.length > MAX_LOOPS)
        fail(`识别得到超过 ${MAX_LOOPS} 个闭合环；请裁剪图片、只选主体或加强去噪。`);
    }
    return loops;
  }

  function signedArea(loop) {
    let area = 0;
    for (let i = 0; i < loop.length; i++) {
      const a = loop[i], b = loop[(i + 1) % loop.length];
      area += a[0] * b[1] - b[0] * a[1];
    }
    return area / 2;
  }

  function removeCollinear(loop) {
    if (loop.length <= 3) return loop.slice();
    const kept = [];
    for (let i = 0; i < loop.length; i++) {
      const a = loop[(i + loop.length - 1) % loop.length];
      const b = loop[i], c = loop[(i + 1) % loop.length];
      if ((b[0] - a[0]) * (c[1] - b[1]) !== (b[1] - a[1]) * (c[0] - b[0])) kept.push(b);
    }
    return kept.length >= 3 ? kept : loop.slice();
  }

  function pointSegmentDistance(point, a, b) {
    const dx = b[0] - a[0], dy = b[1] - a[1];
    if (dx === 0 && dy === 0) return Math.hypot(point[0] - a[0], point[1] - a[1]);
    const t = Math.max(0, Math.min(1, ((point[0] - a[0]) * dx + (point[1] - a[1]) * dy) / (dx * dx + dy * dy)));
    return Math.hypot(point[0] - (a[0] + t * dx), point[1] - (a[1] + t * dy));
  }

  function rdpOpen(points, tolerance) {
    if (points.length <= 2) return [points[0], points[points.length - 1]];
    // Use an explicit stack. A noisy 2M-pixel input can have tens of thousands of
    // corners, and recursive RDP would make the JS call stack an accidental limit.
    const keep = new Uint8Array(points.length);
    keep[0] = 1; keep[points.length - 1] = 1;
    const stack = [[0, points.length - 1]];
    while (stack.length) {
      const [first, last] = stack.pop();
      let farthest = -1, split = -1;
      for (let i = first + 1; i < last; i++) {
        const distance = pointSegmentDistance(points[i], points[first], points[last]);
        if (distance > farthest) { farthest = distance; split = i; }
      }
      if (farthest > tolerance) {
        keep[split] = 1;
        stack.push([first, split], [split, last]);
      }
    }
    return points.filter((point, index) => keep[index]);
  }

  function simplifyClosed(loop, tolerance) {
    const exact = removeCollinear(loop);
    if (!(tolerance > 0) || exact.length <= 3) return exact;
    let split = 1, farthest = -1;
    for (let i = 1; i < exact.length; i++) {
      const dx = exact[i][0] - exact[0][0], dy = exact[i][1] - exact[0][1];
      const distance2 = dx * dx + dy * dy;
      if (distance2 > farthest) { farthest = distance2; split = i; }
    }
    const first = exact.slice(0, split + 1);
    const second = exact.slice(split).concat([exact[0]]);
    const simplified = rdpOpen(first, tolerance).slice(0, -1).concat(rdpOpen(second, tolerance).slice(0, -1));
    return removeCollinear(simplified);
  }

  function orientation(a, b, c) {
    const value = (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0]);
    return value === 0 ? 0 : value > 0 ? 1 : -1;
  }

  function onSegment(a, b, p) {
    return orientation(a, b, p) === 0 && p[0] >= Math.min(a[0], b[0]) && p[0] <= Math.max(a[0], b[0])
      && p[1] >= Math.min(a[1], b[1]) && p[1] <= Math.max(a[1], b[1]);
  }

  function segmentsIntersect(a, b, c, d) {
    const o1 = orientation(a, b, c), o2 = orientation(a, b, d);
    const o3 = orientation(c, d, a), o4 = orientation(c, d, b);
    return (o1 * o2 < 0 && o3 * o4 < 0) || (o1 === 0 && onSegment(a, b, c)) ||
      (o2 === 0 && onSegment(a, b, d)) || (o3 === 0 && onSegment(c, d, a)) ||
      (o4 === 0 && onSegment(c, d, b));
  }

  function samePoint(a, b) { return a[0] === b[0] && a[1] === b[1]; }

  function simpleLoop(loop) {
    if (loop.length < 3 || signedArea(loop) === 0) return false;
    for (let i = 0; i < loop.length; i++) for (let j = i + 1; j < loop.length; j++) {
      if (j === i + 1 || (i === 0 && j + 1 === loop.length)) continue;
      if (segmentsIntersect(loop[i], loop[(i + 1) % loop.length], loop[j], loop[(j + 1) % loop.length])) return false;
    }
    return true;
  }

  function loopsDoNotCross(loops) {
    for (let a = 0; a < loops.length; a++) for (let b = a + 1; b < loops.length; b++) {
      const left = loops[a], right = loops[b];
      for (let i = 0; i < left.length; i++) for (let j = 0; j < right.length; j++) {
        const p = left[i], q = left[(i + 1) % left.length];
        const r = right[j], s = right[(j + 1) % right.length];
        if (!segmentsIntersect(p, q, r, s)) continue;
        const sharedEndpoint = samePoint(p, r) || samePoint(p, s) || samePoint(q, r) || samePoint(q, s);
        if (!sharedEndpoint) return false;
      }
    }
    return true;
  }

  function simplifySafely(rawLoops, tolerance, foregroundArea, warnings) {
    const exact = rawLoops.map(removeCollinear);
    const exactPointCount = exact.reduce((sum, loop) => sum + loop.length, 0);
    if (exactPointCount > MAX_EXACT_POINTS)
      fail(`像素轮廓含 ${exactPointCount} 个有效拐点，细节过多；请裁剪、去噪或调整阈值。`);
    if (!(tolerance > 0)) return exact;
    const candidate = exact.map(loop => simplifyClosed(loop, tolerance));
    const candidatePointCount = candidate.reduce((sum, loop) => sum + loop.length, 0);
    // Topology validation below is pairwise. Do not let a photograph with tens of
    // thousands of corners enter quadratic work when it cannot satisfy the output
    // contract in any case.
    if (candidatePointCount > MAX_POINTS)
      fail(`识别轮廓仍有 ${candidatePointCount} 个点，超过 ${MAX_POINTS} 点上限；请裁剪、去噪或提高简化容差。`);
    const signsPreserved = candidate.every((loop, index) =>
      Math.sign(signedArea(loop)) === Math.sign(signedArea(exact[index])));
    const simple = candidate.every(simpleLoop) && loopsDoNotCross(candidate);
    const candidateArea = candidate.reduce((sum, loop) => sum + signedArea(loop), 0);
    const areaLimit = Math.max(1, foregroundArea * 0.02);
    if (!signsPreserved || !simple || Math.abs(candidateArea - foregroundArea) > areaLimit) {
      warnings.push('轮廓简化会改变孔洞、相交关系或面积，已回退到无损直线合并。');
      return exact;
    }
    return candidate;
  }

  function extractContours(image, options) {
    options = options || {};
    const warnings = [];
    const input = validateInput(image, options);
    const segmentation = makeMask(input, warnings);
    const chosen = selectComponents(segmentation, input, warnings);
    let filled = { holes: 0, area: 0 };
    if (input.fillHoles) {
      filled = fillEnclosedBackground(chosen.mask, input.width, input.height);
      if (filled.holes)
        warnings.push(`已按显式设置填充 ${filled.holes} 个内孔，共 ${filled.area} 像素；输出轮廓将忽略这些孔。`);
    }
    const foregroundArea = chosen.selectedPixels + filled.area;
    const rawLoops = tracePixelBoundaries(chosen.mask, input.width, input.height);
    const rawArea = rawLoops.reduce((sum, loop) => sum + signedArea(loop), 0);
    if (Math.abs(rawArea - foregroundArea) > 1e-9)
      fail('像素轮廓面积校验失败；请调整阈值或去噪。');
    const loops = simplifySafely(rawLoops, input.simplifyTolerance, foregroundArea, warnings);
    const pointCount = loops.reduce((sum, loop) => sum + loop.length, 0);
    if (pointCount > MAX_POINTS)
      fail(`识别轮廓仍有 ${pointCount} 个点，超过 ${MAX_POINTS} 点上限；请裁剪、去噪或提高简化容差。`);
    return {
      loops,
      width: input.width,
      height: input.height,
      threshold: segmentation.threshold,
      mode: segmentation.mode,
      stats: {
        mode: segmentation.mode,
        threshold: segmentation.threshold,
        input_pixels: input.width * input.height,
        foreground_pixels: foregroundArea,
        initial_foreground_pixels: segmentation.foregroundPixels,
        components: chosen.components.length,
        eligible_components: chosen.keptCount,
        selected_components: chosen.selected.length,
        removed_components: chosen.removed.length,
        removed_area_px: chosen.removedPixels,
        filled_holes: filled.holes,
        filled_area_px: filled.area,
        loops: loops.length,
        points: pointCount,
        raw_points: rawLoops.reduce((sum, loop) => sum + loop.length, 0),
        simplify_tolerance_px: input.simplifyTolerance,
        background_rgb: segmentation.backgroundRgb
      },
      warnings
    };
  }

  return { extractContours };
}));
