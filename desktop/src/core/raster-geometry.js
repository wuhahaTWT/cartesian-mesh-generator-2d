'use strict';

const UNITS = { m: 1, cm: .01, mm: .001 };

function calibrateRasterLoops({ loops, pixelWidth, pixelHeight, calibration }) {
  if (!Number.isInteger(pixelWidth) || !Number.isInteger(pixelHeight) || pixelWidth < 2 || pixelHeight < 2 ||
      pixelWidth > 2048 || pixelHeight > 2048 || pixelWidth * pixelHeight > 2097152)
    throw new Error('图片处理尺寸无效，请重新导入。');
  if (!Array.isArray(loops) || !loops.length || loops.length > 64)
    throw new Error('需要至少一个完整轮廓；复杂图片请先裁剪。');
  let count = 0, minX = Infinity, maxX = -Infinity, minY = Infinity, maxY = -Infinity;
  const edges = [];
  for (const loop of loops) {
    if (!Array.isArray(loop) || loop.length < 3) throw new Error('轮廓至少需要三个点。');
    count += loop.length;
    if (count > 4096) throw new Error('轮廓过于复杂，请裁剪图片或调整识别设置。');
    let area = 0;
    for (let i = 0; i < loop.length; i++) {
      const p = loop[i], q = loop[(i + 1) % loop.length];
      for (const point of [p, q]) {
        if (!Array.isArray(point) || point.length !== 2 || !point.every(Number.isFinite) ||
            point[0] < 0 || point[0] > pixelWidth || point[1] < 0 || point[1] > pixelHeight)
          throw new Error('轮廓坐标不在图片范围内。');
      }
      if (Math.hypot(p[0] - q[0], p[1] - q[1]) < 1e-8) throw new Error('轮廓包含重复点。');
      area += p[0] * q[1] - q[0] * p[1];
      minX = Math.min(minX, p[0]); maxX = Math.max(maxX, p[0]);
      minY = Math.min(minY, p[1]); maxY = Math.max(maxY, p[1]);
      edges.push({ p, q, loop, index: i, minX: Math.min(p[0], q[0]), maxX: Math.max(p[0], q[0]),
        minY: Math.min(p[1], q[1]), maxY: Math.max(p[1], q[1]) });
    }
    if (Math.abs(area) < 1e-8) throw new Error('轮廓面积为零。');
  }
  // Reject intersections/touches rather than silently changing the identified shape.
  const cross = (a, b, p) => (b[0] - a[0]) * (p[1] - a[1]) - (b[1] - a[1]) * (p[0] - a[0]);
  const sign = x => Math.abs(x) < 1e-8 ? 0 : Math.sign(x);
  for (let i = 0; i < edges.length; i++) for (let j = i + 1; j < edges.length; j++) {
    const a = edges[i], b = edges[j];
    if (a.loop === b.loop && (j === i + 1 || Math.abs(a.index - b.index) === a.loop.length - 1)) continue;
    if (a.maxX < b.minX || b.maxX < a.minX || a.maxY < b.minY || b.maxY < a.minY) continue;
    if (sign(cross(a.p, a.q, b.p)) * sign(cross(a.p, a.q, b.q)) <= 0 &&
        sign(cross(b.p, b.q, a.p)) * sign(cross(b.p, b.q, a.q)) <= 0)
      throw new Error('识别轮廓存在交叉或接触，请调整阈值或裁剪图片。');
  }
  const width = Number(calibration?.width), unit = calibration?.unit;
  const physicalWidth = width * UNITS[unit];
  if (!Number.isFinite(physicalWidth) || !(physicalWidth > 0) || !(maxX > minX))
    throw new Error('请填写轮廓的实际宽度，并选择 m、cm 或 mm。');
  const scale = physicalWidth / (maxX - minX);
  const cx = (minX + maxX) / 2, cy = (minY + maxY) / 2;
  const converted = loops.map(loop => loop.map(([x, y]) => [(x - cx) * scale, (cy - y) * scale]));
  if (!converted.flat().flat().every(Number.isFinite)) throw new Error('标定尺寸超出可处理范围。');
  return { loops: converted, metresPerPixel: scale, physicalWidth, physicalHeight: (maxY - minY) * scale,
    pixelBounds: { minX, maxX, minY, maxY }, pointCount: count, outputUnits: 'm' };
}
// Check encoded dimensions before asking Chromium to decode potentially huge images.
function inspectRaster(bytes) {
  let format, width, height;
  if (bytes.length >= 24 && bytes.subarray(0,8).equals(Buffer.from([137,80,78,71,13,10,26,10])) && bytes.toString('ascii',12,16) === 'IHDR') {
    format = 'png'; width = bytes.readUInt32BE(16); height = bytes.readUInt32BE(20);
  } else if (bytes.length >= 4 && bytes[0] === 255 && bytes[1] === 216) {
    format = 'jpeg'; let i = 2;
    while (i + 3 < bytes.length) {
      if (bytes[i++] !== 255) break;
      while (bytes[i] === 255) i++;
      const marker = bytes[i++];
      if (marker === 217 || marker === 218) break;
      if (marker === 1 || (marker >= 208 && marker <= 216)) continue;
      if (i + 2 > bytes.length) break;
      const size = bytes.readUInt16BE(i);
      if (size < 2 || i + size > bytes.length) break;
      if ([192,193,194,195,197,198,199,201,202,203,205,206,207].includes(marker) && size >= 8) {
        height = bytes.readUInt16BE(i + 3); width = bytes.readUInt16BE(i + 5); break;
      }
      i += size;
    }
  }
  if (!width || !height || !format) throw new Error('文件内容不是可识别的 PNG / JPEG 图片。');
  if (width > 20000 || height > 20000 || width * height > 40000000)
    throw new Error('图片超过 4000 万像素或单边 20000 像素，请先裁剪或缩小。');
  return { format, width, height };
}
module.exports = { calibrateRasterLoops, inspectRaster };
