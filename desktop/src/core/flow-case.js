'use strict';

const { createHash } = require('node:crypto');
const { parseCm2d } = require('./cm2d');
const { validateFlowRequest } = require('./flow');
const { validateBoundaryMesh } = require('./flow-boundaries');

const FORMAT = 'cartmesh2d-flow-case-v1';
const MAX_BYTES = 32 * 1024 * 1024;
const fail = message => { throw new Error(`流动工况：${message}`); };
function keys(value, allowed) {
  if (!value || typeof value !== 'object' || Array.isArray(value) ||
      Object.keys(value).some(key => !allowed.includes(key))) fail('存在未知或不支持的字段。');
}
function meshIdentity(source) {
  if (!(typeof source === 'string' || Buffer.isBuffer(source))) fail('缺少当前最终网格。');
  const mesh = parseCm2d(source.toString());
  return { mesh, identity: { sha256: createHash('sha256').update(source).digest('hex'),
    cells: mesh.cells.length, faces: mesh.edges.length, vertices: mesh.vertices.length } };
}
function normalize(request, mesh) {
  if (!request || request.resume) fail('保存的是从零起算的设置，请先取消续算选项；已接受状态请保存检查点。');
  const normalized = validateFlowRequest(request);
  if (normalized.case === 'custom')
    normalized.boundaryDefinition = validateBoundaryMesh(normalized.boundaryDefinition, mesh, normalized.speed);
  return normalized;
}
function createFlowCaseDocument(request, meshSource) {
  const { mesh, identity } = meshIdentity(meshSource);
  return { format: FORMAT, mesh: identity, request: normalize(request, mesh) };
}
function serializeFlowCase(document) {
  const text = JSON.stringify(document, null, 2) + '\n';
  if (Buffer.byteLength(text) > MAX_BYTES) fail('文件超过32MB。');
  return text;
}
function parseFlowCaseDocument(text, meshSource) {
  if (typeof text !== 'string' || Buffer.byteLength(text) > MAX_BYTES) fail('文件内容无效或超过32MB。');
  let document;
  try { document = JSON.parse(text); } catch { fail('不是有效的JSON文件。'); }
  keys(document, ['format', 'mesh', 'request']);
  if (document.format !== FORMAT) fail('不支持的文件版本。');
  keys(document.mesh, ['sha256', 'cells', 'faces', 'vertices']);
  const { mesh, identity } = meshIdentity(meshSource);
  if (Object.entries(identity).some(([key, value]) => document.mesh[key] !== value))
    fail('与当前最终网格不匹配。请先生成保存该工况时的同一网格。');
  const request = normalize(document.request, mesh);
  // A saved case is canonical: refusing extra/inactive fields prevents silently
  // accepting settings from a newer implementation and then ignoring them.
  keys(document.request, Object.keys(request));
  if (Object.keys(request).some(key => !Object.hasOwn(document.request, key))) fail('缺少求解参数。');
  for (const [key, value] of Object.entries(request)) {
    if (key === 'boundaryDefinition' || key === 'initialVortex') continue;
    if (document.request[key] !== value) fail('求解参数类型无效。');
  }
  if (request.initialVortex) keys(document.request.initialVortex, ['centre', 'radius', 'peakSpeed']);
  if (request.boundaryDefinition) {
    keys(document.request.boundaryDefinition, ['cells', 'faces', 'records']);
    for (const row of document.request.boundaryDefinition.records)
      keys(row, ['face', 'owner', 'x', 'y', 'sx', 'sy', 'type', 'name', 'u', 'v', 'p']);
  }
  return { format: FORMAT, mesh: identity, request };
}

module.exports = { MAX_BYTES, createFlowCaseDocument, serializeFlowCase, parseFlowCaseDocument };
