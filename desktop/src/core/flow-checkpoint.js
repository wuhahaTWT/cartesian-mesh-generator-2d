'use strict';

const fs = require('node:fs');
const fsp = fs.promises;
const path = require('node:path');
const readline = require('node:readline');
const { quotedTokens, condition, normalizeBoundaryDefinition } = require('./flow-boundaries');

const SUPPORTED_CASES = new Set(['external', 'channel', 'duct', 'cavity', 'custom']);
const SUPPORTED_CONVECTION = new Set(['upwind', 'limited-linear', 'face-limited-linear']);
const SUPPORTED_STRESS = new Set(['symmetric']);
const SUPPORTED_OUTLET_BACKFLOW = new Set(['reject', 'normal-inlet']);

function invalid(message) {
  throw new Error(`Invalid flow checkpoint metadata: ${message}`);
}

function numberToken(value, name, { positive = false } = {}) {
  if (!value || !/^[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?$/.test(value)) {
    invalid(`${name} is not a finite number`);
  }
  const number = Number(value);
  if (!Number.isFinite(number) || (positive && number <= 0)) {
    invalid(`${name} must be finite and positive`);
  }
  return number;
}

function countToken(value, name) {
  if (!/^\d+$/.test(value)) invalid(`${name} count is invalid`);
  const count = Number(value);
  if (!Number.isSafeInteger(count) || count <= 0) invalid(`${name} count must be positive`);
  return count;
}

function fields(line) {
  return line.trim().split(/\s+/);
}

function parseConfig(line, version) {
  // std::quoted output uses backslash escaping. The supported case names do not
  // contain escapes, but parsing them explicitly keeps malformed metadata out.
  const match = /^CONFIG\s+"((?:[^"\\]|\\.)*)"\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)\s+(\S+)(?:\s+(\S+))?$/.exec(line.trim());
  if (!match) invalid('CONFIG line is malformed');
  const scenario = match[1];
  if (!SUPPORTED_CASES.has(scenario)) invalid(`unsupported case '${scenario}'`);
  const nu = numberToken(match[2], 'nu', { positive: true });
  const speed = numberToken(match[3], 'speed', { positive: true });
  const convection = match[4];
  const viscousStress = match[5];
  const slope = numberToken(match[6], 'pressure slope');
  if (!SUPPORTED_CONVECTION.has(convection)) invalid(`unsupported convection '${convection}'`);
  if (!SUPPORTED_STRESS.has(viscousStress)) invalid(`unsupported viscous stress '${viscousStress}'`);
  if (slope !== 0) invalid('manufactured/nonzero pressure slope is unsupported');
  const outletBackflow = match[7] === undefined ? 'reject' : match[7];
  if (version === 1 && match[7] !== undefined) invalid('legacy v1 CONFIG has unexpected outlet backflow mode');
  if (version >= 2 && match[7] === undefined) invalid('missing outlet backflow mode');
  if (!SUPPORTED_OUTLET_BACKFLOW.has(outletBackflow)) invalid(`unsupported outlet backflow mode '${outletBackflow}'`);
  if ((version === 4) !== (scenario === 'custom') || (scenario === 'custom' && outletBackflow !== 'reject'))
    invalid('explicit boundaries require v4 and reject backflow');
  return { case: scenario, nu, speed, convection, viscousStress, outletBackflow };
}

async function readCheckpointMetadata(filePath) {
  if (typeof filePath !== 'string' || filePath.length === 0) invalid('path is empty');
  let stat;
  try {
    stat = await fsp.stat(filePath);
  } catch (error) {
    throw new Error(`Invalid flow checkpoint metadata: cannot stat file (${error.message})`);
  }
  if (!stat.isFile()) invalid('path is not a regular file');
  if (stat.size === 0) invalid('file is empty');

  const stream = fs.createReadStream(filePath, { encoding: 'utf8' });
  const input = readline.createInterface({ input: stream, crlfDelay: Infinity });
  let metadata = null, version;
  let cells = false;
  let faces = false;
  let time;
  let viscosity = false, boundaryCount, cellCount, faceCount;
  const boundaries = new Map();
  try {
    for await (const rawLine of input) {
      const line = rawLine.trim();
      if (line === '') continue;
      if (!metadata) {
        if (!/^CARTMESH2D_FLOW_CHECKPOINT [124]$/.test(line)) invalid('unsupported header');
        version = Number(line.at(-1));
        metadata = {};
        continue;
      }
      if (!metadata.discretization) {
        if (line !== 'DISCRETIZATION Euler-RC-v2') invalid('unsupported discretization');
        metadata.discretization = true;
        continue;
      }
      if (!metadata.config) {
        Object.assign(metadata, parseConfig(line, version));
        metadata.config = true;
        continue;
      }
      if (version === 4 && !viscosity) {
        if (line !== 'FACE_VISCOSITY 0') invalid('desktop custom flow requires constant viscosity');
        viscosity = true;
        continue;
      }
      if (version === 4 && boundaryCount === undefined) {
        const parts=fields(line);
        if (parts.length!==2 || parts[0]!=='BOUNDARIES') invalid('missing BOUNDARIES');
        boundaryCount=countToken(parts[1],'boundaries');
        continue;
      }
      if (version === 4 && boundaries.size < boundaryCount) {
        const parts=quotedTokens(line);
        if (parts.length!==7 || parts[0]!=='BOUNDARY') invalid('malformed BOUNDARY');
        const record=condition({face:Number(parts[1]),type:parts[2],name:parts[3],
          u:numberToken(parts[4],'u'),v:numberToken(parts[5],'v'),p:numberToken(parts[6],'p')});
        if (boundaries.has(record.face)) invalid('duplicate boundary face');
        boundaries.set(record.face,record);
        continue;
      }
      if (!cells) {
        const parts = fields(line);
        if (parts[0] === 'CELLS') {
          if (parts.length !== 2) invalid('CELLS metadata is malformed');
          cellCount=countToken(parts[1], 'cells');
          cells = true;
        } else {
          invalid('CELLS metadata is missing');
        }
        continue;
      }
      if (!faces) {
        const parts = fields(line);
        if (parts[0] === 'FACES') {
          if (parts.length !== 2) invalid('FACES metadata is malformed');
          faceCount=countToken(parts[1], 'faces');
          faces = true;
        }
        continue;
      }
      if (version === 4 && line.startsWith('FACE ')) {
        const parts=fields(line), face=Number(parts[1]);
        if (boundaries.has(face)) {
          if (parts.length!==13 || parts[3]!=='-') invalid('custom face geometry is malformed');
          const b=boundaries.get(face);
          if (b.owner !== undefined) invalid('duplicate custom face geometry');
          Object.assign(b,{owner:Number(parts[2]),x:numberToken(parts[5],'face x'),y:numberToken(parts[6],'face y'),
            sx:numberToken(parts[7],'face sx'),sy:numberToken(parts[8],'face sy')});
        }
      }
      if (line.startsWith('TIME ')) {
        const parts = fields(line);
        if (parts.length !== 2) invalid('TIME metadata is malformed');
        time = numberToken(parts[1], 'time');
        if (time < 0) invalid('time must be nonnegative');
        break;
      }
    }
  } catch (error) {
    if (error && error.message && error.message.startsWith('Invalid flow checkpoint metadata:')) throw error;
    throw new Error(`Invalid flow checkpoint metadata: cannot read file (${error.message})`);
  } finally {
    input.close();
    stream.destroy();
  }
  if (!metadata || !metadata.discretization || !metadata.config || !cells || !faces)
    invalid('truncated metadata');
  if (time === undefined) invalid('missing TIME');
  return {
    ...(version === 4 ? {boundaryDefinition:normalizeBoundaryDefinition({cells:cellCount,faces:faceCount,records:[...boundaries.values()]})} : {}),
    case: metadata.case,
    nu: metadata.nu,
    speed: metadata.speed,
    convection: metadata.convection,
    viscousStress: metadata.viscousStress,
    outletBackflow: metadata.outletBackflow,
    time,
    fileName: path.basename(filePath)
  };
}

module.exports = { readCheckpointMetadata };
