'use strict';

const fs = require('node:fs');
const fsp = fs.promises;
const path = require('node:path');
const readline = require('node:readline');

const SUPPORTED_CASES = new Set(['external', 'channel', 'cavity']);
const SUPPORTED_CONVECTION = new Set(['upwind', 'limited-linear']);
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
  try {
    for await (const rawLine of input) {
      const line = rawLine.trim();
      if (line === '') continue;
      if (!metadata) {
        if (line !== 'CARTMESH2D_FLOW_CHECKPOINT 1' && line !== 'CARTMESH2D_FLOW_CHECKPOINT 2') invalid('unsupported header');
        version = line.endsWith(' 2') ? 2 : 1;
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
      if (!cells) {
        const parts = fields(line);
        if (parts[0] === 'CELLS') {
          if (parts.length !== 2) invalid('CELLS metadata is malformed');
          countToken(parts[1], 'cells');
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
          countToken(parts[1], 'faces');
          faces = true;
        }
        continue;
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
