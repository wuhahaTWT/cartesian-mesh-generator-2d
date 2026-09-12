'use strict';

const fs = require('node:fs');
const path = require('node:path');

const SCHEMA_VERSION = 1;
const TOOL_NAMES = Object.freeze([
  'cartmesh2d_cli',
  'cartmesh2d_hybrid_cli',
  'cartmesh2d_dxf_cli'
]);

const SAMPLE_SOURCES = Object.freeze([
  'acceptance/circle.xy',
  'complex/naca2412_dense.xy',
  'complex/thick_cambered_airfoil.xy',
  'complex/two_obstacles.xy',
  'complex/superellipse_24.xy',
  'complex/serpentine_body.xy',
  'complex/gear_star.xy',
  'complex/annulus.xy',
  'complex/nozzle_profile.xy',
  'h4_3/narrow_gap.xy',
  'h4_3/sharp_trailing_edge.xy'
]);

const FORMAT_BY_PLATFORM = Object.freeze({
  darwin: 'mach-o',
  linux: 'elf',
  win32: 'pe'
});

const archFromCpuType = cpuType => ({
  7: 'ia32',
  0x01000007: 'x64',
  0x0100000c: 'arm64'
})[cpuType];

const archFromElfMachine = machine => ({
  3: 'ia32',
  62: 'x64',
  183: 'arm64'
})[machine];

const archFromPeMachine = machine => ({
  0x014c: 'ia32',
  0x8664: 'x64',
  0xaa64: 'arm64'
})[machine];

function inspectBinary(file) {
  const handle = fs.openSync(file, 'r');
  try {
    const stat = fs.fstatSync(handle);
    const header = Buffer.alloc(Math.min(stat.size, 4096));
    fs.readSync(handle, header, 0, header.length, 0);

    if (header.length >= 20 && header.subarray(0, 4).equals(Buffer.from([0x7f, 0x45, 0x4c, 0x46]))) {
      const littleEndian = header[5] === 1;
      const machine = littleEndian ? header.readUInt16LE(18) : header.readUInt16BE(18);
      return { format: 'elf', architectures: [archFromElfMachine(machine) || `elf-machine-${machine}`] };
    }

    if (header.length >= 64 && header[0] === 0x4d && header[1] === 0x5a) {
      const peOffset = header.readUInt32LE(0x3c);
      if (peOffset + 6 > header.length) {
        throw new Error(`${file}: PE header is outside the inspected prefix`);
      }
      if (!header.subarray(peOffset, peOffset + 4).equals(Buffer.from('PE\0\0', 'binary'))) {
        throw new Error(`${file}: invalid PE signature`);
      }
      const machine = header.readUInt16LE(peOffset + 4);
      return { format: 'pe', architectures: [archFromPeMachine(machine) || `pe-machine-${machine}`] };
    }

    if (header.length >= 8) {
      const magic = header.readUInt32BE(0);
      const thinMach = new Map([
        [0xfeedface, 'be'],
        [0xfeedfacf, 'be'],
        [0xcefaedfe, 'le'],
        [0xcffaedfe, 'le']
      ]);
      if (thinMach.has(magic)) {
        const cpuType = thinMach.get(magic) === 'le' ? header.readUInt32LE(4) : header.readUInt32BE(4);
        return { format: 'mach-o', architectures: [archFromCpuType(cpuType) || `mach-cpu-${cpuType}`] };
      }

      const fatMach = new Map([
        [0xcafebabe, { endian: 'be', stride: 20 }],
        [0xcafebabf, { endian: 'be', stride: 32 }],
        [0xbebafeca, { endian: 'le', stride: 20 }],
        [0xbfbafeca, { endian: 'le', stride: 32 }]
      ]);
      if (fatMach.has(magic)) {
        const { endian, stride } = fatMach.get(magic);
        const read32 = offset => endian === 'le' ? header.readUInt32LE(offset) : header.readUInt32BE(offset);
        const count = read32(4);
        if (count < 1 || 8 + count * stride > header.length) {
          throw new Error(`${file}: invalid Mach-O universal header`);
        }
        const architectures = [];
        for (let index = 0; index < count; index += 1) {
          const cpuType = read32(8 + index * stride);
          architectures.push(archFromCpuType(cpuType) || `mach-cpu-${cpuType}`);
        }
        return { format: 'mach-o', architectures: [...new Set(architectures)] };
      }
    }

    throw new Error(`${file}: unsupported native executable format`);
  } finally {
    fs.closeSync(handle);
  }
}

function executableFile(tool, platform = process.platform) {
  return platform === 'win32' ? `${tool}.exe` : tool;
}

function expectedFormat(platform) {
  const format = FORMAT_BY_PLATFORM[platform];
  if (!format) throw new Error(`unsupported platform: ${platform}`);
  return format;
}

function ensureRegularFile(file, label) {
  let stat;
  try {
    stat = fs.statSync(file);
  } catch (error) {
    throw new Error(`${label} is missing: ${file}`);
  }
  if (!stat.isFile()) throw new Error(`${label} is not a regular file: ${file}`);
  return stat;
}

function verifyRuntime({ runtimeDir, expectedPlatform = process.platform, expectedArch = process.arch }) {
  if (expectedPlatform !== process.platform) {
    throw new Error(
      `refusing cross-platform packaging on ${process.platform}; ${expectedPlatform} packages must be built on ${expectedPlatform}`
    );
  }
  if (expectedArch !== process.arch) {
    throw new Error(
      `refusing cross-architecture packaging on ${process.arch}; requested runtime architecture is ${expectedArch}`
    );
  }

  const manifestPath = path.join(runtimeDir, 'manifest.json');
  ensureRegularFile(manifestPath, 'runtime manifest');
  const manifest = JSON.parse(fs.readFileSync(manifestPath, 'utf8'));
  if (manifest.schemaVersion !== SCHEMA_VERSION) {
    throw new Error(`unsupported runtime manifest schema: ${manifest.schemaVersion}`);
  }
  if (manifest.platform !== expectedPlatform || manifest.arch !== expectedArch) {
    throw new Error(
      `runtime is for ${manifest.platform}/${manifest.arch}, expected ${expectedPlatform}/${expectedArch}; run npm run build:native on this host`
    );
  }

  const expectedExecutables = TOOL_NAMES.map(tool => executableFile(tool, expectedPlatform));
  const listedExecutables = manifest.executables.map(entry => entry.file);
  if (JSON.stringify(listedExecutables) !== JSON.stringify(expectedExecutables)) {
    throw new Error(`runtime executable list does not match: ${listedExecutables.join(', ')}`);
  }

  const format = expectedFormat(expectedPlatform);
  for (const entry of manifest.executables) {
    const file = path.join(runtimeDir, 'bin', entry.file);
    const stat = ensureRegularFile(file, `runtime executable ${entry.file}`);
    if (expectedPlatform !== 'win32' && (stat.mode & 0o111) === 0) {
      throw new Error(`runtime executable is not executable: ${file}`);
    }
    const inspected = inspectBinary(file);
    if (entry.format !== inspected.format || inspected.format !== format) {
      throw new Error(`${entry.file} has ${inspected.format} format, expected ${format}`);
    }
    if (!inspected.architectures.includes(expectedArch)) {
      throw new Error(`${entry.file} is for ${inspected.architectures.join('/')}, expected ${expectedArch}`);
    }
    if (JSON.stringify(entry.architectures) !== JSON.stringify(inspected.architectures)) {
      throw new Error(`${entry.file} no longer matches its runtime manifest architecture`);
    }
  }

  const expectedSamples = SAMPLE_SOURCES.map(source => path.basename(source));
  if (JSON.stringify(manifest.samples) !== JSON.stringify(expectedSamples)) {
    throw new Error('runtime sample list does not match the desktop sample contract');
  }
  for (const sample of expectedSamples) {
    ensureRegularFile(path.join(runtimeDir, 'samples', sample), `runtime sample ${sample}`);
  }
  const sampleModule = path.join(path.dirname(runtimeDir), 'src', 'core', 'samples.js');
  ensureRegularFile(sampleModule, 'desktop sample catalog');
  const catalogSamples = require(sampleModule).SAMPLES.map(sample => sample.file);
  if (new Set(catalogSamples).size !== catalogSamples.length
      || JSON.stringify([...catalogSamples].sort()) !== JSON.stringify([...expectedSamples].sort())) {
    throw new Error('runtime sample list does not match desktop/src/core/samples.js');
  }

  return manifest;
}

module.exports = {
  FORMAT_BY_PLATFORM,
  SAMPLE_SOURCES,
  SCHEMA_VERSION,
  TOOL_NAMES,
  executableFile,
  expectedFormat,
  inspectBinary,
  verifyRuntime
};
