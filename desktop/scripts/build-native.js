#!/usr/bin/env node
'use strict';

const childProcess = require('node:child_process');
const fs = require('node:fs');
const path = require('node:path');
const {
  SAMPLE_SOURCES,
  SCHEMA_VERSION,
  TOOL_NAMES,
  executableFile,
  expectedFormat,
  inspectBinary,
  verifyRuntime
} = require('./runtime');

const desktopDir = path.resolve(__dirname, '..');
const projectDir = path.resolve(desktopDir, '..');
const buildDir = path.resolve(process.env.CARTMESH2D_BUILD_DIR || path.join(projectDir, 'build'));
const runtimeDir = path.join(desktopDir, 'runtime');
const jobsText = process.env.CARTMESH2D_BUILD_JOBS || '2';

if (!/^\d+$/.test(jobsText) || Number(jobsText) < 1) {
  throw new Error(`CARTMESH2D_BUILD_JOBS must be a positive integer, got ${jobsText}`);
}
if (!['darwin', 'linux', 'win32'].includes(process.platform)) {
  throw new Error(`unsupported build host: ${process.platform}`);
}
if (!['x64', 'arm64', 'ia32'].includes(process.arch)) {
  throw new Error(`unsupported build architecture: ${process.arch}`);
}

function run(command, args) {
  console.log(`> ${command} ${args.join(' ')}`);
  const result = childProcess.spawnSync(command, args, {
    cwd: projectDir,
    env: process.env,
    shell: false,
    stdio: 'inherit'
  });
  if (result.error) throw result.error;
  if (result.status !== 0) {
    throw new Error(`${command} failed with exit code ${result.status}`);
  }
}

function configureArguments() {
  const args = ['-S', projectDir, '-B', buildDir, '-DCMAKE_BUILD_TYPE=Release'];
  if (process.platform === 'darwin') {
    const compiler = process.env.CARTMESH2D_CXX || '/usr/bin/clang++';
    if (!fs.existsSync(compiler)) throw new Error(`macOS C++ compiler not found: ${compiler}`);
    args.push(`-DCMAKE_CXX_COMPILER=${compiler}`);
  }
  if (process.platform === 'win32') {
    const cmakeArch = { x64: 'x64', arm64: 'ARM64', ia32: 'Win32' }[process.arch];
    args.push('-A', cmakeArch);
  }
  return args;
}

function findBuiltExecutable(tool) {
  const filename = executableFile(tool);
  const candidates = process.platform === 'win32'
    ? [path.join(buildDir, 'Release', filename), path.join(buildDir, filename)]
    : [path.join(buildDir, filename), path.join(buildDir, 'Release', filename)];
  const found = candidates.find(candidate => fs.existsSync(candidate));
  if (!found) throw new Error(`CMake did not produce ${filename} in ${buildDir}`);
  return found;
}

function verifyMacDependencies(file) {
  if (process.platform !== 'darwin') return;
  const result = childProcess.spawnSync('otool', ['-L', file], {
    encoding: 'utf8',
    shell: false
  });
  if (result.error) throw result.error;
  if (result.status !== 0) throw new Error(`otool failed for ${file}`);
  if (result.stdout.includes('@rpath/lib/libstdc++')) {
    throw new Error(`${file} links a non-system libstdc++; it would not run inside the app`);
  }
}

function prepareRuntime() {
  const stagingDir = fs.mkdtempSync(path.join(desktopDir, '.runtime-staging-'));
  const binDir = path.join(stagingDir, 'bin');
  const samplesDir = path.join(stagingDir, 'samples');
  fs.mkdirSync(binDir);
  fs.mkdirSync(samplesDir);

  try {
    const executables = TOOL_NAMES.map(tool => {
      const source = findBuiltExecutable(tool);
      const file = executableFile(tool);
      const target = path.join(binDir, file);
      const inspected = inspectBinary(source);
      const requiredFormat = expectedFormat(process.platform);
      if (inspected.format !== requiredFormat || !inspected.architectures.includes(process.arch)) {
        throw new Error(
          `${source} is ${inspected.format}/${inspected.architectures.join('+')}, expected ${requiredFormat}/${process.arch}`
        );
      }
      fs.copyFileSync(source, target);
      if (process.platform !== 'win32') fs.chmodSync(target, 0o755);
      verifyMacDependencies(target);
      return { name: tool, file, format: inspected.format, architectures: inspected.architectures };
    });

    for (const sourceName of SAMPLE_SOURCES) {
      const source = path.join(projectDir, 'examples', ...sourceName.split('/'));
      const target = path.join(samplesDir, path.basename(sourceName));
      if (!fs.existsSync(source)) throw new Error(`sample source is missing: ${source}`);
      fs.copyFileSync(source, target);
    }

    const manifest = {
      schemaVersion: SCHEMA_VERSION,
      platform: process.platform,
      arch: process.arch,
      executables,
      samples: SAMPLE_SOURCES.map(source => path.basename(source))
    };
    fs.writeFileSync(path.join(stagingDir, 'manifest.json'), `${JSON.stringify(manifest, null, 2)}\n`);

    fs.rmSync(runtimeDir, { recursive: true, force: true });
    fs.renameSync(stagingDir, runtimeDir);
    verifyRuntime({ runtimeDir });
  } catch (error) {
    fs.rmSync(stagingDir, { recursive: true, force: true });
    throw error;
  }
}

run('cmake', configureArguments());
run('cmake', [
  '--build', buildDir,
  '--config', 'Release',
  '--target', ...TOOL_NAMES,
  '--parallel', jobsText
]);
prepareRuntime();
console.log(`prepared native runtime for ${process.platform}/${process.arch}: ${runtimeDir}`);
