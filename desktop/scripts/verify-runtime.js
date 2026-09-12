#!/usr/bin/env node
'use strict';

const path = require('node:path');
const { verifyRuntime } = require('./runtime');

function option(name) {
  const index = process.argv.indexOf(name);
  return index === -1 ? undefined : process.argv[index + 1];
}

try {
  const desktopDir = path.resolve(__dirname, '..');
  const manifest = verifyRuntime({
    runtimeDir: path.join(desktopDir, 'runtime'),
    expectedPlatform: option('--platform') || process.platform,
    expectedArch: option('--arch') || process.arch
  });
  console.log(
    `verified ${manifest.executables.length} native tools and ${manifest.samples.length} samples for ${manifest.platform}/${manifest.arch}`
  );
} catch (error) {
  console.error(error.message);
  process.exitCode = 1;
}
