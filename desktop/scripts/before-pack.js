'use strict';

const path = require('node:path');
const { Arch } = require('builder-util');
const { verifyRuntime } = require('./runtime');

exports.default = async function beforePack(context) {
  const targetArch = Arch[context.arch];
  if (typeof targetArch !== 'string') {
    throw new Error(`electron-builder supplied an unknown target architecture: ${context.arch}`);
  }
  const runtimeDir = path.join(context.packager.projectDir, 'runtime');
  const manifest = verifyRuntime({
    runtimeDir,
    expectedPlatform: context.electronPlatformName,
    expectedArch: targetArch
  });
  console.log(`packaging verified runtime ${manifest.platform}/${manifest.arch}`);
};
