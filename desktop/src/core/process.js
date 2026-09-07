'use strict';
const { spawn } = require('node:child_process');

function run(command, args, onLine = () => {}, signal) {
  return new Promise((resolve, reject) => {
    if (signal?.aborted) return reject(new Error('操作已取消'));
    const child = spawn(command, args, { windowsHide: true });
    let killTimer;
    const cancel = () => {
      child.kill('SIGTERM');
      killTimer = setTimeout(() => child.kill('SIGKILL'), 1500);
    };
    signal?.addEventListener('abort', cancel, { once: true });
    const cleanup = () => { clearTimeout(killTimer); signal?.removeEventListener('abort', cancel); };
    let stdout = '';
    let stderr = '';
    const consume = (chunk, isError) => {
      const text = chunk.toString();
      if (isError) stderr += text; else stdout += text;
      text.split(/\r?\n/).filter(Boolean).forEach(onLine);
    };
    child.stdout.on('data', chunk => consume(chunk, false));
    child.stderr.on('data', chunk => consume(chunk, true));
    child.on('error', error => { cleanup(); reject(error); });
    child.on('close', code => {
      cleanup();
      if (signal?.aborted) return reject(new Error('操作已取消'));
      // Both CLIs are fail-closed: a non-zero exit means no mesh was committed, and
      // the reason is on stderr.  Surfacing stdout as the fallback keeps the size
      // field's refusal readable even when it printed its diagnosis first.
      if (code === 0) resolve({ stdout, stderr, code });
      else reject(Object.assign(new Error(stderr.trim() || stdout.trim() || `退出码 ${code}`),
                                { stdout, stderr, code }));
    });
  });
}

module.exports = { run };
