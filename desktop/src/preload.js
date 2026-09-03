const { contextBridge, ipcRenderer } = require('electron');

contextBridge.exposeInMainWorld('cartmesh', {
  selectDxf: () => ipcRenderer.invoke('select-dxf'),
  selectOutput: () => ipcRenderer.invoke('select-output'),
  examplePath: () => ipcRenderer.invoke('example-path'),
  smokeConfig: () => ipcRenderer.invoke('smoke-config'),
  smokeCapture: () => ipcRenderer.invoke('smoke-capture'),
  capabilities: () => ipcRenderer.invoke('mesh-capabilities'),
  estimate: request => ipcRenderer.invoke('estimate-mesh-job', request),
  generate: request => ipcRenderer.invoke('generate-mesh', request),
  openPath: target => ipcRenderer.invoke('open-path', target),
  onGenerationLine: callback => {
    ipcRenderer.removeAllListeners('generation-line');
    ipcRenderer.on('generation-line', (_event, line) => callback(line));
  }
});
