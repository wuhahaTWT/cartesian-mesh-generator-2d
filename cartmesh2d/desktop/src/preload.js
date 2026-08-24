const { contextBridge, ipcRenderer } = require('electron');

contextBridge.exposeInMainWorld('cartmesh', {
  selectDxf: () => ipcRenderer.invoke('select-dxf'),
  selectOutput: () => ipcRenderer.invoke('select-output'),
  examplePath: () => ipcRenderer.invoke('example-path'),
  smokeConfig: () => ipcRenderer.invoke('smoke-config'),
  smokeCapture: () => ipcRenderer.invoke('smoke-capture'),
  generate: request => ipcRenderer.invoke('generate-mesh', request),
  openPath: target => ipcRenderer.invoke('open-path', target),
  onGenerationLine: callback => {
    ipcRenderer.removeAllListeners('generation-line');
    ipcRenderer.on('generation-line', (_event, line) => callback(line));
  }
});
