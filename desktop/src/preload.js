const { contextBridge, ipcRenderer } = require('electron');

contextBridge.exposeInMainWorld('cartmesh', {
  catalog: () => ipcRenderer.invoke('catalog'),
  pickGeometry: () => ipcRenderer.invoke('pick-geometry'),
  exportResult: () => ipcRenderer.invoke('export-result'),
  cancel: () => ipcRenderer.invoke('cancel'),
  previewGeometry: request => ipcRenderer.invoke('preview-geometry', request),
  probeSizing: request => ipcRenderer.invoke('probe-sizing', request),
  generate: request => ipcRenderer.invoke('generate', request),
  openPath: target => ipcRenderer.invoke('open-path', target),
  onRunLine: callback => {
    ipcRenderer.removeAllListeners('run-line');
    ipcRenderer.on('run-line', (_event, line) => callback(line));
  }
});
