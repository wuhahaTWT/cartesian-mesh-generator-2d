'use strict';

/* The extraction code is shared with the main process tests.  Its browser build
 * attaches the same API to self.CartMeshRaster. */
importScripts('../core/raster.js');

self.onmessage = event => {
  const { id, image, options } = event.data || {};
  if (!Number.isInteger(id)) return;
  try {
    if (!self.CartMeshRaster || typeof self.CartMeshRaster.extractContours !== 'function') {
      throw new Error('图片轮廓提取模块没有正确加载。');
    }
    const data = image && image.data instanceof Uint8ClampedArray
      ? image.data
      : new Uint8ClampedArray(image && image.data);
    const result = self.CartMeshRaster.extractContours({
      width: image.width,
      height: image.height,
      data
    }, options || {});
    self.postMessage({ id, result });
  } catch (error) {
    self.postMessage({ id, error: error && error.message ? error.message : String(error) });
  }
};
