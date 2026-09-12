'use strict';
const test = require('node:test');
const assert = require('node:assert/strict');
const { calibrateRasterLoops } = require('../src/core/raster-geometry');
const fixture = () => ({pixelWidth:100,pixelHeight:80,loops:[[[10,20],[90,20],[90,60],[10,60]]],calibration:{width:20,unit:'mm'}});
test('raster calibration uses extracted width, flips image y and preserves aspect', () => {
 const input=fixture(), original=JSON.stringify(input), result=calibrateRasterLoops(input);
 assert.equal(result.physicalWidth,.02);assert.equal(result.physicalHeight,.01);
 assert.deepEqual(result.loops[0][0],[-.01,.005]);assert.equal(result.metresPerPixel,.00025);
 assert.equal(JSON.stringify(input),original);
});
test('raster requires meaningful calibration and rejects self-intersection and contour contact', () => {
 for(const width of [0,-1,NaN,Infinity]) assert.throws(()=>calibrateRasterLoops({...fixture(),calibration:{width,unit:'m'}}));
 assert.throws(()=>calibrateRasterLoops({...fixture(),loops:[[[10,10],[90,60],[10,60],[80,10]]]}),/交叉/);
 const f=fixture();f.loops.push([[40,20],[50,30],[30,30]]);
 assert.throws(()=>calibrateRasterLoops(f),/交叉|接触/);
 assert.throws(()=>calibrateRasterLoops({...fixture(),loops:[[[10,20],[110,20],[90,60]]]}),/范围/);
});
test('raster calibration keeps independent components and interior holes', () => {
 const f=fixture();f.loops.push([[30,30],[30,40],[50,40],[50,30]]);
 assert.equal(calibrateRasterLoops(f).loops.length,2);
});
test('encoded image dimensions are bounded before decode and JPEG metadata is scanned', () => {
 const {inspectRaster}=require('../src/core/raster-geometry');
 const png=Buffer.alloc(24);Buffer.from([137,80,78,71,13,10,26,10]).copy(png);png.write('IHDR',12);png.writeUInt32BE(100,16);png.writeUInt32BE(80,20);
 assert.deepEqual(inspectRaster(png),{format:'png',width:100,height:80});
 png.writeUInt32BE(100000,16);assert.throws(()=>inspectRaster(png),/超过/);
 const jpg=Buffer.from([255,216,255,224,0,4,0,0,255,194,0,8,8,0,80,0,100,1]);
 assert.deepEqual(inspectRaster(jpg),{format:'jpeg',width:100,height:80});
 assert.throws(()=>inspectRaster(Buffer.from('<svg/>')),/不是/);
 assert.throws(()=>inspectRaster(jpg.subarray(0,15)),/不是/);
});
