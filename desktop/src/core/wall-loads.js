'use strict';

const DEFINITION = 'fluid-on-wall / density / depth; shared-face pressure and selected viscous flux; torque positive counterclockwise';
const VALUES = ['length','pressureForceX','pressureForceY','viscousForceX','viscousForceY',
  'forceX','forceY','pressureTorque','viscousTorque','torque'];
function validateWallLoads(summary) {
  const keys = ['namedWallLoads','wallLoadReference','wallLoadDefinition'];
  if (!keys.some(k => summary[k] !== undefined)) return; // Legacy result.
  const fail = () => { throw new Error('命名壁面受力或力矩数据不一致。'); };
  const near = (a,b,scale=1) => Number.isFinite(a) && Number.isFinite(b)
    && Math.abs(a-b) <= 1e-10*Math.max(1,Math.abs(a),Math.abs(b),scale);
  if (summary.case !== 'custom' || summary.wallLoadDefinition !== DEFINITION
      || !Array.isArray(summary.wallLoadReference) || summary.wallLoadReference.length !== 2
      || summary.wallLoadReference.some(x => x !== 0) || !Array.isArray(summary.namedWallLoads)
      || !Array.isArray(summary.boundaryConditions)) fail();
  const expected = new Map();
  for (const b of summary.boundaryConditions) if (['wall','moving-wall'].includes(b.type))
    expected.set(b.name,(expected.get(b.name)||0)+1);
  const seen = new Set();
  for (const load of summary.namedWallLoads) {
    if (!load || typeof load !== 'object' || seen.has(load.name) || !expected.has(load.name)
        || load.faces !== expected.get(load.name) || Object.keys(load).length !== VALUES.length+2
        || !VALUES.every(k => typeof load[k] === 'number' && Number.isFinite(load[k])) || !(load.length>0)) fail();
    seen.add(load.name);
    for (const [sum,a,b] of [['forceX','pressureForceX','viscousForceX'],['forceY','pressureForceY','viscousForceY'],
      ['torque','pressureTorque','viscousTorque']])
      if (!near(load[sum],load[a]+load[b],Math.abs(load[a])+Math.abs(load[b]))) fail();
  }
  if (seen.size !== expected.size) fail();
  for (const [total,part] of [['wallForceX','forceX'],['wallForceY','forceY'],
    ['wallViscousForceX','viscousForceX'],['wallViscousForceY','viscousForceY']]) {
    const sum=summary.namedWallLoads.reduce((s,p)=>s+p[part],0);
    const scale=summary.namedWallLoads.reduce((s,p)=>s+Math.abs(p[part]),0);
    if (!near(summary[total],sum,scale)) fail();
  }
}
module.exports = { validateWallLoads, DEFINITION };
