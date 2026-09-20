'use strict';
const test=require('node:test');
const assert=require('node:assert/strict');
const {validateWallLoads,DEFINITION}=require('../src/core/wall-loads');
const sample=()=>({case:'custom',wallLoadReference:[0,0],wallLoadDefinition:DEFINITION,
  boundaryConditions:[{face:1,type:'moving-wall',name:'转子'},{face:2,type:'wall',name:'壳体'},
    {face:3,type:'velocity-inlet',name:'inlet'}],
  namedWallLoads:[{name:'转子',faces:1,length:2,pressureForceX:1,pressureForceY:2,viscousForceX:3,viscousForceY:4,
    forceX:4,forceY:6,pressureTorque:.5,viscousTorque:-1,torque:-.5},
    {name:'壳体',faces:1,length:3,pressureForceX:-1,pressureForceY:-2,viscousForceX:-3,viscousForceY:-4,
    forceX:-4,forceY:-6,pressureTorque:-.5,viscousTorque:1,torque:.5}],
  wallForceX:0,wallForceY:0,wallViscousForceX:0,wallViscousForceY:0});
test('wall loads retain signed patch forces and moments with explicit units/reference',()=>{
  const result=sample();validateWallLoads(result);assert.equal(result.namedWallLoads[0].torque,-.5);
  validateWallLoads({case:'custom'});
});
test('wall load validation rejects missing patches, altered reference, force and moment mismatches',()=>{
  for(const mutate of [s=>s.namedWallLoads.pop(),s=>s.namedWallLoads.push({...s.namedWallLoads[0]}),
    s=>s.wallLoadReference=[1,0],s=>delete s.wallLoadDefinition,s=>s.namedWallLoads[0].faces++,
    s=>s.namedWallLoads[0].torque=.5,s=>s.namedWallLoads[0].forceX=5,s=>s.wallForceX=1,
    s=>s.namedWallLoads[0].viscousTorque=NaN,s=>s.namedWallLoads[0].length=0]){
    const value=sample();mutate(value);assert.throws(()=>validateWallLoads(value),/壁面受力/);
  }
});
