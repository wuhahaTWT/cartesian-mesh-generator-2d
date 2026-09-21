'use strict';
const test=require('node:test');
const assert=require('node:assert/strict');
const {validateBoundaryFluxes,DEFINITION}=require('../src/core/boundary-flux');
const sample=()=>({case:'custom',globalImbalance:0,boundaryFluxDefinition:DEFINITION,
 boundaryConditions:[{face:0,name:'入口',type:'pressure-opening'},{face:1,name:'出口',type:'pressure-opening'},
  {face:2,name:'wall',type:'wall'},{face:3,name:'symmetry',type:'symmetry'}],
 namedBoundaryFluxes:[{name:'入口',faces:1,length:2,inflow:.5,outflow:0,net:-.5,normalMeanVelocity:-.25},
  {name:'出口',faces:1,length:1,inflow:0,outflow:.5,net:.5,normalMeanVelocity:.5},
  {name:'wall',faces:1,length:4,inflow:0,outflow:0,net:0,normalMeanVelocity:0},
  {name:'symmetry',faces:1,length:4,inflow:0,outflow:0,net:0,normalMeanVelocity:0}]});
test('named flux uses signed outward volume per depth and preserves local inflow/outflow',()=>{
 validateBoundaryFluxes(sample());validateBoundaryFluxes({case:'custom'});
 const mixed={case:'custom',globalImbalance:0,boundaryFluxDefinition:DEFINITION,
  boundaryConditions:[{face:0,name:'opening',type:'pressure-opening'},{face:1,name:'opening',type:'pressure-opening'}],
  namedBoundaryFluxes:[{name:'opening',faces:2,length:2,inflow:.2,outflow:.2,net:0,normalMeanVelocity:0}]};
 validateBoundaryFluxes(mixed);assert.equal(mixed.namedBoundaryFluxes[0].inflow,.2);
});
test('named flux rejects omitted patches, incorrect units, signs, mean velocity and balances',()=>{
 for(const mutate of [s=>s.namedBoundaryFluxes.pop(),s=>s.namedBoundaryFluxes.push({...s.namedBoundaryFluxes[0]}),
  s=>s.boundaryFluxDefinition='kg/s',s=>delete s.namedBoundaryFluxes,s=>s.namedBoundaryFluxes[0].faces++,
  s=>s.namedBoundaryFluxes[0].inflow=-.5,s=>s.namedBoundaryFluxes[0].net=.5,
  s=>s.namedBoundaryFluxes[0].normalMeanVelocity=1,s=>s.globalImbalance=.1,
  s=>s.namedBoundaryFluxes[0].length=NaN,s=>s.namedBoundaryFluxes[0].length=0,
  s=>{Object.assign(s.namedBoundaryFluxes[3],{inflow:.1,outflow:.1});}]){
  const value=sample();mutate(value);assert.throws(()=>validateBoundaryFluxes(value),/命名边界流量/);
 }
});
