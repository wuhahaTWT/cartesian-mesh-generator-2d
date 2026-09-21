'use strict';
const DEFINITION='m2/s per unit depth; outward positive; inflow/outflow nonnegative; corrected face flux';
const VALUES=['length','inflow','outflow','net','normalMeanVelocity'];
function validateBoundaryFluxes(summary) {
  if(summary.namedBoundaryFluxes===undefined && summary.boundaryFluxDefinition===undefined)return; // Legacy output.
  const fail=()=>{throw new Error('命名边界流量与类型、单位或守恒汇总不一致。');};
  const finite=x=>typeof x==='number' && Number.isFinite(x);
  const near=(a,b,scale=0)=>finite(a)&&finite(b)&&Math.abs(a-b)<=1e-12+1e-10*Math.max(Math.abs(a),Math.abs(b),scale);
  if(summary.case!=='custom' || summary.boundaryFluxDefinition!==DEFINITION || !Array.isArray(summary.namedBoundaryFluxes)
    || !Array.isArray(summary.boundaryConditions))fail();
  const groups=new Map();
  for(const b of summary.boundaryConditions){
    if(!groups.has(b.name))groups.set(b.name,{faces:0,type:b.type});
    const g=groups.get(b.name);if(g.type!==b.type)fail();g.faces++;
  }
  const seen=new Set();let total=0,scale=0;
  for(const row of summary.namedBoundaryFluxes){
    const g=groups.get(row?.name);
    if(!g || seen.has(row.name) || row.faces!==g.faces || Object.keys(row).length!==VALUES.length+2
      || !VALUES.every(k=>finite(row[k])) || row.length<=0 || row.inflow<0 || row.outflow<0)fail();
    if(!near(row.net,row.outflow-row.inflow,row.inflow+row.outflow) || !near(row.normalMeanVelocity,row.net/row.length))fail();
    if(['wall','moving-wall','smooth-moving-wall','symmetry'].includes(g.type) && (row.inflow!==0 || row.outflow!==0))fail();
    seen.add(row.name);total+=row.net;scale+=row.inflow+row.outflow;
  }
  if(seen.size!==groups.size || !near(total,summary.globalImbalance,scale))fail();
}
module.exports={validateBoundaryFluxes,DEFINITION};
