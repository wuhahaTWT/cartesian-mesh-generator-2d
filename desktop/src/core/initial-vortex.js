'use strict';

function normalizeInitialVortex(value) {
  if (!value || typeof value!=='object' || Array.isArray(value) ||
      !Array.isArray(value.centre) || value.centre.length!==2)
    throw new Error('初始局部涡需要二维中心、半径和带符号峰值速度。');
  const values=[...value.centre,value.radius,value.peakSpeed];
  if (values.some(v=>typeof v!=='number'||!Number.isFinite(v)) || !(value.radius>0))
    throw new Error('初始局部涡的数值必须有限，半径必须大于 0。');
  return {centre:[...value.centre],radius:value.radius,peakSpeed:value.peakSpeed};
}
function validateInitialVortexOutput(summary,request=null,startTime=0) {
  const value=summary.initialVortex;
  if(value!==undefined){
    if (!summary.temporalDiscretization || startTime!==0 ||
        value?.definition!=='compact-cubic-v1' || value?.checkpointSuffix!=='.initial.checkpoint')
      throw new Error('初始局部涡摘要只适用于从零时刻开始的非定常计算。');
    normalizeInitialVortex(value);
  }
  if (request) {
    const expected=request.initialVortex;
    if(Boolean(value)!==Boolean(expected) || (value && JSON.stringify(normalizeInitialVortex(value))!==JSON.stringify(expected)))
      throw new Error('初始局部涡结果与请求不一致。');
  }
}
module.exports={normalizeInitialVortex,validateInitialVortexOutput};
