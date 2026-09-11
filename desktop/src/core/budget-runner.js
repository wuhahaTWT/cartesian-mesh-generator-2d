'use strict';

const { planBudget, refineBudget, assessBudget } = require('./cell-budget');

// Keep every attempt's provenance and select the nearest successfully exported
// solver mesh. A later failed attempt must not destroy an earlier valid result.
async function runBudget(request, frame, { generate, progress = () => {}, signal }) {
  let choice = planBudget(request, frame);
  const target = choice.budgetPlan.targetCells;
  const attempts = [];
  const originalFirstLayer = choice.firstLayerRelativeSize;
  let layerRetries = 0;
  let best = null;
  let lastError;
  const seen = new Set();
  for (let i = 0; i < 3; i++) {
    if (signal?.aborted) throw new Error('操作已取消');
    const key = JSON.stringify([choice.wallRelativeSize, choice.backgroundRelativeSize,
      choice.firstLayerRelativeSize, choice.cellsPerLevel]);
    if (seen.has(key)) break;
    seen.add(key);
    progress({ attempt: i + 1, maximum: 3, parameters: choice });
    const started = Date.now();
    try {
      const payload = await generate(choice);
      if (signal?.aborted) throw new Error('操作已取消');
      if (payload.incomplete) throw new Error(payload.incomplete);
      if (choice.method === 'hybrid' && payload.result.actualMethod !== 'hybrid') {
        const error = new Error('未保留所选混合网格方法，纯网格回退结果不作为成功结果。');
        error.layerFallback = true; throw error;
      }
      lastError = null;
      const actual = payload.mesh.cells.length;
      const assessment = assessBudget(actual, target);
      attempts.push({ parameters: choice, seconds: (Date.now() - started) / 1000,
        success: true, actualCells: actual, assessment });
      const distance = Math.abs(actual - target) / target;
      if (!best || distance < best.distance) best = { payload, choice, assessment, distance };
      if (assessment.withinRange) break;
      choice = refineBudget(choice, actual, target);
      if (!choice) break;
    } catch (error) {
      lastError = error;
      attempts.push({ parameters: choice, seconds: (Date.now() - started) / 1000,
        success: false, reason: error.message });
      if (signal?.aborted) throw new Error('操作已取消');
      // A failed layer construction may be retried with a thinner/thicker first
      // layer. This bounded geometric adjustment never changes the quality gates,
      // number of layers, domain, or requested method.
      if (error.layerFallback && i < 2) {
        layerRetries++;
        choice = { ...choice, firstLayerRelativeSize: originalFirstLayer * (layerRetries === 1 ? 0.65 : 1.5),
          budgetPlan: { ...choice.budgetPlan, attempt: i+2, correctionKind: 'first-layer-retry' } };
        continue;
      }
      break;
    }
  }
  if (!best) throw lastError || new Error('未得到可导出的网格。');
  const payload = best.payload;
  payload.automatic = true;
  payload.selectedRequest = best.choice;
  payload.attempts = attempts;
  payload.cellBudget = { ...best.assessment, targetCells: target,
    actualCells: payload.mesh.cells.length, attempts: attempts.length,
    reached: best.assessment.withinRange, stoppedReason: lastError?.message || null };
  return payload;
}

module.exports = { runBudget };
