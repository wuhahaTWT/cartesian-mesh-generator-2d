'use strict';

// `--smoke=<sample-id> [--out=<dir>] [--method=<id>] [--flow-convection=<scheme>] [--flow-pressure-preconditioner=<id>] [--shot=<png>]` drives the real
// renderer through a full run and writes a screenshot.  Clicking the actual controls
// is the only check that covers the renderer, the IPC surface and the CLI together.
async function runSmoke({ app, fs, mainWindow, getCurrentResult, getSessionDirectory, exportPackage }) {
  const argument = name => {
    const found = process.argv.find(item => item.startsWith(`--${name}=`));
    return found ? found.slice(name.length + 3) : null;
  };
  const sampleId = argument('smoke');
  const outputDirectory = argument('out') || '';
  const method = argument('method') || 'cutcell';
  const shot = argument('shot');
  // Optional physical sizing for repeatable large-mesh smoke runs. These
  // populate the real form after loading the geometry, through its events.
  const sizingInputs = Object.fromEntries([
    ['backgroundLevel','background-level'],
    ['backgroundMinimumLevel','background-minimum-level'],
    ['backgroundPadding','background-padding'],
    ['backgroundMode','background-mode'],
    ['wallRelativeSize', 'wall-relative-size'],
    ['backgroundRelativeSize', 'background-relative-size'],
    ['referenceLength', 'reference-length'],
    ['relativePadding', 'padding-relative-size'],
    ['relativeBandCells', 'band-cells'],
    ['firstLayerRelativeSize', 'first-layer-relative-size'],
    ['nLayers', 'layers'],
    ['growthRatio', 'growth-ratio'],
    ['extrusionRelativeSize', 'extrusion-relative-size'],
    ['smallAlpha', 'small-alpha'],
    ['autoPadding', 'auto-padding'],
    ['customTargetCells', 'target-cells'],
  ].map(([id, flag]) => [id, argument(flag)]).filter(([, value]) => value !== null));
  if (outputDirectory) await fs.mkdir(outputDirectory, { recursive: true });
  if (argument('small-window')) mainWindow.setSize(800, 560);

  // The renderer's init awaits the catalog over IPC, so the hook appears a moment
  // after the page finishes loading.
  for (let attempt = 0; attempt < 100; attempt++) {
    if (await mainWindow.webContents.executeJavaScript('Boolean(window.__smoke)')) break;
    await new Promise(resolve => setTimeout(resolve, 100));
  }

  await mainWindow.webContents.executeJavaScript(`(async () => {
    const loadedFonts = await document.fonts.load('13px "CartMesh UI"', '生成网格');
    await document.fonts.ready;
    if (!loadedFonts.length || loadedFonts.some(font => font.status !== 'loaded'))
      throw new Error('Bundled Chinese interface font did not load');
    const smoke = window.__smoke;
    smoke.setSidebarCollapsed(false);
    smoke.setOutput(${JSON.stringify(outputDirectory)});
    smoke.selectMethod(${JSON.stringify(method)});
    document.getElementById('controlMode').value = ${JSON.stringify(argument('control') || 'auto')};
    document.getElementById('controlMode').dispatchEvent(new Event('change'));
    document.getElementById('density').value = ${JSON.stringify(argument('target-cells') ? 'custom' : argument('density') || '5000')};
    if (${JSON.stringify(Boolean(argument('regions')))}) {
      smoke.addRegion();
      smoke.addRegion();
      smoke.state.regions[1].xmin = 6; smoke.state.regions[1].xmax = 14;
      smoke.state.regions[1].ymin = -1.6; smoke.state.regions[1].ymax = 1.6;
      smoke.state.regions[1].levelsBelowWall = 6;
      smoke.renderRegions();
    }
    const sample = smoke.state.catalog.samples.find(item => item.id === ${JSON.stringify(sampleId)});
    if (!sample) throw new Error('unknown sample ' + ${JSON.stringify(sampleId)});
    document.getElementById('sample').value = sample.id;
    await smoke.chooseGeometry(sample.path, sample.label, sample);
    const sampleField = sample.fluidRegion === 'interior' && sample.interiorSizeField || sample.sizeField;
    if (sampleField?.farFieldSpans !== undefined &&
        Number(document.getElementById('autoPadding').value) !== sampleField.farFieldSpans)
      throw new Error('Loading sample did not apply its physical domain to automatic mode');
    if (${JSON.stringify(Boolean(argument('image')))}) {
      const previousPath = smoke.state.geometryPath;
      const pendingImport = smoke.importGeometryFile(${JSON.stringify(argument('image'))});
      const firstResult = await window.__rasterSmoke.ready;
      if (${JSON.stringify(Boolean(argument('image-reject')))}) {
        if (firstResult || !window.__rasterSmoke.state().confirmDisabled) throw new Error('Invalid image was accepted');
        const reason = window.__rasterSmoke.state().status;
        window.__rasterSmoke.cancel(); await pendingImport;
        if (smoke.state.geometryPath !== previousPath) throw new Error('Cancel discarded previous geometry');
        return { rasterPreview: true, rejected: true, reason, previousGeometryPreserved: true };
      }
      if (!firstResult) throw new Error(window.__rasterSmoke.state().status);
      if (!window.__rasterSmoke.state().confirmDisabled) throw new Error('Image accepted without physical calibration');
      const options = ${JSON.stringify({ mode: argument('image-mode') || 'auto', fillHoles: argument('image-fill') === 'true' })};
      await window.__rasterSmoke.setOptions(options);
      const extracted = window.__rasterSmoke.getResult();
      if (!extracted) throw new Error(window.__rasterSmoke.state().status);
      smoke.state.rasterEvidence = { loopCount: extracted.loops.length,
        vertices: extracted.loops.reduce((n,loop)=>n+loop.length,0), stats: extracted.stats, warnings: extracted.warnings };
      if (${JSON.stringify(Boolean(argument('image-preview-only')))})
        return { rasterPreview: true, ...smoke.state.rasterEvidence, calibrationRequired: window.__rasterSmoke.state().confirmDisabled };
      window.__rasterSmoke.setCalibration(Number(${JSON.stringify(argument('image-width') || '200')}), 'mm');
      const committed = await window.__rasterSmoke.confirm();
      if (!committed) throw new Error(window.__rasterSmoke.state()?.status || 'Raster commit failed');
      await pendingImport;
      const expectedWidth = Number(${JSON.stringify(argument('image-width') || '200')}) / 1000;
      if (Math.abs(smoke.state.frame.width - expectedWidth) > expectedWidth * 1e-10 ||
          !document.getElementById('sourceUnits').disabled) throw new Error('Image calibration was not preserved');
    }
    if (${JSON.stringify(argument('verified-preset') === 'true')}) await smoke.loadVerifiedPreset();
    const sizingInputs = ${JSON.stringify(sizingInputs)};
    if ('referenceLength' in sizingInputs) {
      const referenceMode = document.getElementById('referenceMode');
      referenceMode.value = 'explicit';
      referenceMode.dispatchEvent(new Event('change'));
    }
    for (const [id, value] of Object.entries(sizingInputs)) {
      const input = document.getElementById(id);
      input.value = value;
      if (!input.checkValidity() || (input.type==='number' && !Number.isFinite(Number(input.value))))
        throw new Error('Invalid smoke sizing input: ' + id);
      input.dispatchEvent(new Event('input'));
      input.dispatchEvent(new Event('change'));
    }
    if (${JSON.stringify(argument('allow-unsafe') === 'true')}) {
      const input = document.getElementById('relativeAllowUnsafe');
      input.checked = true;
      input.dispatchEvent(new Event('change'));
    }
    document.getElementById('density').dispatchEvent(new Event('change'));
    const theme = document.getElementById('appTheme');
    theme.value = ${JSON.stringify(argument('theme') || 'modern')};
    theme.dispatchEvent(new Event('change'));
    const regionInput = document.querySelector('#regionList input');
    if (regionInput) {
      regionInput.focus();
      if (!regionInput.isConnected || document.activeElement !== regionInput) throw new Error('Region input lost focus');
    }
    if (${JSON.stringify(Boolean(argument('welcome-shot')))}) return { welcomeOnly: true };
    if (${JSON.stringify(Boolean(argument('cutcell-strategy')))}) {
      const strat = document.getElementById('cutcellStrategy');
      strat.value = ${JSON.stringify(argument('cutcell-strategy') || 'legacy')};
      strat.dispatchEvent(new Event('change'));
    }
    const pending = smoke.generate();
    if (!document.getElementById('sample').disabled || document.getElementById('cancel').hidden)
      throw new Error('Parameters are not locked during generation');
    await pending;
    if (${JSON.stringify(Boolean(argument('outlet-backflow')))}) {
      document.getElementById('flowOutletBackflow').value = ${JSON.stringify(argument('outlet-backflow') || 'reject')};
      document.getElementById('flowOutletBackflow').dispatchEvent(new Event('change'));
    }
    if (${JSON.stringify(Boolean(argument('flow')))}) {
      document.getElementById('flowCase').value = ${JSON.stringify(argument('flow') || 'external')};
      document.getElementById('flowCase').dispatchEvent(new Event('change'));
      document.getElementById('flowMaxIterations').value = ${JSON.stringify(argument('flow-max-iterations') || '20')};
      document.getElementById('flowTolerance').value = ${JSON.stringify(argument('flow-tolerance') || '0.000001')};
      document.getElementById('flowNu').value = ${JSON.stringify(argument('flow-nu') || '0.01')};
      document.getElementById('flowSpeed').value = ${JSON.stringify(argument('flow-speed') || '1')};
      document.getElementById('flowConvection').value = ${JSON.stringify(argument('flow-convection') || 'upwind')};
      document.getElementById('flowConvection').dispatchEvent(new Event('change'));
      document.getElementById('flowPressurePreconditioner').value = ${JSON.stringify(argument('flow-pressure-preconditioner') || 'ic0')};
      document.getElementById('flowPressurePreconditioner').dispatchEvent(new Event('change'));
      document.getElementById('flowSteadyAcceleration').value = ${JSON.stringify(argument('flow-steady-acceleration') || 'none')};
      document.getElementById('flowLinearPolicy').value = ${JSON.stringify(argument('flow-linear-policy') || 'strict')};
      document.getElementById('flowVelocityRelaxation').value = ${JSON.stringify(argument('flow-velocity-relaxation') || '.6')};
      document.getElementById('flowPressureCorrections').value = ${JSON.stringify(argument('flow-pressure-corrections') || '4')};
      document.getElementById('flowCoupling').value = ${JSON.stringify(argument('flow-coupling') || 'simple')};
      document.getElementById('flowPressureRelaxation').value = ${JSON.stringify(argument('flow-pressure-relaxation') || '.25')};
      document.getElementById('flowMomentumCorrectors').value = ${JSON.stringify(argument('flow-momentum-correctors') || '2')};
      document.getElementById('flowCoupling').dispatchEvent(new Event('change'));
      document.getElementById('flowSteadyAcceleration').dispatchEvent(new Event('change'));
      if (${JSON.stringify(argument('flow') === 'custom')}) {
        await smoke.prepareFlowBoundaries(${JSON.stringify(argument('flow-boundary-preset') || 'duct')});
        if (!smoke.state.flowBoundaryDefinition) throw new Error('Custom boundary editor did not receive definition');
        const pressure=document.querySelector('[data-boundary-key="p"][data-patch-name="outlet"]') || document.querySelector('[data-boundary-key="p"]');
        if (${JSON.stringify(argument('flow-outlet-pressure') !== null)}) {
          if (!pressure) throw new Error('No explicit pressure patch in editor');
          pressure.value=${JSON.stringify(argument('flow-outlet-pressure') || '0')};
          pressure.dispatchEvent(new Event('input'));
        }
        if (${JSON.stringify(argument('flow-opening-pressure') !== null)}) {
          const inlet=document.querySelector('[data-boundary-key="p"][data-patch-name="inlet"]');
          if(!inlet)throw new Error('No pressure-driven inlet in editor');
          inlet.value=${JSON.stringify(argument('flow-opening-pressure') || '0')};
          inlet.dispatchEvent(new Event('input'));
        }
        document.querySelector('#flowBoundaryPatches button').click();
        if (!smoke.view.boundaryHighlight.length) throw new Error('Patch location was not highlighted');
      }
      if (${JSON.stringify(Boolean(argument('flow-dt')))}) {
        document.getElementById('flowMode').value=${JSON.stringify(argument('flow-end-time') ? 'adaptive' : 'transient')};
        document.getElementById('flowMode').dispatchEvent(new Event('change'));
        document.getElementById('flowDt').value=${JSON.stringify(argument('flow-dt') || '.01')};
        document.getElementById('flowSteps').value=${JSON.stringify(argument('flow-steps') || '2')};
        document.getElementById('flowEndTime').value=${JSON.stringify(argument('flow-end-time') || '1')};
        document.getElementById('flowMaxCourant').value=${JSON.stringify(argument('flow-max-courant') || '1')};
      }
      if (${JSON.stringify(argument('flow-vortex-speed')!==null)}) {
        document.getElementById('flowInitialSettings').open=true;
        document.getElementById('flowInitialVortex').checked=true;
        document.getElementById('flowInitialVortex').dispatchEvent(new Event('change'));
        for(const [id,value] of Object.entries({flowVortexX:${JSON.stringify(argument('flow-vortex-x') || '3')},
          flowVortexY:${JSON.stringify(argument('flow-vortex-y') || '0')},flowVortexRadius:${JSON.stringify(argument('flow-vortex-radius') || '1')},
          flowVortexSpeed:${JSON.stringify(argument('flow-vortex-speed') || '.01')}}))document.getElementById(id).value=value;
        if(document.getElementById('flowVortexRadius').disabled)throw new Error('Initial vortex controls are disabled for fresh transient run');
      }
      if (${JSON.stringify(argument('flow-case-check')==='true')}) {
        const saved=await smoke.saveFlowCase();
        if(!saved)throw new Error('App did not save flow case');
        document.getElementById('flowNu').value='4';
        document.getElementById('flowTolerance').value='0.000001';
        document.getElementById('flowSpeed').value='2';
        document.getElementById('flowMode').value='steady';
        document.getElementById('flowMode').dispatchEvent(new Event('change'));
        document.getElementById('flowConvection').value='upwind';
        document.getElementById('flowSteadyAcceleration').value='none';
        document.getElementById('flowLinearPolicy').value='strict';
        document.getElementById('flowVelocityRelaxation').value='.4';
        document.getElementById('flowPressureCorrections').value='2';
        document.getElementById('flowCoupling').value='simple';
        document.getElementById('flowPressureRelaxation').value='.25';
        document.getElementById('flowMomentumCorrectors').value='2';
        document.getElementById('flowInitialVortex').checked=false;
        document.getElementById('flowVortexSpeed').value='7';
        smoke.state.flowBoundaryDefinition=null;
        const loaded=await smoke.loadFlowCase();
        const repeated=await smoke.saveFlowCase();
        if(!loaded || !repeated || JSON.stringify(saved.document)!==JSON.stringify(repeated.document))
          throw new Error('Saved flow settings did not survive editor modification and reload');
        if(document.getElementById('flowResume').checked || document.getElementById('thermalResume').checked)
          throw new Error('Loading a zero-time case silently enabled resume');
        smoke.state.flowCaseSmoke={mesh:saved.document.mesh,request:saved.document.request,settingsRoundTrip:true};
      }
      await smoke.runFlow();
      if (${JSON.stringify(argument('flow-case-check')==='true')}) {
        if(!smoke.state.flow)throw new Error('Restored case failed to calculate');
        if(document.getElementById('flowMode').value!=='steady' && !document.getElementById('saveFlowCase').disabled)
          throw new Error('Resume settings could be saved as a fresh case');
        const fields=JSON.stringify(smoke.state.flow.fields),history=JSON.stringify(smoke.state.flow.history);
        if(!await smoke.loadFlowCase())throw new Error('Case reload after calculation failed');
        if(smoke.state.flow || !document.getElementById('flowResult').hidden)
          throw new Error('Old result was incorrectly displayed as loaded case output');
        await smoke.runFlow();
        if(!smoke.state.flow || fields!==JSON.stringify(smoke.state.flow.fields) || history!==JSON.stringify(smoke.state.flow.history))
          throw new Error('Repeated saved case changed its actual field or history');
        smoke.state.flowCaseSmoke.repeatedCalculationIdentical=true;
      }
      if (${JSON.stringify(argument('flow-vortex-speed')!==null)}) {
        if(!smoke.state.flow?.summary.initialVortex || !smoke.state.flow.files.initialCheckpoint)
          throw new Error('Initial vortex summary/checkpoint missing from real App result');
        if(!document.getElementById('flowInitialVortex').disabled || !document.getElementById('flowVortexRadius').disabled)
          throw new Error('Resume would allow duplicate vortex initialization');
        smoke.state.initialVortexSmoke=smoke.state.flow.summary.initialVortex;
      }
      if (${JSON.stringify(Boolean(argument('flow-end-time')))}) {
        if (!smoke.state.flow || smoke.state.flow.summary.timeStepControl!=='adaptive-cfl-retry') throw new Error('Missing adaptive App result');
        const initial=smoke.state.flow.summary;
        const before=initial.acceptedTime;
        const target=before+Number(document.getElementById('flowDt').value);
        document.getElementById('flowEndTime').value=String(target);
        document.getElementById('flowResume').checked=true;
        document.getElementById('flowResume').dispatchEvent(new Event('change'));
        await smoke.runFlow();
        if(smoke.state.flow?.summary.acceptedTime!==target) throw new Error('Adaptive absolute-time resume failed');
        const completed=smoke.state.flow;
        const iterations=document.getElementById('flowMaxIterations').value;
        document.getElementById('flowMaxIterations').value='1';
        document.getElementById('flowMaxRetries').value='0';
        document.getElementById('flowEndTime').value=String(target+.04);
        await smoke.runFlow();
        if(smoke.state.flow?.summary.acceptedTime!==target || smoke.state.flowRestart?.time!==target)
          throw new Error('Adaptive failure replaced the previous accepted result');
        document.getElementById('flowMaxIterations').value=iterations;
        document.getElementById('flowMaxRetries').value='10';
        document.getElementById('flowEndTime').value=String(target+100);
        const pending=smoke.runFlow(),deadline=Date.now()+90000;
        while(smoke.state.busy && !smoke.state.flowHistory.length && Date.now()<deadline)
          await new Promise(resolve=>setTimeout(resolve,25));
        if(!smoke.state.busy || !smoke.state.flowHistory.length) throw new Error('Adaptive cancellation had no accepted live step');
        document.getElementById('cancel').click();await pending;
        const saved=smoke.state.flowRestart?.time;
        if(!(saved>target) || smoke.state.flow?.summary.acceptedTime!==completed.summary.acceptedTime)
          throw new Error('Adaptive cancellation lost accepted checkpoint or earlier result');
        document.getElementById('flowEndTime').value=String(saved+.02);
        await smoke.runFlow();
        if(smoke.state.flow?.summary.acceptedTime!==saved+.02) throw new Error('Adaptive cancelled restart used wrong physical time');
        smoke.state.adaptiveSmoke={initial,failedAcceptedTime:target,cancelledAcceptedTime:saved,
          resumedTime:smoke.state.flow.summary.acceptedTime};
      }
      if (${JSON.stringify(Boolean(argument('flow-resume-steps')))}) {
        const before=smoke.state.flow?.summary.acceptedTime;
        document.getElementById('flowSteps').value=${JSON.stringify(argument('flow-resume-steps') || '2')};
        document.getElementById('flowResume').checked=true;
        document.getElementById('flowResume').dispatchEvent(new Event('change'));
        if (${JSON.stringify(Boolean(argument('flow-resume-pressure-preconditioner')))}) {
          const method=document.getElementById('flowPressurePreconditioner');
          if (method.disabled) throw new Error('Pressure solver selector is locked during idle resume');
          method.value=${JSON.stringify(argument('flow-resume-pressure-preconditioner') || 'ic0')};
          method.dispatchEvent(new Event('change'));
        }
        await smoke.runFlow();
        const expected=before+Number(document.getElementById('flowDt').value)*Number(document.getElementById('flowSteps').value);
        if (!smoke.state.flow || Math.abs(smoke.state.flow.summary.acceptedTime-expected)>1e-10)
          throw new Error('Desktop resume did not advance accepted physical time');
      }
      if (${JSON.stringify(argument('flow-failure-check') === 'true')}) {
        const accepted=smoke.state.flow.summary.acceptedTime;
        const maxIterations=document.getElementById('flowMaxIterations').value;
        document.getElementById('flowMaxIterations').value='1';
        document.getElementById('flowSteps').value='2';
        await smoke.runFlow();
        if (smoke.state.flow?.summary.acceptedTime!==accepted || smoke.state.flowRestart?.time!==accepted)
          throw new Error('Rejected time step replaced the complete result or accepted checkpoint');
        document.getElementById('flowMaxIterations').value=maxIterations;
      }
      if (${JSON.stringify(argument('flow-cancel-check') === 'true')}) {
        const completed=smoke.state.flow;
        document.getElementById('flowSteps').value='10000';
        const pendingFlow=smoke.runFlow();
        const deadline=Date.now()+90000;
        while (smoke.state.busy && !smoke.state.flowHistory.length && Date.now()<deadline)
          await new Promise(resolve=>setTimeout(resolve,50));
        if (!smoke.state.busy || !smoke.state.flowHistory.length) throw new Error('Could not observe live accepted step for cancellation');
        document.getElementById('cancel').click();
        await pendingFlow;
        if (!smoke.state.flowRestart || smoke.state.flowRestart.time<=completed.summary.acceptedTime
            || smoke.state.flow?.summary.acceptedTime!==completed.summary.acceptedTime)
          throw new Error('Cancellation lost checkpoint or previous complete result');
        const resumeTime=smoke.state.flowRestart.time;
        document.getElementById('flowSteps').value='2';
        await smoke.runFlow();
        if (Math.abs(smoke.state.flow.summary.acceptedTime-resumeTime-2*Number(document.getElementById('flowDt').value))>1e-10)
          throw new Error('Resume from cancelled calculation did not use the saved time');
      }
      if (!smoke.state.flow || document.getElementById('flowSpeedOption').hidden ||
          document.getElementById('displayMode').value !== 'speed')
        throw new Error('Native flow result did not reach the renderer');
      if (${JSON.stringify(argument('adapt-flow') === 'true')}) {
        await smoke.suggestFlowRefinement();
        if (!smoke.state.regions.length)
          throw new Error('Adaptive refinement did not populate refinement regions');
        if (document.getElementById('adaptiveIndicatorOption').hidden ||
            document.getElementById('displayMode').value !== 'adaptive')
          throw new Error('Adaptive gradient heatmap did not activate in displayMode');
        smoke.state.adaptiveRefinementSmoke = {
          regionCount: smoke.state.regions.length,
          displayMode: document.getElementById('displayMode').value,
          indicatorRange: smoke.view.fieldRange,
          note: document.getElementById('flowAdaptNote').textContent
        };
      }
    }
    if (${JSON.stringify(argument('thermal') === 'true')}) {
      for(const [id,value] of Object.entries({flowCase:'external',flowNu:'.1',flowSpeed:'1',flowConvection:${JSON.stringify(argument('flow-convection') || 'limited-linear')},flowPressurePreconditioner:'aggregation',flowMaxIterations:'1500',flowDt:'.05',flowSteps:'2',thermalDiffusivity:'.1'})) document.getElementById(id).value=value;
      await smoke.runThermal();
      if(!smoke.state.thermal || smoke.state.thermal.summary.time!==.1)throw new Error('Thermal result did not reach renderer');
      if(smoke.state.thermal.summary.outletBackflow!==${JSON.stringify(argument('outlet-backflow') || 'reject')})
        throw new Error('Thermal outlet mode did not reach the native solver');
      const before=smoke.state.thermal;
      await smoke.runThermal();
      if(!smoke.state.thermal || Math.abs(smoke.state.thermal.summary.time-.2)>1e-12)throw new Error('Thermal resume did not advance');
      const completed=smoke.state.thermal;
      document.getElementById('flowMaxIterations').value='1';
      await smoke.runThermal();
      if(smoke.state.thermal.summary.time!==completed.summary.time)throw new Error('Failed thermal run replaced complete result');
      document.getElementById('flowMaxIterations').value='1500';
      document.getElementById('flowSteps').value='10000';
      const pending=smoke.runThermal();
      const deadline=Date.now()+30000;
      while(smoke.state.busy&&!smoke.state.thermalHistory.length&&Date.now()<deadline)await new Promise(r=>setTimeout(r,30));
      if(!smoke.state.busy||!smoke.state.thermalHistory.length)throw new Error('No live thermal step observed before cancel');
      await window.cartmesh.cancel();await pending;
      const restartTime=smoke.state.thermalRestart?.time;
      if(!(restartTime>completed.summary.time)||smoke.state.thermal.summary.time!==completed.summary.time)throw new Error('Thermal cancellation lost accepted state');
      document.getElementById('flowSteps').value='1';await smoke.runThermal();
      if(Math.abs(smoke.state.thermal.summary.time-restartTime-.05)>1e-12)throw new Error('Thermal cancel/resume clock differs');
      if(document.getElementById('displayMode').value!=='temperature'||document.getElementById('thermalOption').hidden)throw new Error('Temperature map not displayed');
      document.getElementById('thermalBlock').scrollIntoView({block:'start'});
    }
    if (${JSON.stringify(Boolean(argument('interaction-check')))}) {
      const mesh = smoke.state.mesh;
      const display = document.getElementById('displayMode').value;
      const sidebar = document.getElementById('toggleSidebar');
      const initiallyCollapsed = document.body.classList.contains('sidebar-collapsed');
      sidebar.click();
      if (document.body.classList.contains('sidebar-collapsed') === initiallyCollapsed ||
          sidebar.getAttribute('aria-expanded') !== String(!document.body.classList.contains('sidebar-collapsed')) ||
          localStorage.getItem('cartmesh2d-sidebar-collapsed') !== String(document.body.classList.contains('sidebar-collapsed'))) {
        throw new Error('Sidebar toggle did not update state or persistence');
      }
      sidebar.click();
      if (document.body.classList.contains('sidebar-collapsed') !== initiallyCollapsed)
        throw new Error('Sidebar toggle did not restore state');
      for (const value of ['duet','modern']) {
        window.CartMeshTheme.set(value);
        if (smoke.state.mesh !== mesh || document.getElementById('displayMode').value !== display ||
            localStorage.getItem(window.CartMeshTheme.storageKey) !== value)
          throw new Error('Theme switch changed mesh/display state or failed persistence');
      }
      window.CartMeshTheme.set(${JSON.stringify(argument('theme') || 'modern')});
      sidebar.focus();
      window.dispatchEvent(new KeyboardEvent('keydown', { key: 'b', metaKey: true }));
      if (document.body.classList.contains('sidebar-collapsed') === initiallyCollapsed)
        throw new Error('Sidebar Command+B shortcut did not toggle');
      window.dispatchEvent(new KeyboardEvent('keydown', { key: 'b', ctrlKey: true }));
      if (document.body.classList.contains('sidebar-collapsed') !== initiallyCollapsed)
        throw new Error('Sidebar Control+B shortcut did not restore state');
      const actual = smoke.state.selectedRequest;
      if (actual) {
      const previousMode = document.getElementById('controlMode').value;
      document.getElementById('actualToManual').click();
      if (document.getElementById('controlMode').value !== 'manual' ||
          Number(document.getElementById('wallRelativeSize').value) !== actual.wallRelativeSize ||
          Number(document.getElementById('relativeBandCells').value) !== actual.cellsPerLevel ||
          Number(document.getElementById('relativePadding').value) !== actual.farFieldSpans)
        throw new Error('Actual-to-manual transfer lost parameters');
      document.getElementById('controlMode').value = previousMode;
      document.getElementById('controlMode').dispatchEvent(new Event('change'));
      } else if (document.getElementById('controlMode').value !== 'manual') {
        throw new Error('Missing selected automatic parameters');
      }
    }
    if (${JSON.stringify(Boolean(argument('repeat')))}) await smoke.generate();
    const mode = ${JSON.stringify(argument('mode') || 'level')};
    if (mode !== 'level') {
      const select = document.getElementById('displayMode');
      select.value = mode;
      select.dispatchEvent(new Event('change'));
    }
    return {
      status: document.getElementById('statusTitle').textContent,
      detail: document.getElementById('statusText').textContent,
      counters: document.getElementById('counters').innerText,
      gates: document.getElementById('gates').innerText,
      histogram: document.getElementById('histogram').innerText,
      log: document.getElementById('log').textContent,
      cellBudget: smoke.state.cellBudget,
      raster: smoke.state.rasterEvidence,
      adaptive: smoke.state.adaptiveSmoke || null,
      adaptiveRefinement: smoke.state.adaptiveRefinementSmoke || null,
      initialVortex: smoke.state.initialVortexSmoke || null,
      flowCase: smoke.state.flowCaseSmoke || null,
      bundledChineseFontLoaded: true,
      theme: document.documentElement.dataset.theme,
      interactionChecks: ${JSON.stringify(Boolean(argument('interaction-check')))},
      thermal: smoke.state.thermal ? {summary:smoke.state.thermal.summary,request:smoke.state.thermal.request,files:smoke.state.thermal.files,
        fieldCells:smoke.state.thermal.fields.cells.length,historyRows:smoke.state.thermal.history.length,
        resultText:document.getElementById('thermalResult').innerText,restart:smoke.state.thermalRestart,
        displayMode:document.getElementById('displayMode').value,monitorVisible:!document.getElementById('thermalTimeline').hidden,
        failureAndCancelResumeChecked:true} : null,
      flow: smoke.state.flow ? {
        summary: smoke.state.flow.summary,
        fieldCells: smoke.state.flow.fields.cells.length,
        displayMode: document.getElementById('displayMode').value,
        resultText: document.getElementById('flowResult').innerText,
        historyRows: smoke.state.flow.history?.length || 0,
        monitorVisible: !document.getElementById('flowTimeline').hidden,
        restart: smoke.state.flowRestart,
        cancellationChecked: ${JSON.stringify(argument('flow-cancel-check') === 'true')},
        failureChecked: ${JSON.stringify(argument('flow-failure-check') === 'true')}
      } : null
    };
  })()`).then(async report => {
    if (report.rasterPreview) {
      await new Promise(resolve => setTimeout(resolve, 200));
      report.dialogLayout = await mainWindow.webContents.executeJavaScript(`(() => {
        const dialog = document.querySelector('.raster-import');
        if (!dialog) return null;
        const confirm = dialog.querySelector('[data-raster="confirm"]');
        const box = dialog.getBoundingClientRect();
        const button = confirm?.getBoundingClientRect();
        return { width: innerWidth, height: innerHeight, visible: box.top >= 0 && box.bottom <= innerHeight,
          confirmReachable: Boolean(button && button.top >= 0 && button.bottom <= innerHeight) };
      })()`);
      if (shot) {
        await fs.writeFile(shot, (await mainWindow.webContents.capturePage()).toPNG());
        await fs.writeFile(shot + '.json', JSON.stringify(report, null, 2));
      }
      if (argument('flow-require-converged')==='true' && !report.flow?.summary?.converged)
      throw new Error('App acceptance requires a converged flow; diagnostic output is retained.');
    console.log(JSON.stringify(report, null, 2));
      await mainWindow.webContents.executeJavaScript('window.__rasterSmoke.cancel()');
      await fs.rm(getSessionDirectory(), { recursive: true, force: true });
      app.exit(0); return;
    }
    if (report.welcomeOnly) {
      await mainWindow.webContents.executeJavaScript(`window.__smoke.state.mesh = null; document.getElementById('empty').hidden = false; window.__smoke.view.clear();`);
      await mainWindow.webContents.executeJavaScript(`(async () => {
        const artwork = new Image(); artwork.src = 'assets/duet-workshop.png'; await artwork.decode();
      })()`);
      // Wake the compositor before capturing a newly loaded CSS wallpaper.
      await mainWindow.webContents.capturePage(undefined, { stayAwake: true });
      await new Promise(resolve => setTimeout(resolve, 1000));
      await fs.writeFile(argument('welcome-shot'), (await mainWindow.webContents.capturePage(undefined, { stayAwake: true })).toPNG());
      await fs.rm(getSessionDirectory(), { recursive: true, force: true });
      app.exit(0); return;
    }
    if(getCurrentResult()?.background) {
      report.background = await mainWindow.webContents.executeJavaScript(`(async()=>{
        const rejected=[];
        for(const name of ['runFlow','runThermal']) {
          try {await window.cartmesh[name]({}); throw new Error('Background accepted by '+name);}
          catch(error){if(!error.message.includes('完整笛卡尔背景网格'))throw error;rejected.push(name);}
        }
        return {rejected,displayMode:document.getElementById('displayMode').value,
          panelsHidden:['flowBlock','thermalBlock'].every(id=>document.getElementById(id).hidden),
          legend:document.getElementById('legend').innerText};
      })()`);
      if(!report.background.panelsHidden)throw new Error('Background grid exposed flow controls');
    }
    if (argument('export')) report.exported = await exportPackage(argument('export'));
    if (argument('flow-dt')) {
      await mainWindow.webContents.executeJavaScript("document.getElementById('flowTimeline').scrollIntoView({block:'nearest'}); document.getElementById('flowMode').scrollIntoView({block:'start'});");
      await mainWindow.webContents.capturePage(undefined, { stayAwake: true });
      await new Promise(resolve => setTimeout(resolve, 300));
    }
    const captureFreshPng = async filePath => {
      await mainWindow.webContents.executeJavaScript('new Promise(r => requestAnimationFrame(() => requestAnimationFrame(r)))');
      await mainWindow.webContents.capturePage(undefined, { stayAwake: true });
      await new Promise(resolve => setTimeout(resolve, 180));
      await fs.writeFile(filePath, (await mainWindow.webContents.capturePage(undefined, { stayAwake: true })).toPNG());
    };
    if (argument('mesh-shot')) await captureFreshPng(argument('mesh-shot'));
    if (argument('adapt-remesh') === 'true') {
      await mainWindow.webContents.executeJavaScript(`(async () => {
        await window.__smoke.generate();
      })()`);
      const afterCells = await mainWindow.webContents.executeJavaScript('window.__smoke.state.result?.counts?.cells || 0');
      report.adaptiveRemeshCells = afterCells;
      if (argument('remesh-shot')) await captureFreshPng(argument('remesh-shot'));
      if (argument('remesh-export')) report.remeshExported = await exportPackage(argument('remesh-export'));
    }
    if (argument('home-check')) {
      report.previewTiming = await mainWindow.webContents.executeJavaScript(`(() => {
        const view = window.__smoke.view;
        const original = view.scale;
        const elapsed = [];
        for (let i = 0; i < 6; i++) {
          view.scale = original * (1 + i * .01);
          const start = performance.now(); view.draw();
          elapsed.push(performance.now() - start);
        }
        view.scale = original; view.draw();
        return { drawSubmissionMs: elapsed, note: 'Canvas command submission; excludes asynchronous GPU completion' };
      })()`);
      report.home = await mainWindow.webContents.executeJavaScript(`(async () => {
        await window.__smoke.returnToStart();
        const { state, view } = window.__smoke;
        return { emptyVisible: !document.getElementById('empty').hidden,
          released: !state.mesh && !state.result && !view.mesh && !view.meshCache && !view.outline,
          geometryCleared: !state.geometryPath,
          exportAvailable: !document.getElementById('exportResult').hidden };
      })()`);
      report.home.mainPreviewReleased = getCurrentResult()?.mesh === null;
      if (!Object.values(report.home).every(Boolean)) throw new Error('Return to start did not release preview or preserve export');
      if (argument('export')) report.home.exported = await exportPackage(argument('export').replace(/\.zip$/, '') + '-home.zip');
    }
    if(report.thermal && shot){
      mainWindow.setSize(1320,900);
      await new Promise(resolve=>setTimeout(resolve,200));
      await mainWindow.webContents.executeJavaScript("document.getElementById('thermalBlock').scrollIntoView({block:'start'}); document.querySelector('.results').scrollTop=document.querySelector('.results').scrollHeight;");
      await new Promise(resolve=>setTimeout(resolve,200));
      await fs.writeFile(shot.replace(/\.png$/, '-thermal.png'),(await mainWindow.webContents.capturePage()).toPNG());
    }
    mainWindow.setSize(800, 560);
    await new Promise(resolve => setTimeout(resolve, 200));
    report.layout = await mainWindow.webContents.executeJavaScript(`(() => {
      const panel = document.querySelector('.panel');
      panel.scrollTop = panel.scrollHeight;
      const button = document.getElementById('generate').getBoundingClientRect();
      const bodyRect = document.body.getBoundingClientRect();
      const shellRect = document.querySelector('.shell').getBoundingClientRect();
      const canvasRect = document.getElementById('canvas').getBoundingClientRect();
      const geometry = { bodyHeight: bodyRect.height, innerHeight,
        canvasHeight: canvasRect.height, shellBottom: shellRect.bottom };
      return { bottomReachable: button.bottom <= innerHeight && button.top >= 0,
        pageHeight: document.documentElement.scrollHeight, windowHeight: innerHeight,
        geometry,
        wallpaper: {
          image: getComputedStyle(document.body).backgroundImage,
          canvasBackground: getComputedStyle(document.getElementById('canvasWrap')).backgroundColor,
          welcomeHidden: document.getElementById('empty').hidden
        },
      previewCells: window.__smoke.state.mesh?.cells.length,
      exportedCells: window.__smoke.state.result?.openFoam.cells,
      sidebarCollapsed: document.body.classList.contains('sidebar-collapsed'),
      sidebarExpanded: document.getElementById('toggleSidebar')?.getAttribute('aria-expanded') === 'true' };
    })()`);
    if (argument('flow-dt')) {
      report.flow.smallWindowLayout = await mainWindow.webContents.executeJavaScript(`(() => {
        const a=document.getElementById('flowResult').getBoundingClientRect();
        const b=document.getElementById('flowTimeline').getBoundingClientRect();
        document.querySelector('.results').scrollTop=document.querySelector('.results').scrollHeight;
        return { resultBottom:a.bottom, monitorTop:b.top, separated:b.top>=a.bottom-1 };
      })()`);
      if (!report.flow.smallWindowLayout.separated) throw new Error('Transient monitor overlaps the result summary');
    }
    if (argument('flow-require-converged')==='true' && !report.flow?.summary?.converged)
      throw new Error('App acceptance requires a converged flow; diagnostic output is retained.');
    if (argument('flow-linear-policy') || argument('flow-velocity-relaxation') || argument('flow-pressure-corrections')) {
      report.flow.linearControls=await mainWindow.webContents.executeJavaScript(`(() => {
        const summary=window.__smoke.state.flow?.summary,result=document.getElementById('flowResult');
        const policy=document.getElementById('flowLinearPolicy').value;
        const relaxation=Number(document.getElementById('flowVelocityRelaxation').value);
        const pressureCorrections=Number(document.getElementById('flowPressureCorrections').value);
        if(!summary || summary.linearPolicy!==policy || summary.velocityRelaxation!==relaxation || summary.pressureCorrectionPasses!==pressureCorrections)
          throw new Error('Linear efficiency controls do not match actual result');
        if(policy==='adaptive' && (!summary.strictLinearFinal || !result.innerText.includes('严格复核通过')))
          throw new Error('Adaptive linear solve lacks strict final certification in App');
        if(!result.innerText.includes('速度松弛系数'))throw new Error('Relaxation not shown in actual App');
        if(!result.innerText.includes('压力校正次数'))throw new Error('Pressure corrections not shown in App');
        return {policy,relaxation,pressureCorrections,strictLinearFinal:summary.strictLinearFinal,strictAcceptedSteps:summary.strictAcceptedSteps};
      })()`);
    }
    if (argument('flow-coupling')) {
      report.flow.couplingControls=await mainWindow.webContents.executeJavaScript(`(() => {
        const summary=window.__smoke.state.flow?.summary;
        const coupling=document.getElementById('flowCoupling').value;
        const pressureRelaxation=Number(document.getElementById('flowPressureRelaxation').value);
        const momentumCorrectors=Number(document.getElementById('flowMomentumCorrectors').value);
        if(!summary || summary.coupling!==coupling || summary.pressureRelaxation!==pressureRelaxation ||
            (coupling==='pimple' && summary.momentumCorrectors!==momentumCorrectors))
          throw new Error('Coupling controls do not match actual native result');
        if(!document.getElementById('flowResult').innerText.includes(coupling.toUpperCase()))
          throw new Error('Coupling algorithm is absent from actual App result');
        return {coupling,pressureRelaxation,momentumCorrectors:summary.momentumCorrectors};
      })()`);
    }
    console.log(JSON.stringify(report, null, 2));
    if (!report.layout.bottomReachable) throw new Error('Sidebar bottom is inaccessible');
    if (Math.abs(report.layout.geometry.bodyHeight - report.layout.geometry.innerHeight) > 1 ||
        report.layout.geometry.canvasHeight <= 0 || report.layout.geometry.shellBottom > report.layout.geometry.innerHeight + 1)
      throw new Error('Window geometry is not constrained to the viewport');
    if (shot) {
      await mainWindow.webContents.executeJavaScript("document.querySelector('.panel').scrollTop=0");
      if (argument('flow-adaptive-shot')) {
        mainWindow.setSize(1320,900);
        await new Promise(resolve=>setTimeout(resolve,200));
        await mainWindow.webContents.executeJavaScript("document.getElementById('flowMode').scrollIntoView({block:'start'}); document.querySelector('.results').scrollTop=0;");
        await new Promise(resolve=>setTimeout(resolve,200));
      }
      if (argument('flow-initial-shot')) {
        await mainWindow.webContents.executeJavaScript("document.getElementById('flowInitialSettings').open=true; document.getElementById('flowInitialSettings').scrollIntoView({block:'start'}); document.querySelector('.results').scrollTop=0;");
      }
      if (argument('flow-case-check')) {
        await mainWindow.webContents.executeJavaScript("document.getElementById('saveFlowCase').scrollIntoView({block:'start'}); document.querySelector('.results').scrollTop=0;");
      }
      if (argument('flow-load-shot')) {
        mainWindow.setSize(1320,900);
        await new Promise(resolve=>setTimeout(resolve,200));
        await mainWindow.webContents.executeJavaScript(`(() => {
          const result=document.getElementById('flowResult');
          const loads=window.__smoke.state.flow?.summary.namedWallLoads;
          if (!loads?.length || !result.innerText.includes('力矩')) throw new Error('Wall loads did not reach the actual App');
          for (const load of loads) if (!result.innerText.includes(load.name)) throw new Error('Missing named wall load in App');
          result.scrollIntoView({block:'end'});
        })()`);
      }
      if (argument('flow-flux-shot')) {
        mainWindow.setSize(1320,900);
        await new Promise(resolve=>setTimeout(resolve,200));
        await mainWindow.webContents.executeJavaScript(`(() => {
          const result=document.getElementById('flowResult'), summary=window.__smoke.state.flow?.summary;
          if(!summary?.namedBoundaryFluxes?.length || !result.innerText.includes('净流量'))throw new Error('Boundary flux did not reach App');
          for(const flux of summary.namedBoundaryFluxes)if(!result.innerText.includes(flux.name+' · 净流量'))throw new Error('Missing named boundary flux in App');
          [...result.querySelectorAll('div')].find(row=>row.innerText.includes('边界体积流量'))?.scrollIntoView({block:'start'});
        })()`);
      }
      if (argument('flow-boundary-shot')) {
        mainWindow.setSize(1100,800);
        await new Promise(resolve=>setTimeout(resolve,200));
        await mainWindow.webContents.executeJavaScript("document.getElementById('flowBoundarySettings').scrollIntoView({block:'start'});");
      }
      if (argument('flow-acceleration-shot')) {
        mainWindow.setSize(1320,900);
        await new Promise(resolve=>setTimeout(resolve,200));
        await mainWindow.webContents.executeJavaScript(`(() => {
          const result=document.getElementById('flowResult'),summary=window.__smoke.state.flow?.summary;
          if(summary?.steadyAcceleration!=='anderson' || !(summary.accelerationAccepted>0) || !result.innerText.includes('历史迭代'))
            throw new Error('Steady acceleration did not reach the actual App result');
          if(document.getElementById('flowSteadyAcceleration').value!=='anderson')throw new Error('Acceleration control differs from result');
          document.getElementById('flowSteadyAcceleration').scrollIntoView({block:'center'});
          result.scrollIntoView({block:'start'});
        })()`);
      }
      if (argument('flow-linear-shot')) {
        mainWindow.setSize(1320,900);
        await new Promise(resolve=>setTimeout(resolve,200));
        await mainWindow.webContents.executeJavaScript("document.getElementById('flowLinearPolicy').scrollIntoView({block:'center'}); document.getElementById('flowResult').scrollIntoView({block:'start'});");
      }
      if (argument('flow-coupling-shot')) {
        mainWindow.setSize(1320,900);
        await mainWindow.webContents.executeJavaScript("document.getElementById('flowCoupling').scrollIntoView({block:'start'}); document.getElementById('flowResult').scrollIntoView({block:'start'});");
      }
      await new Promise(resolve => setTimeout(resolve, 400));
      await fs.writeFile(shot, (await mainWindow.webContents.capturePage()).toPNG());
      console.log(`screenshot=${shot}`);
      await fs.writeFile(shot + '.json', JSON.stringify(report, null, 2));
    }
    await fs.rm(getSessionDirectory(), { recursive: true, force: true });
    app.exit(/失败/.test(report.status) ? 1 : 0);
  }).catch(error => {
    console.error(error);
    app.exit(1);
  });
}

module.exports = { runSmoke };
