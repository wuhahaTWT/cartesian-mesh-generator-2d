'use strict';

(() => {
  const MAX_PIXELS = 40_000_000;
  const MAX_DIMENSION = 20_000;
  const PREVIEW_LIMIT = 1024;
  let active = null;

  const el = (tag, className, text) => {
    const node = document.createElement(tag);
    if (className) node.className = className;
    if (text !== undefined) node.textContent = text;
    return node;
  };

  const control = (tag, attrs = {}) => {
    const node = document.createElement(tag);
    Object.entries(attrs).forEach(([key, value]) => {
      if (key === 'text') node.textContent = value;
      else if (key === 'checked') node.checked = value;
      else node.setAttribute(key, value);
    });
    return node;
  };

  function buildDialog(name) {
    const dialog = el('dialog', 'raster-import');
    dialog.setAttribute('aria-labelledby', 'rasterImportTitle');
    const form = el('div', 'raster-import__form');
    const head = el('header', 'raster-import__head');
    const heading = el('div');
    const title = el('h1', '', '从图片提取二维轮廓');
    title.id = 'rasterImportTitle';
    heading.append(title, el('p', '', name || 'PNG / JPG'));
    const close = control('button', { type: 'button', class: 'raster-import__close', text: '取消' });
    head.append(heading, close);

    const body = el('div', 'raster-import__body');
    const preview = el('section', 'raster-import__preview');
    const stage = el('div', 'raster-import__stage');
    const canvas = el('canvas', 'raster-import__canvas');
    canvas.setAttribute('aria-label', '原图与提取轮廓预览；点击可选择前景区域');
    stage.append(canvas);
    const caption = el('p', 'raster-import__caption', '正在读取图片…');
    preview.append(stage, caption);

    const controls = el('aside', 'raster-import__controls');
    const extraction = el('section', 'raster-import__group');
    extraction.append(el('h2', '', '轮廓识别'));
    const modeLabel = el('label', 'raster-import__field', '前景判断');
    const mode = control('select', { 'data-raster': 'mode' });
    [['auto','自动判断'], ['dark','深色物体 / 线条'], ['light','浅色物体'],
      ['background','与背景不同的区域'], ['alpha','透明区域边缘']].forEach(([value, text]) => {
      const option = control('option', { value, text }); mode.append(option);
    });
    modeLabel.append(mode);

    const autoLabel = el('label', 'raster-import__check');
    const autoThreshold = control('input', { type: 'checkbox', checked: true, 'data-raster': 'auto-threshold' });
    autoLabel.append(autoThreshold, document.createTextNode(' 自动选择阈值'));
    const thresholdRow = el('div', 'raster-import__row');
    const threshold = control('input', { type: 'range', min: '0', max: '255', step: '1', value: '128', disabled: '', 'data-raster': 'threshold' });
    const thresholdValue = control('input', { type: 'number', min: '0', max: '255', step: '1', value: '128', disabled: '', 'aria-label': '阈值数值' });
    thresholdRow.append(threshold, thresholdValue);

    const selectionLabel = el('label', 'raster-import__field', '保留区域');
    const selection = control('select', { 'data-raster': 'selection' });
    [['largest','只保留最大区域'], ['all','保留全部有效区域']].forEach(([value, text]) => selection.append(control('option', { value, text })));
    selectionLabel.append(selection);

    const fillLabel = el('label', 'raster-import__check');
    const fillHoles = control('input', { type: 'checkbox', 'data-raster': 'fill-holes' });
    fillLabel.append(fillHoles, document.createTextNode(' 填充内部（闭合线稿，忽略内孔）'));
    const seedRow = el('div', 'raster-import__seed');
    const seedText = el('span', '', '可点击左侧图片选择一个前景区域');
    const clearSeed = control('button', { type: 'button', text: '清除点选', disabled: '' });
    seedRow.append(seedText, clearSeed);
    extraction.append(modeLabel, autoLabel, thresholdRow, selectionLabel, fillLabel, seedRow,
      el('p', 'raster-import__help', '适合正视、平面轮廓和背景干净的图片。复杂场景可能需要先裁图或清理背景。'));

    const calibrationGroup = el('section', 'raster-import__group');
    calibrationGroup.append(el('h2', '', '真实尺寸（必填）'));
    const calibrationLabel = el('label', 'raster-import__field', '所提取轮廓整体包围盒的实际宽度');
    const calibrationRow = el('div', 'raster-import__row');
    const calibration = control('input', { type: 'number', min: '0', step: 'any', placeholder: '例如 120', 'data-raster': 'calibration' });
    const unit = control('select', { 'data-raster': 'unit' });
    [['mm','mm'], ['cm','cm'], ['m','m']].forEach(([value, text]) => unit.append(control('option', { value, text })));
    calibrationRow.append(calibration, unit);
    calibrationLabel.append(calibrationRow);
    calibrationGroup.append(calibrationLabel,
      el('p', 'raster-import__help', '这里填写轮廓本身的宽度，不是整张照片的宽度。图片没有可靠的物理比例，因此不会自动猜尺寸。'));

    const limits = el('section', 'raster-import__group');
    limits.append(el('h2', '', '识别结果'));
    const status = el('p', 'raster-import__status', '等待图片解码。');
    status.setAttribute('role', 'status');
    status.setAttribute('aria-live', 'polite');
    limits.append(status,
      el('p', 'raster-import__help', '只生成二维闭合轮廓。照片无法自动恢复三维形状、被遮挡细节或真实尺寸。'));
    controls.append(extraction, calibrationGroup, limits);
    body.append(preview, controls);

    const foot = el('footer', 'raster-import__foot');
    const footnote = el('p', 'raster-import__help', '确认后将保存原图、识别设置和叠加预览。');
    const actions = el('div', 'raster-import__actions');
    const cancel = control('button', { type: 'button', text: '取消' });
    const confirm = control('button', { type: 'button', class: 'raster-import__confirm', text: '使用此轮廓', disabled: '', 'data-raster': 'confirm' });
    actions.append(cancel, confirm);
    foot.append(footnote, actions);
    form.append(head, body, foot);
    dialog.append(form);

    return { dialog, canvas, caption, mode, autoThreshold, threshold, thresholdValue,
      selection, fillHoles, seedText, clearSeed, calibration, unit, status, confirm,
      closeButtons: [close, cancel] };
  }

  function describeStats(result) {
    const stats = result.stats || {};
    const loopCount = Array.isArray(result.loops) ? result.loops.length : 0;
    const vertices = Array.isArray(result.loops)
      ? result.loops.reduce((sum, loop) => sum + loop.length, 0) : 0;
    const bits = [`找到 ${loopCount} 个闭合环，${vertices} 个顶点`];
    if (Number.isFinite(stats.components)) bits.push(`${stats.components} 个候选区域`);
    if (Number.isFinite(stats.threshold)) bits.push(`阈值 ${Math.round(stats.threshold)}`);
    return bits.join('；');
  }

  function open(filePath) {
    if (active) active.finish(null);
    return new Promise(resolve => {
      const ui = buildDialog(filePath && filePath.split(/[\\/]/).pop());
      document.body.append(ui.dialog);
      const session = {
        filePath, ui, resolve, worker: null, image: null, rgba: null, result: null,
        requestId: 0, timer: null, seed: null, busy: false, closed: false,
        source: null, sampledWidth: 0, sampledHeight: 0, originalWidth: 0, originalHeight: 0,
        pending: new Map(), readySettled: false, committing: false
      };
      session.ready = new Promise(readyResolve => { session.readyResolve = readyResolve; });
      active = session;

      session.finish = value => {
        if (session.closed) return;
        session.closed = true;
        clearTimeout(session.timer);
        if (session.worker) session.worker.terminate();
        session.pending.forEach(done => done(null));
        session.pending.clear();
        if (!session.readySettled) { session.readySettled = true; session.readyResolve(null); }
        if (session.image) { session.image.onload = null; session.image.onerror = null; session.image.src = ''; }
        if (ui.dialog.open) ui.dialog.close();
        ui.dialog.remove();
        if (active === session) active = null;
        resolve(value);
      };

      const setStatus = (message, kind = '') => {
        ui.status.textContent = message;
        ui.status.dataset.kind = kind;
      };

      const calibrationValue = () => Number(ui.calibration.value);
      const canConfirm = () => !session.busy && session.result && session.result.loops &&
        session.result.loops.length > 0 && Number.isFinite(calibrationValue()) && calibrationValue() > 0;
      const updateConfirm = () => { ui.confirm.disabled = !canConfirm(); };

      const draw = (showSeed = true) => {
        if (!session.image || !session.sampledWidth || !session.sampledHeight) return;
        const ctx = ui.canvas.getContext('2d');
        ctx.clearRect(0, 0, ui.canvas.width, ui.canvas.height);
        ctx.drawImage(session.image, 0, 0, ui.canvas.width, ui.canvas.height);
        if (!session.result || !Array.isArray(session.result.loops)) return;
        ctx.save();
        ctx.strokeStyle = '#f03b46';
        ctx.lineWidth = Math.max(1.5, Math.min(ui.canvas.width, ui.canvas.height) / 420);
        ctx.lineJoin = 'round';
        ctx.lineCap = 'round';
        ctx.shadowColor = 'rgba(255,255,255,.85)';
        ctx.shadowBlur = ctx.lineWidth * 1.1;
        for (const loop of session.result.loops) {
          if (!loop.length) continue;
          ctx.beginPath();
          ctx.moveTo(loop[0][0], loop[0][1]);
          for (let i = 1; i < loop.length; i += 1) ctx.lineTo(loop[i][0], loop[i][1]);
          ctx.closePath();
          ctx.stroke();
        }
        if (showSeed && session.seed) {
          ctx.shadowBlur = 0;
          ctx.fillStyle = '#f03b46';
          ctx.beginPath(); ctx.arc(session.seed.x, session.seed.y, Math.max(3, ctx.lineWidth * 2), 0, Math.PI * 2); ctx.fill();
        }
        ctx.restore();
      };

      const options = () => {
        const value = {
          mode: ui.mode.value,
          selection: ui.selection.value,
          fillHoles: ui.fillHoles.checked
        };
        if (!ui.autoThreshold.checked) value.threshold = Number(ui.threshold.value);
        if (session.seed) value.seed = { x: session.seed.x, y: session.seed.y };
        return value;
      };

      const process = () => {
        if (session.closed || !session.rgba) return Promise.resolve(null);
        clearTimeout(session.timer);
        session.result = null;
        session.busy = true;
        ui.dialog.setAttribute('aria-busy', 'true');
        updateConfirm();
        setStatus('正在提取轮廓…');
        const id = ++session.requestId;
        if (!session.worker) {
          session.worker = new Worker('raster-worker.js');
          session.worker.onmessage = event => {
            const done = session.pending.get(event.data.id);
            session.pending.delete(event.data.id);
            if (session.closed) { if (done) done(null); return; }
            if (event.data.id !== session.requestId) { if (done) done(null); return; }
            session.busy = false;
            ui.dialog.setAttribute('aria-busy', 'false');
            if (event.data.error) {
              setStatus(event.data.error, 'error');
              draw(); updateConfirm();
              if (!session.readySettled) { session.readySettled = true; session.readyResolve(null); }
              if (done) done(null);
              return;
            }
            session.result = event.data.result;
            const warnings = Array.isArray(session.result.warnings) ? session.result.warnings : [];
            const message = [describeStats(session.result), ...warnings].join('\n');
            setStatus(message, warnings.length ? 'warning' : '');
            const automaticThreshold = Number.isFinite(session.result.threshold)
              ? session.result.threshold : session.result.stats && session.result.stats.threshold;
            if (ui.autoThreshold.checked && Number.isFinite(automaticThreshold)) {
              const detected = Math.round(automaticThreshold);
              ui.threshold.value = String(detected);
              ui.thresholdValue.value = String(detected);
            }
            draw(); updateConfirm();
            if (!session.readySettled) { session.readySettled = true; session.readyResolve(session.result); }
            if (done) done(session.result);
          };
          session.worker.onerror = event => {
            if (session.closed) return;
            session.busy = false;
            ui.dialog.setAttribute('aria-busy', 'false');
            setStatus(event.message || '轮廓处理线程失败。', 'error');
            updateConfirm();
            const failedId = session.requestId;
            const done = session.pending.get(failedId);
            session.pending.delete(failedId);
            if (!session.readySettled) { session.readySettled = true; session.readyResolve(null); }
            if (done) done(null);
          };
        }
        return new Promise(done => {
          session.pending.set(id, done);
          const copy = new Uint8ClampedArray(session.rgba);
          session.worker.postMessage({ id, image: {
            width: session.sampledWidth, height: session.sampledHeight, data: copy
          }, options: options() }, [copy.buffer]);
        });
      };

      const schedule = () => {
        clearTimeout(session.timer);
        // Invalidate the visible result immediately.  During the debounce window
        // the controls already describe new settings, so the old contour must not
        // remain confirmable or be paired with those settings at commit time.
        session.requestId += 1;
        session.result = null;
        session.busy = true;
        ui.dialog.setAttribute('aria-busy', 'true');
        setStatus('正在提取轮廓…');
        draw();
        updateConfirm();
        session.timer = setTimeout(process, 100);
      };

      const setControlBusy = busy => {
        session.busy = busy;
        session.committing = busy;
        ui.dialog.setAttribute('aria-busy', String(busy));
        ui.closeButtons.forEach(button => { button.disabled = busy; });
        ui.mode.disabled = busy;
        ui.autoThreshold.disabled = busy;
        ui.threshold.disabled = busy || ui.autoThreshold.checked;
        ui.thresholdValue.disabled = busy || ui.autoThreshold.checked;
        ui.selection.disabled = busy;
        ui.fillHoles.disabled = busy;
        ui.clearSeed.disabled = busy || !session.seed;
        ui.calibration.disabled = busy;
        ui.unit.disabled = busy;
        updateConfirm();
      };

      ui.closeButtons.forEach(button => button.addEventListener('click', () => session.finish(null)));
      ui.dialog.addEventListener('cancel', event => { event.preventDefault(); if (!session.committing) session.finish(null); });
      [ui.mode, ui.selection, ui.fillHoles].forEach(node => node.addEventListener('change', schedule));
      ui.autoThreshold.addEventListener('change', () => {
        ui.threshold.disabled = ui.autoThreshold.checked;
        ui.thresholdValue.disabled = ui.autoThreshold.checked;
        schedule();
      });
      ui.threshold.addEventListener('input', () => {
        ui.thresholdValue.value = ui.threshold.value; schedule();
      });
      ui.thresholdValue.addEventListener('input', () => {
        const number = Math.max(0, Math.min(255, Number(ui.thresholdValue.value)));
        if (Number.isFinite(number)) ui.threshold.value = String(number);
        schedule();
      });
      ui.calibration.addEventListener('input', updateConfirm);
      ui.unit.addEventListener('change', updateConfirm);
      ui.clearSeed.addEventListener('click', () => {
        session.seed = null; ui.seedText.textContent = '可点击左侧图片选择一个前景区域';
        ui.clearSeed.disabled = true; schedule();
      });
      ui.canvas.addEventListener('click', event => {
        if (!session.rgba || session.busy) return;
        const rect = ui.canvas.getBoundingClientRect();
        if (!(rect.width > 0 && rect.height > 0)) return;
        session.seed = {
          x: Math.max(0, Math.min(session.sampledWidth - 1, (event.clientX - rect.left) * session.sampledWidth / rect.width)),
          y: Math.max(0, Math.min(session.sampledHeight - 1, (event.clientY - rect.top) * session.sampledHeight / rect.height))
        };
        ui.seedText.textContent = `已点选 (${Math.round(session.seed.x)}, ${Math.round(session.seed.y)})`;
        ui.clearSeed.disabled = false; schedule();
      });
      const commit = async () => {
        if (!canConfirm()) return;
        setControlBusy(true);
        ui.confirm.textContent = '正在保存…';
        setStatus('正在保存原图、轮廓和识别记录…');
        try {
          draw(false);
          const overlayDataUrl = ui.canvas.toDataURL('image/png');
          draw(true);
          const saved = await window.cartmesh.commitRaster({
            sourcePath: session.filePath,
            loops: session.result.loops,
            pixelWidth: session.sampledWidth,
            pixelHeight: session.sampledHeight,
            calibration: { width: calibrationValue(), unit: ui.unit.value },
            settings: {
              ...options(),
              sourcePixelWidth: session.originalWidth,
              sourcePixelHeight: session.originalHeight,
              processingScale: session.sampledWidth / session.originalWidth
            },
            stats: session.result.stats || {},
            warnings: session.result.warnings || [],
            overlayDataUrl
          });
          session.finish(saved && { geometryPath: saved.geometryPath, label: saved.label });
          return saved;
        } catch (error) {
          setControlBusy(false);
          ui.confirm.textContent = '使用此轮廓';
          setStatus(error && error.message ? error.message : String(error), 'error');
          return null;
        }
      };
      ui.confirm.addEventListener('click', commit);

      ui.dialog.showModal();
      ui.closeButtons[0].focus();

      (async () => {
        try {
          if (!window.cartmesh || typeof window.cartmesh.readRaster !== 'function' || typeof window.cartmesh.commitRaster !== 'function') {
            throw new Error('当前桌面后端尚未提供图片导入接口。');
          }
          session.source = await window.cartmesh.readRaster(filePath);
          if (session.closed) return;
          const image = new Image();
          session.image = image;
          await new Promise((resolveImage, rejectImage) => {
            image.onload = resolveImage;
            image.onerror = () => rejectImage(new Error('无法解码这张图片，请确认文件未损坏。'));
            image.src = session.source.dataUrl;
          });
          if (session.closed) return;
          const width = image.naturalWidth;
          const height = image.naturalHeight;
          if (!(width > 0 && height > 0)) throw new Error('图片没有有效的像素尺寸。');
          if (width > MAX_DIMENSION || height > MAX_DIMENSION || width * height > MAX_PIXELS) {
            throw new Error(`图片为 ${width} × ${height}，超过 4000 万像素或单边 20000 像素限制。请先裁剪或缩小。`);
          }
          session.originalWidth = width; session.originalHeight = height;
          const scale = Math.min(1, PREVIEW_LIMIT / width, PREVIEW_LIMIT / height);
          session.sampledWidth = Math.max(1, Math.round(width * scale));
          session.sampledHeight = Math.max(1, Math.round(height * scale));
          ui.canvas.width = session.sampledWidth;
          ui.canvas.height = session.sampledHeight;
          const ctx = ui.canvas.getContext('2d', { willReadFrequently: true });
          ctx.drawImage(image, 0, 0, session.sampledWidth, session.sampledHeight);
          session.rgba = ctx.getImageData(0, 0, session.sampledWidth, session.sampledHeight).data;
          const ratio = width / session.sampledWidth;
          ui.caption.textContent = scale < 1
            ? `原图 ${width} × ${height} px；按 ${session.sampledWidth} × ${session.sampledHeight} px 识别（每个处理像素约对应原图 ${ratio.toFixed(2)} px）。缩小不会恢复或推断细节。`
            : `按原始分辨率 ${width} × ${height} px 识别；红线为待导入轮廓。`;
          await process();
        } catch (error) {
          if (session.closed) return;
          session.busy = false;
          setStatus(error && error.message ? error.message : String(error), 'error');
          ui.caption.textContent = '图片读取失败。';
          updateConfirm();
          if (!session.readySettled) { session.readySettled = true; session.readyResolve(null); }
        }
      })();

      session.process = process;
      session.confirm = commit;
    });
  }

  window.CartMeshRasterImport = { open };
  const smoke = {
    open,
    setCalibration(width, unit = 'mm') {
      if (!active) return;
      active.ui.calibration.value = String(width);
      active.ui.unit.value = unit;
      active.ui.calibration.dispatchEvent(new Event('input'));
    },
    async setOptions(values = {}) {
      if (!active) return null;
      if (values.mode !== undefined) active.ui.mode.value = values.mode;
      if (values.selection !== undefined) active.ui.selection.value = values.selection;
      if (values.fillHoles !== undefined) active.ui.fillHoles.checked = Boolean(values.fillHoles);
      if (values.threshold === undefined || values.threshold === null) {
        active.ui.autoThreshold.checked = true;
      } else {
        active.ui.autoThreshold.checked = false;
        active.ui.threshold.value = String(values.threshold);
        active.ui.thresholdValue.value = String(values.threshold);
      }
      active.ui.threshold.disabled = active.ui.autoThreshold.checked;
      active.ui.thresholdValue.disabled = active.ui.autoThreshold.checked;
      if (values.seed === null) active.seed = null;
      else if (values.seed) active.seed = { x: Number(values.seed.x), y: Number(values.seed.y) };
      return active.process();
    },
    process() { return active ? active.process() : Promise.resolve(null); },
    confirm() { return active ? active.confirm() : undefined; },
    cancel() { return active ? active.finish(null) : undefined; },
    getResult() { return active ? active.result : null; },
    state() { return active ? { busy: active.busy, seed: active.seed, result: active.result,
      status: active.ui.status.textContent, confirmDisabled: active.ui.confirm.disabled } : null; }
  };
  Object.defineProperty(smoke, 'ready', { get: () => active ? active.ready : Promise.resolve(null) });
  window.__rasterSmoke = smoke;
})();
