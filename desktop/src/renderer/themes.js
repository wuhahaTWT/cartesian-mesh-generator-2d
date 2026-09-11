'use strict';

// The theme switch is intentionally independent from mesh state and viewport
// displayMode. It can be loaded before or after app.js and remains harmless when
// the host has not added the optional selector yet.
(function () {
  const STORAGE_KEY = 'cartmesh2d-theme';
  const THEMES = new Set(['modern', 'duet']);

  function validTheme(value) {
    return THEMES.has(value) ? value : 'modern';
  }

  function readSavedTheme() {
    try {
      return validTheme(window.localStorage.getItem(STORAGE_KEY));
    } catch (_) {
      return 'modern';
    }
  }

  function saveTheme(theme) {
    try { window.localStorage.setItem(STORAGE_KEY, theme); } catch (_) { /* private mode */ }
  }

  function applyTheme(value, persist = true) {
    const theme = validTheme(value);
    document.documentElement.dataset.theme = theme;
    if (document.body) document.body.dataset.theme = theme;
    if (persist) saveTheme(theme);
    const selector = document.getElementById('appTheme');
    if (selector && selector.value !== theme) selector.value = theme;
    window.dispatchEvent(new CustomEvent('cartmesh2d-themechange', {
      detail: { theme }
    }));
    return theme;
  }

  function init() {
    const selector = document.getElementById('appTheme');
    const theme = applyTheme(readSavedTheme(), false);
    if (!selector) return;

    // A native select is already keyboard operable. These labels keep the
    // control discoverable to VoiceOver and other accessibility trees.
    selector.setAttribute('aria-label', '界面主题');
    selector.setAttribute('title', '选择界面主题');
    selector.dataset.themeControl = 'true';
    selector.value = theme;
    selector.addEventListener('change', event => applyTheme(event.target.value));
  }

  window.CartMeshTheme = Object.freeze({
    current: () => validTheme(document.documentElement.dataset.theme),
    set: value => applyTheme(value),
    storageKey: STORAGE_KEY
  });

  if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', init, { once: true });
  } else {
    init();
  }
})();
