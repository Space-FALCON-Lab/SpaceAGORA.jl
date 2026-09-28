// Time-history plots. Every quantity in the selection panel, the face panel
// and the ensemble panel is clickable; the click opens this panel with the
// quantity's history over the run drawn on a canvas, the playback cursor on
// it, and the current value in the header. Clicking or dragging in the plot
// seeks the timeline. One quantity shows at a time; clicking the same
// quantity again closes the panel.
//
// A plot spec is { key, title, unit, t: Float64Array, series: [{ name, y: Float64Array }], log }.
// `log` is true, false or 'auto' (log axis when every value is positive and
// the history spans more than two and a half decades).
import { elapsedString } from 'viewer/timeline.js';

const SERIES_COLORS = ['#f2b950', '#7fd1ff', '#8bd88b', '#ff8c69'];
const MARGIN = { left: 60, right: 12, top: 10, bottom: 24 };

// Largest "nice" step (1, 2, 5 x 10^n) that gives at most `target` intervals over `span`.
function niceStep(span, target) {
  const raw = span / Math.max(1, target);
  const p = Math.pow(10, Math.floor(Math.log10(raw)));
  for (const m of [1, 2, 5, 10]) if (m * p >= raw) return m * p;
  return 10 * p;
}

// Elapsed-time axis unit for a run of `spanS` seconds.
function timeUnit(spanS) {
  if (spanS >= 3 * 86400) return { label: 'd', scale: 1 / 86400 };
  if (spanS >= 3 * 3600) return { label: 'h', scale: 1 / 3600 };
  if (spanS >= 3 * 60) return { label: 'min', scale: 1 / 60 };
  return { label: 's', scale: 1 };
}

export function formatQuantity(v, digits = 4) {
  if (!Number.isFinite(v)) return '–';
  if (v === 0) return '0';
  const a = Math.abs(v);
  if (a >= 1e5 || a < 1e-3) return v.toExponential(digits - 1);
  return v.toPrecision(digits);
}

function tickLabel(v) {
  const a = Math.abs(v);
  if (v === 0) return '0';
  if (a >= 1e5 || a < 1e-3) return v.toExponential(1);
  if (a >= 100) return v.toFixed(0);
  if (a >= 10) return v.toFixed(1);
  if (a >= 1) return v.toFixed(2);
  return v.toFixed(3);
}

export function createPlotPanel(container, timeline, frames) {
  const panel = document.createElement('div');
  panel.className = 'sa-panel sa-plot';
  panel.hidden = true;
  panel.innerHTML = `
    <div class="sa-plot-head">
      <h1 data-role="title"></h1>
      <span class="sa-plot-value" data-role="value"></span>
      <label title="Logarithmic value axis"><input type="checkbox" data-role="log"> log</label>
      <button data-role="close" title="Close the plot">×</button>
    </div>
    <canvas data-role="canvas"></canvas>
    <p class="sa-hint">Click in the plot to seek. Click the quantity again to close.</p>`;
  container.appendChild(panel);
  const canvas = panel.querySelector('[data-role="canvas"]');
  const ctx = canvas.getContext('2d');
  const titleEl = panel.querySelector('[data-role="title"]'), valueEl = panel.querySelector('[data-role="value"]');
  const logBox = panel.querySelector('[data-role="log"]');
  panel.querySelector('[data-role="close"]').addEventListener('click', () => hide());
  logBox.addEventListener('change', () => { if (spec) { useLog = logBox.checked; layout(); draw(timeline.t, true); } });

  let spec = null, useLog = false, lastDraw = -Infinity, lastT = NaN;
  let width = 0, height = 0, dpr = 1;
  let axis = null; // { x0, x1, y0, y1, xs, ys, ticksX, ticksY, unit }

  function values(t) {
    // Sampled value of each series at t (linear between samples).
    if (!spec || spec.t.length === 0) return [];
    const T = spec.t;
    let lo = 0, hi = T.length - 1;
    if (t <= T[0]) hi = 0; else if (t >= T[hi]) lo = hi; else { while (hi - lo > 1) { const m = (lo + hi) >> 1; if (T[m] <= t) lo = m; else hi = m; } }
    const f = hi > lo ? (t - T[lo]) / (T[hi] - T[lo]) : 0;
    return spec.series.map((s) => s.y[lo] + f * (s.y[hi] - s.y[lo]));
  }

  function layout() {
    if (!spec) return;
    const rect = canvas.getBoundingClientRect();
    dpr = Math.min(window.devicePixelRatio || 1, 2);
    width = Math.max(120, rect.width); height = Math.max(80, rect.height);
    canvas.width = Math.round(width * dpr); canvas.height = Math.round(height * dpr);
    const T = spec.t;
    const t0 = T.length ? T[0] : 0, t1 = T.length ? T[T.length - 1] : 1;
    // Zeros (heating above the atmosphere, say) are gaps on a log axis; a
    // negative value rules the log axis out.
    let lo = Infinity, hi = -Infinity, loPos = Infinity, anyNegative = false, anyPositive = false;
    for (const s of spec.series) for (let k = 0; k < s.y.length; k++) {
      const v = s.y[k];
      if (!Number.isFinite(v)) continue;
      if (v < 0) anyNegative = true;
      if (v > 0) { anyPositive = true; if (v < loPos) loPos = v; }
      if (v < lo) lo = v; if (v > hi) hi = v;
    }
    if (!Number.isFinite(lo)) { lo = 0; hi = 1; }
    const logAllowed = anyPositive && !anyNegative;
    if (spec.log === 'auto' && !logBox.dataset.touched) useLog = logAllowed && hi / loPos > Math.pow(10, 2.5);
    else if (spec.log !== 'auto' && !logBox.dataset.touched) useLog = !!spec.log;
    if (useLog && !logAllowed) useLog = false;
    logBox.checked = useLog;
    logBox.disabled = !logAllowed;
    let y0, y1, ticksY;
    if (useLog) {
      y0 = Math.log10(loPos); y1 = Math.log10(hi);
      if (y1 - y0 < 1e-9) { y0 -= 0.5; y1 += 0.5; }
      const pad = 0.04 * (y1 - y0); y0 -= pad; y1 += pad;
      ticksY = [];
      const step = Math.max(1, Math.ceil((y1 - y0) / 6));
      for (let e = Math.ceil(y0); e <= Math.floor(y1); e += step) ticksY.push({ v: e, label: `1e${e}` });
    } else {
      y0 = lo; y1 = hi;
      if (y1 - y0 < 1e-12 * Math.max(1, Math.abs(y0))) { const d = Math.abs(y0) > 0 ? 0.05 * Math.abs(y0) : 0.5; y0 -= d; y1 += d; }
      const pad = 0.06 * (y1 - y0); y0 -= pad; y1 += pad;
      const step = niceStep(y1 - y0, 6);
      ticksY = [];
      for (let v = Math.ceil(y0 / step) * step; v <= y1 + 1e-9 * step; v += step) ticksY.push({ v, label: tickLabel(Math.abs(v) < 1e-12 * step ? 0 : v) });
    }
    const unit = timeUnit(t1 - t0);
    const xStep = niceStep((t1 - t0) * unit.scale, 6);
    const ticksX = [];
    for (let v = Math.ceil(t0 * unit.scale / xStep) * xStep; v <= t1 * unit.scale + 1e-9; v += xStep) ticksX.push({ v: v / unit.scale, label: `${tickLabel(v - t0 * unit.scale)} ${unit.label}` });
    const plotW = width - MARGIN.left - MARGIN.right, plotH = height - MARGIN.top - MARGIN.bottom;
    axis = {
      x0: t0, x1: Math.max(t1, t0 + 1e-9), y0, y1, unit, ticksX, ticksY,
      xs: (t) => MARGIN.left + (t - t0) / Math.max(1e-9, t1 - t0) * plotW,
      ys: (v) => { const yv = useLog ? Math.log10(Math.max(v, 1e-300)) : v; return MARGIN.top + (1 - (yv - y0) / (y1 - y0)) * plotH; },
      plotW, plotH,
    };
  }

  function draw(t, force = false) {
    if (!spec || panel.hidden || !axis) return;
    const now = performance.now();
    if (!force && now - lastDraw < 50 && t === lastT) return;
    lastDraw = now; lastT = t;
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    ctx.clearRect(0, 0, width, height);
    ctx.font = '10px ui-monospace, Menlo, monospace';
    ctx.fillStyle = 'rgba(255,255,255,0.03)';
    ctx.fillRect(MARGIN.left, MARGIN.top, axis.plotW, axis.plotH);
    // grid and ticks
    ctx.strokeStyle = 'rgba(143,161,184,0.22)'; ctx.lineWidth = 1;
    ctx.fillStyle = '#8fa1b8'; ctx.textAlign = 'right'; ctx.textBaseline = 'middle';
    for (const tk of axis.ticksY) {
      const y = Math.round(axis.ys(useLog ? Math.pow(10, tk.v) : tk.v)) + 0.5;
      ctx.beginPath(); ctx.moveTo(MARGIN.left, y); ctx.lineTo(width - MARGIN.right, y); ctx.stroke();
      ctx.fillText(tk.label, MARGIN.left - 6, y);
    }
    ctx.textAlign = 'center'; ctx.textBaseline = 'top';
    for (const tk of axis.ticksX) {
      const x = Math.round(axis.xs(tk.v)) + 0.5;
      ctx.beginPath(); ctx.moveTo(x, MARGIN.top); ctx.lineTo(x, height - MARGIN.bottom); ctx.stroke();
      ctx.fillText(tk.label, x, height - MARGIN.bottom + 5);
    }
    ctx.strokeStyle = 'rgba(143,161,184,0.6)';
    ctx.strokeRect(MARGIN.left + 0.5, MARGIN.top + 0.5, axis.plotW - 1, axis.plotH - 1);
    // series
    ctx.save();
    ctx.beginPath(); ctx.rect(MARGIN.left, MARGIN.top, axis.plotW, axis.plotH); ctx.clip();
    spec.series.forEach((s, i) => {
      ctx.strokeStyle = SERIES_COLORS[i % SERIES_COLORS.length]; ctx.lineWidth = 1.5;
      ctx.beginPath();
      let pen = false;
      for (let k = 0; k < spec.t.length; k++) {
        const v = s.y[k];
        if (!Number.isFinite(v) || (useLog && v <= 0)) { pen = false; continue; }
        const x = axis.xs(spec.t[k]), y = axis.ys(v);
        if (pen) ctx.lineTo(x, y); else { ctx.moveTo(x, y); pen = true; }
      }
      ctx.stroke();
    });
    ctx.restore();
    // cursor
    const vals = values(t);
    const cx = axis.xs(Math.min(axis.x1, Math.max(axis.x0, t)));
    ctx.strokeStyle = 'rgba(242,185,80,0.9)'; ctx.lineWidth = 1;
    ctx.beginPath(); ctx.moveTo(Math.round(cx) + 0.5, MARGIN.top); ctx.lineTo(Math.round(cx) + 0.5, height - MARGIN.bottom); ctx.stroke();
    vals.forEach((v, i) => {
      if (!Number.isFinite(v) || (useLog && v <= 0)) return;
      ctx.fillStyle = SERIES_COLORS[i % SERIES_COLORS.length];
      ctx.beginPath(); ctx.arc(cx, axis.ys(v), 3, 0, 2 * Math.PI); ctx.fill();
    });
    // legend for multi-series
    if (spec.series.length > 1) {
      ctx.textAlign = 'left'; ctx.textBaseline = 'top';
      spec.series.forEach((s, i) => {
        ctx.fillStyle = SERIES_COLORS[i % SERIES_COLORS.length];
        ctx.fillText(s.name, MARGIN.left + 6 + 34 * i, MARGIN.top + 3);
      });
    }
    valueEl.textContent = `${vals.map((v) => formatQuantity(v)).join(', ')}${spec.unit ? ' ' + spec.unit : ''} at ${elapsedString(t - (frames ? frames.tStart : 0))}`;
  }

  function seekFromEvent(e) {
    if (!axis) return;
    const rect = canvas.getBoundingClientRect();
    const f = (e.clientX - rect.left - MARGIN.left) / axis.plotW;
    timeline.seek(axis.x0 + Math.min(1, Math.max(0, f)) * (axis.x1 - axis.x0));
    draw(timeline.t, true);
  }
  let dragging = false;
  canvas.addEventListener('pointerdown', (e) => { dragging = true; canvas.setPointerCapture(e.pointerId); seekFromEvent(e); });
  canvas.addEventListener('pointermove', (e) => { if (dragging) seekFromEvent(e); });
  canvas.addEventListener('pointerup', (e) => { dragging = false; try { canvas.releasePointerCapture(e.pointerId); } catch (err) { /* ignore */ } });
  window.addEventListener('resize', () => { if (!panel.hidden) { layout(); draw(timeline.t, true); } });

  function show(next) {
    spec = next;
    logBox.dataset.touched = '';
    delete logBox.dataset.touched;
    panel.hidden = false;
    titleEl.textContent = spec.title + (spec.unit ? ` (${spec.unit})` : '');
    layout();
    draw(timeline.t, true);
  }
  function hide() { panel.hidden = true; spec = null; }
  logBox.addEventListener('click', () => { logBox.dataset.touched = '1'; });

  return {
    panel,
    show,
    hide,
    toggle(next) { if (spec && spec.key === next.key) hide(); else show(next); },
    update(t) { draw(t); },
    get key() { return spec ? spec.key : null; },
    get visible() { return !panel.hidden; },
    get spec() { return spec; },
  };
}
