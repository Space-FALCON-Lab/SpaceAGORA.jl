// Ensemble mode: many samples of the same configuration shown together.
// Samples are coloured by a per-sample scalar (viridis), every sample's full
// history can be drawn as a faint "spaghetti" line, and a selector walks
// through the samples. The frame payload already carries the samples as
// pseudo-spacecraft (sample i, spacecraft k -> index (i-1)*perSample + k) on a
// common time axis with NaN where a sample has not started or has ended.
import * as THREE from 'three';

import { viridis as viridisRgb } from 'viewer/colormaps.js';

export function viridis(t) { const c = viridisRgb(t); return new THREE.Color(c[0], c[1], c[2]); }

const SPAGHETTI_MAX_POINTS = 600;

export function createEnsemble(spec, frames, craft, container, state, options = {}) {
  const samples = spec.samples;
  const perSample = spec.spacecraft_per_sample || 1;
  const n = samples.length;
  const S = frames.sats;

  // Colours: scalar-mapped when the manifest carries finite scalars, else by index.
  const scalars = samples.map((s) => (Number.isFinite(s.scalar) ? s.scalar : NaN));
  const finite = scalars.filter(Number.isFinite);
  const hasScalar = finite.length > 0;
  const lo = hasScalar ? Math.min(...finite) : 0, hi = hasScalar ? Math.max(...finite) : 1;
  const sampleColors = samples.map((s, i) => {
    if (hasScalar && Number.isFinite(scalars[i])) return viridis(hi > lo ? (scalars[i] - lo) / (hi - lo) : 0.5);
    if (hasScalar) return new THREE.Color(0x777777);
    return new THREE.Color().setHSL((i / Math.max(1, n)) % 1, 0.7, 0.6);
  });
  const colors = new Array(S);
  for (let s = 0; s < S; s++) colors[s] = sampleColors[Math.floor(s / perSample)];
  craft.setColors(colors);

  // Spaghetti: one faint polyline per pseudo-spacecraft over its whole history.
  const group = new THREE.Group();
  group.name = 'ensemble-spaghetti';
  const lines = [];
  const step = Math.max(1, Math.ceil(frames.count / SPAGHETTI_MAX_POINTS));
  for (let s = 0; s < S; s++) {
    const pts = [];
    for (let k = 0; k < frames.count; k += step) {
      const b = (k * S + s) * 3;
      const x = frames.pos[b], y = frames.pos[b + 1], z = frames.pos[b + 2];
      if (Number.isFinite(x) && Number.isFinite(y) && Number.isFinite(z)) pts.push(x, y, z);
    }
    const geometry = new THREE.BufferGeometry();
    geometry.setAttribute('position', new THREE.BufferAttribute(Float32Array.from(pts), 3));
    const line = new THREE.Line(geometry, new THREE.LineBasicMaterial({ color: colors[s].clone(), transparent: true, opacity: 0.35 }));
    line.frustumCulled = false;
    lines.push(line);
    group.add(line);
  }

  const dimMask = new Uint8Array(S);
  let dimOthers = true;
  let spaghetti = true;

  // Panel: legend, selector, toggles.
  const panel = document.createElement('div');
  panel.className = 'sa-panel sa-ensemble';
  panel.style.cssText = 'left: 14px; bottom: 96px; max-width: 40ch;';
  const fmt = (v) => (Math.abs(v) >= 1000 ? v.toFixed(0) : v.toPrecision(4));
  panel.innerHTML = `
    <h1>Ensemble · ${n} samples</h1>
    ${hasScalar ? `<div class="sa-legend"><span>${fmt(lo)}</span><span class="sa-bar"></span><span>${fmt(hi)}</span></div>
    <div class="sa-hint" style="margin:0 0 6px">${spec.scalar_name || 'scalar'}</div>` : '<div class="sa-hint" style="margin:0 0 6px">coloured by sample index</div>'}
    <div class="sa-row" style="gap:6px">
      <button data-role="prev" title="previous sample">◀</button>
      <select data-role="sample" style="max-width: 24ch"></select>
      <button data-role="next" title="next sample">▶</button>
    </div>
    <div class="sa-row sa-toggle" style="margin-top:6px">
      <label><input type="checkbox" data-role="spaghetti" checked> all histories</label>
      <label><input type="checkbox" data-role="dim" checked> dim others</label>
    </div>`;
  const style = document.createElement('style');
  style.textContent = `
    .sa-legend { display: flex; align-items: center; gap: 6px; font-variant-numeric: tabular-nums; }
    .sa-bar { flex: 1; height: 10px; border-radius: 2px; background: linear-gradient(to right,
      rgb(68,1,84), rgb(59,82,139), rgb(33,145,140), rgb(94,201,98), rgb(253,231,37)); }
    .sa-ensemble .sa-row { pointer-events: auto; }
  `;
  document.head.appendChild(style);
  container.appendChild(panel);
  const sel = panel.querySelector('[data-role="sample"]');
  const none = document.createElement('option');
  none.value = '-1'; none.textContent = 'no sample selected';
  sel.appendChild(none);
  samples.forEach((s, i) => {
    const opt = document.createElement('option');
    opt.value = String(i);
    const tag = hasScalar && Number.isFinite(scalars[i]) ? ` · ${fmt(scalars[i])}` : '';
    opt.textContent = `${s.label || `sample ${s.index}`}${s.success ? '' : ' (failed)'}${tag}`;
    sel.appendChild(opt);
  });
  sel.addEventListener('change', () => state.select(Number(sel.value) < 0 ? -1 : Number(sel.value) * perSample));
  panel.querySelector('[data-role="prev"]').addEventListener('click', () => stepSample(-1));
  panel.querySelector('[data-role="next"]').addEventListener('click', () => stepSample(1));
  panel.querySelector('[data-role="spaghetti"]').addEventListener('change', (e) => { spaghetti = e.target.checked; group.visible = spaghetti; });
  panel.querySelector('[data-role="dim"]').addEventListener('change', (e) => { dimOthers = e.target.checked; });

  function currentSample() { return state.selected < 0 ? -1 : Math.floor(state.selected / perSample); }
  function stepSample(d) {
    const cur = currentSample();
    const next = cur < 0 ? (d > 0 ? 0 : n - 1) : (cur + d + n) % n;
    state.select(next * perSample);
  }

  return {
    group,
    lines,
    sampleColors,
    perSample,
    dimMask,
    sampleOf(index) { return index < 0 ? -1 : Math.floor(index / perSample); },
    // Called every frame: dim mask for markers, highlight the selected sample's history.
    update(follow) {
      const cur = currentSample();
      sel.value = String(cur);
      for (let s = 0; s < S; s++) {
        const mine = cur >= 0 && Math.floor(s / perSample) === cur;
        dimMask[s] = (dimOthers && cur >= 0 && !mine) ? 1 : 0;
        const line = lines[s];
        line.material.opacity = cur < 0 ? 0.35 : (mine ? 0.9 : (dimOthers ? 0.08 : 0.35));
        // The followed sample's own history is drawn by its trail at full precision.
        line.visible = spaghetti && !(follow && mine);
      }
    },
  };
}
