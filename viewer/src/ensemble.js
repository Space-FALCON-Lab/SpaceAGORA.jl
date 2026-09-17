// Ensemble mode: many samples of the same configuration shown together.
// The nominal sample (payload.ensemble.nominal, else the first) is drawn as a
// bright white line with a labeled marker; the other samples are the
// Monte Carlo traces, colored by a per-sample scalar (viridis) or by index
// and kept faint so the nominal reads first; and a translucent 3-sigma tube
// follows the nominal (or the sample mean when there is no nominal), its
// cross-section the 3-sigma dispersion of the samples in the radial and
// cross-track directions at each time, with the along-track 3-sigma read out
// in the panel. The frame payload already carries the samples as
// pseudo-spacecraft (sample i, spacecraft k -> index (i-1)*perSample + k) on a
// common time axis with NaN where a sample has not started or has ended.
import * as THREE from 'three';

import { viridis as viridisRgb } from 'viewer/colormaps.js';
import { makeLabelSprite } from 'viewer/spacecraft.js';

export function viridis(t) { const c = viridisRgb(t); return new THREE.Color(c[0], c[1], c[2]); }

const SPAGHETTI_MAX_POINTS = 600;
const TUBE_MAX_RINGS = 400;
const TUBE_SEGMENTS = 20;
const NOMINAL_COLOR = 0xffffff;
const TUBE_COLOR = 0xffb347;

// Per-frame statistics of the first spacecraft of every sample: mean position,
// and the 3-sigma dispersion along the radial, along-track and cross-track
// axes of the reference trajectory (the nominal, or the mean).
function computeDispersion(frames, perSample, n, nominal, ringStep) {
  const S = frames.sats;
  const rings = [];
  const r = new THREE.Vector3(), v = new THREE.Vector3(), R = new THREE.Vector3(), T = new THREE.Vector3(), N = new THREE.Vector3(), d = new THREE.Vector3();
  const mean = new THREE.Vector3();
  const vel = new Float64Array(3);
  for (let k = 0; k < frames.count; k += ringStep) {
    // mean of the present samples
    mean.set(0, 0, 0); let m = 0;
    for (let i = 0; i < n; i++) {
      const b = (k * S + i * perSample) * 3;
      const x = frames.pos[b], y = frames.pos[b + 1], z = frames.pos[b + 2];
      if (!Number.isFinite(x)) continue;
      mean.x += x; mean.y += y; mean.z += z; m++;
    }
    if (m < 2) continue;
    mean.multiplyScalar(1 / m);
    // reference point and axes
    if (nominal >= 0) {
      const b = (k * S + nominal * perSample) * 3;
      if (!Number.isFinite(frames.pos[b])) continue;
      r.set(frames.pos[b], frames.pos[b + 1], frames.pos[b + 2]);
      frames.velocityAt(frames.t[k], nominal * perSample, vel);
    } else {
      r.copy(mean);
      // mean velocity
      vel[0] = vel[1] = vel[2] = 0; let mv = 0; const tmp = new Float64Array(3);
      for (let i = 0; i < n; i++) { const b = (k * S + i * perSample) * 3; if (!Number.isFinite(frames.pos[b])) continue; frames.velocityAt(frames.t[k], i * perSample, tmp); vel[0] += tmp[0]; vel[1] += tmp[1]; vel[2] += tmp[2]; mv++; }
      if (mv) { vel[0] /= mv; vel[1] /= mv; vel[2] /= mv; }
    }
    R.copy(r).normalize();
    v.set(vel[0], vel[1], vel[2]);
    N.copy(r).cross(v);
    if (N.lengthSq() === 0) N.set(0, 0, 1);
    N.normalize();
    T.copy(N).cross(R).normalize();
    // second moments about the reference point, in RTN
    let sRR = 0, sTT = 0, sNN = 0;
    for (let i = 0; i < n; i++) {
      const b = (k * S + i * perSample) * 3;
      const x = frames.pos[b], y = frames.pos[b + 1], z = frames.pos[b + 2];
      if (!Number.isFinite(x)) continue;
      d.set(x - r.x, y - r.y, z - r.z);
      const dr = d.dot(R), dt = d.dot(T), dn = d.dot(N);
      sRR += dr * dr; sTT += dt * dt; sNN += dn * dn;
    }
    const inv = 1 / Math.max(1, m - (nominal >= 0 ? 0 : 1));
    rings.push({ k, t: frames.t[k], center: r.clone(), R: R.clone(), T: T.clone(), N: N.clone(),
      sigmaR: Math.sqrt(sRR * inv), sigmaT: Math.sqrt(sTT * inv), sigmaN: Math.sqrt(sNN * inv), members: m });
  }
  return rings;
}

// Translucent tube: one ring per dispersion sample, radii 3 sigma radial and
// 3 sigma cross-track (km, clamped to a visible minimum), triangulated
// between consecutive rings. Vertices are uploaded relative to the floating
// origin on every frame.
function buildTube(rings) {
  const nr = rings.length, seg = TUBE_SEGMENTS;
  const positions = new Float32Array(nr * seg * 3);
  const indices = [];
  for (let i = 0; i < nr - 1; i++) {
    for (let j = 0; j < seg; j++) {
      const a = i * seg + j, b = i * seg + (j + 1) % seg, c = (i + 1) * seg + j, dd = (i + 1) * seg + (j + 1) % seg;
      indices.push(a, c, b, b, c, dd);
    }
  }
  const geometry = new THREE.BufferGeometry();
  geometry.setAttribute('position', new THREE.BufferAttribute(positions, 3));
  geometry.setIndex(indices);
  const material = new THREE.MeshBasicMaterial({ color: TUBE_COLOR, transparent: true, opacity: 0.18, side: THREE.DoubleSide, depthWrite: false });
  const mesh = new THREE.Mesh(geometry, material);
  mesh.frustumCulled = false;
  mesh.renderOrder = 1;
  // outline rings every few samples so the cone reads even where it is thin
  const outline = new THREE.LineSegments(new THREE.BufferGeometry(), new THREE.LineBasicMaterial({ color: TUBE_COLOR, transparent: true, opacity: 0.45 }));
  outline.frustumCulled = false;
  const outlineEvery = Math.max(1, Math.round(nr / 40));
  const outlinePos = new Float32Array(Math.ceil(nr / outlineEvery) * seg * 2 * 3);
  outline.geometry.setAttribute('position', new THREE.BufferAttribute(outlinePos, 3));
  return { mesh, outline, outlineEvery, positions, outlinePos };
}

export function createEnsemble(spec, frames, craft, container, state, options = {}) {
  const samples = spec.samples;
  const perSample = spec.spacecraft_per_sample || 1;
  const n = samples.length;
  const S = frames.sats;
  const nominal = Number.isInteger(spec.nominal) && spec.nominal >= 0 && spec.nominal < n ? spec.nominal : -1;
  const traces = [];
  for (let i = 0; i < n; i++) if (i !== nominal) traces.push(i);

  // Colors: scalar-mapped when the manifest carries finite scalars, else by index; the nominal is white.
  const scalars = samples.map((s) => (Number.isFinite(s.scalar) ? s.scalar : NaN));
  const finite = scalars.filter(Number.isFinite);
  const hasScalar = finite.length > 0;
  const lo = hasScalar ? Math.min(...finite) : 0, hi = hasScalar ? Math.max(...finite) : 1;
  const sampleColors = samples.map((s, i) => {
    if (i === nominal) return new THREE.Color(NOMINAL_COLOR);
    if (hasScalar && Number.isFinite(scalars[i])) return viridis(hi > lo ? (scalars[i] - lo) / (hi - lo) : 0.5);
    if (hasScalar) return new THREE.Color(0x777777);
    return new THREE.Color().setHSL((i / Math.max(1, n)) % 1, 0.7, 0.6);
  });
  const colors = new Array(S);
  for (let s = 0; s < S; s++) colors[s] = sampleColors[Math.floor(s / perSample)];
  craft.setColors(colors);

  // Spaghetti: one faint polyline per pseudo-spacecraft over its whole history; the nominal's is bright.
  const group = new THREE.Group();
  group.name = 'ensemble';
  const lines = [];
  const step = Math.max(1, Math.ceil(frames.count / SPAGHETTI_MAX_POINTS));
  const lineSources = [];
  for (let s = 0; s < S; s++) {
    const pts = [];
    for (let k = 0; k < frames.count; k += step) {
      const b = (k * S + s) * 3;
      const x = frames.pos[b], y = frames.pos[b + 1], z = frames.pos[b + 2];
      if (Number.isFinite(x) && Number.isFinite(y) && Number.isFinite(z)) pts.push(x, y, z);
    }
    const src = Float64Array.from(pts);
    lineSources.push(src);
    const geometry = new THREE.BufferGeometry();
    geometry.setAttribute('position', new THREE.BufferAttribute(new Float32Array(src.length), 3));
    const isNominal = Math.floor(s / perSample) === nominal;
    const line = new THREE.Line(geometry, new THREE.LineBasicMaterial({ color: colors[s].clone(), transparent: true, opacity: isNominal ? 1.0 : 0.3 }));
    line.frustumCulled = false;
    line.renderOrder = isNominal ? 3 : 0;
    lines.push(line);
    group.add(line);
  }
  const nominalLabel = nominal >= 0 ? makeLabelSprite('nominal', '#ffffff') : null;
  if (nominalLabel) group.add(nominalLabel);

  // 3-sigma tube around the reference trajectory.
  const ringStep = Math.max(1, Math.ceil(frames.count / TUBE_MAX_RINGS));
  const rings = n >= 2 ? computeDispersion(frames, perSample, n, nominal, ringStep) : [];
  const tube = rings.length >= 2 ? buildTube(rings) : null;
  if (tube) { group.add(tube.mesh); group.add(tube.outline); }
  let tubeVisible = true;
  let tubeScale = 3;   // sigma multiple

  const dimMask = new Uint8Array(S);
  let dimOthers = true;
  let spaghetti = true;
  // Many traces make their labels a cloud: keep only the nominal's (and the selected sample's) label.
  const labelMask = new Uint8Array(S);
  if (n > 6) for (let s = 0; s < S; s++) labelMask[s] = Math.floor(s / perSample) === nominal ? 0 : 1;
  craft.setLabelMask(labelMask);

  // Panel: legend, selector, toggles, dispersion readout.
  const panel = document.createElement('div');
  panel.className = 'sa-panel sa-ensemble';
  panel.style.cssText = 'right: 14px; top: 46%; max-width: 42ch;';
  const fmt = (v) => (Math.abs(v) >= 1000 ? v.toFixed(0) : v.toPrecision(4));
  panel.innerHTML = `
    <h1>Ensemble · ${n} samples${nominal >= 0 ? ' + nominal' : ''}</h1>
    ${hasScalar ? `<div class="sa-legend"><span>${fmt(lo)}</span><span class="sa-bar"></span><span>${fmt(hi)}</span></div>
    <div class="sa-hint" data-role="scalar-label" style="margin:0 0 6px"></div>` : '<div class="sa-hint" style="margin:0 0 6px">traces colored by sample index</div>'}
    <div class="sa-key"><span class="sa-swatch" style="background:#fff"></span> nominal <span class="sa-swatch" style="background:#ffb347;opacity:.6;margin-left:8px"></span> <span data-role="tube-label">3σ tube</span></div>
    <dl class="sa-disp"><dt data-key="sr" title="Click for the time history">3σ radial</dt><dd data-key="sr" data-role="sr">–</dd><dt data-key="st" title="Click for the time history">3σ along-track</dt><dd data-key="st" data-role="st">–</dd><dt data-key="sn" title="Click for the time history">3σ cross-track</dt><dd data-key="sn" data-role="sn">–</dd><dt data-key="members" title="Click for the time history">samples present</dt><dd data-key="members" data-role="members">–</dd></dl>
    <div class="sa-row" style="gap:6px">
      <button data-role="prev" title="previous sample">◀</button>
      <select data-role="sample" style="max-width: 24ch"></select>
      <button data-role="next" title="next sample">▶</button>
    </div>
    <div class="sa-row sa-toggle" style="margin-top:6px">
      <label><input type="checkbox" data-role="spaghetti" checked> traces</label>
      <label><input type="checkbox" data-role="dim" checked> dim others</label>
      <label><input type="checkbox" data-role="tube" checked> tube</label>
      <label>σ× <select data-role="sigma"><option value="1">1</option><option value="2">2</option><option value="3" selected>3</option></select></label>
    </div>`;
  if (hasScalar) panel.querySelector('[data-role="scalar-label"]').textContent = spec.scalar_name || 'scalar';
  const style = document.createElement('style');
  style.textContent = `
    .sa-legend { display: flex; align-items: center; gap: 6px; font-variant-numeric: tabular-nums; }
    .sa-bar { flex: 1; height: 10px; border-radius: 2px; background: linear-gradient(to right,
      rgb(68,1,84), rgb(59,82,139), rgb(33,145,140), rgb(94,201,98), rgb(253,231,37)); }
    .sa-ensemble .sa-row { pointer-events: auto; }
    .sa-key { display: flex; align-items: center; gap: 6px; margin: 4px 0 6px; }
    .sa-swatch { display: inline-block; width: 14px; height: 4px; border-radius: 2px; }
    .sa-disp { margin: 0 0 6px; display: grid; grid-template-columns: max-content 1fr; gap: 1px 10px; font-variant-numeric: tabular-nums; }
    .sa-disp dt { color: #8fa1b8; } .sa-disp dd { margin: 0; }
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
    opt.textContent = `${i === nominal ? 'nominal' : (s.label || `sample ${s.index}`)}${s.success ? '' : ' (failed)'}${tag}`;
    sel.appendChild(opt);
  });
  sel.addEventListener('change', () => state.select(Number(sel.value) < 0 ? -1 : Number(sel.value) * perSample));
  panel.querySelector('[data-role="prev"]').addEventListener('click', () => stepSample(-1));
  panel.querySelector('[data-role="next"]').addEventListener('click', () => stepSample(1));
  panel.querySelector('[data-role="spaghetti"]').addEventListener('change', (e) => { spaghetti = e.target.checked; });
  panel.querySelector('[data-role="dim"]').addEventListener('change', (e) => { dimOthers = e.target.checked; });
  panel.querySelector('[data-role="tube"]').addEventListener('change', (e) => { tubeVisible = e.target.checked; });
  panel.querySelector('[data-role="sigma"]').addEventListener('change', (e) => { tubeScale = Number(e.target.value); panel.querySelector('[data-role="tube-label"]').textContent = `${tubeScale}σ tube`; tubeDirty = true; });
  const readout = { sr: panel.querySelector('[data-role="sr"]'), st: panel.querySelector('[data-role="st"]'), sn: panel.querySelector('[data-role="sn"]'), members: panel.querySelector('[data-role="members"]') };
  // Dispersion readouts open their history over the run (state.plotSeries is wired by main.js).
  panel.querySelector('.sa-disp').addEventListener('click', (e) => {
    const el = e.target.closest('[data-key]');
    if (!el || !state.plotSeries || rings.length === 0) return;
    const key = el.dataset.key;
    const t = Float64Array.from(rings, (r) => r.t);
    const pick = { sr: ['sigmaR', 'radial'], st: ['sigmaT', 'along-track'], sn: ['sigmaN', 'cross-track'] }[key];
    if (pick) {
      state.plotSeries({ key: `ensemble:${key}`, title: `${tubeScale}σ ${pick[1]}`, unit: 'km', t, series: [{ name: '', y: Float64Array.from(rings, (r) => tubeScale * r[pick[0]]) }], log: 'auto' });
    } else {
      state.plotSeries({ key: 'ensemble:members', title: 'samples present', unit: '', t, series: [{ name: '', y: Float64Array.from(rings, (r) => r.members) }], log: false });
    }
  });

  function currentSample() { return state.selected < 0 ? -1 : Math.floor(state.selected / perSample); }
  function stepSample(d) {
    const cur = currentSample();
    const next = cur < 0 ? (d > 0 ? 0 : n - 1) : (cur + d + n) % n;
    state.select(next * perSample);
  }

  const anchor = new Float64Array(3);
  let lastAnchor = [NaN, NaN, NaN];
  let tubeDirty = true;
  const p = new THREE.Vector3(), q = new THREE.Vector3();
  const fmtKm = (km) => (km >= 10 ? `${km.toFixed(1)} km` : km >= 0.01 ? `${(1000 * km).toFixed(0)} m` : `${(1000 * km).toFixed(2)} m`);

  function uploadLines() {
    for (let s = 0; s < S; s++) {
      const src = lineSources[s], arr = lines[s].geometry.getAttribute('position').array;
      for (let i = 0; i < src.length; i += 3) { arr[i] = src[i] - anchor[0]; arr[i + 1] = src[i + 1] - anchor[1]; arr[i + 2] = src[i + 2] - anchor[2]; }
      lines[s].geometry.getAttribute('position').needsUpdate = true;
      lines[s].geometry.setDrawRange(0, src.length / 3);
    }
  }

  function uploadTube() {
    if (!tube) return;
    const seg = TUBE_SEGMENTS, pos = tube.positions, out = tube.outlinePos;
    let o = 0;
    const minKm = 0.0005; // half a meter so a converged ensemble still shows a thread
    rings.forEach((ring, i) => {
      const a = Math.max(tubeScale * ring.sigmaR, minKm), b = Math.max(tubeScale * ring.sigmaN, minKm);
      for (let j = 0; j < seg; j++) {
        const th = 2 * Math.PI * j / seg;
        p.copy(ring.center).addScaledVector(ring.R, a * Math.cos(th)).addScaledVector(ring.N, b * Math.sin(th));
        const idx = (i * seg + j) * 3;
        pos[idx] = p.x - anchor[0]; pos[idx + 1] = p.y - anchor[1]; pos[idx + 2] = p.z - anchor[2];
      }
      if (i % tube.outlineEvery === 0 && o + seg * 6 <= out.length) {
        for (let j = 0; j < seg; j++) {
          const th0 = 2 * Math.PI * j / seg, th1 = 2 * Math.PI * ((j + 1) % seg) / seg;
          p.copy(ring.center).addScaledVector(ring.R, a * Math.cos(th0)).addScaledVector(ring.N, b * Math.sin(th0));
          q.copy(ring.center).addScaledVector(ring.R, a * Math.cos(th1)).addScaledVector(ring.N, b * Math.sin(th1));
          out[o++] = p.x - anchor[0]; out[o++] = p.y - anchor[1]; out[o++] = p.z - anchor[2];
          out[o++] = q.x - anchor[0]; out[o++] = q.y - anchor[1]; out[o++] = q.z - anchor[2];
        }
      }
    });
    tube.mesh.geometry.getAttribute('position').needsUpdate = true;
    tube.outline.geometry.getAttribute('position').needsUpdate = true;
    tube.outline.geometry.setDrawRange(0, o / 3);
    tube.mesh.geometry.computeBoundingSphere();
  }

  function ringAt(t) {
    if (rings.length === 0) return null;
    let lo = 0, hi = rings.length - 1;
    while (hi - lo > 1) { const mid = (lo + hi) >> 1; if (rings[mid].t <= t) lo = mid; else hi = mid; }
    return rings[t - rings[lo].t < rings[hi].t - t ? lo : hi];
  }

  return {
    group,
    lines,
    sampleColors,
    perSample,
    nominal,
    dimMask,
    rings,
    sampleOf(index) { return index < 0 ? -1 : Math.floor(index / perSample); },
    // Called every frame: dim mask for markers, highlight the selected sample's history, place the tube.
    update(follow, t = frames.tStart, anchorKm = null, camera = null, groupMatrix = null) {
      const cur = currentSample();
      sel.value = String(cur);
      if (anchorKm) { anchor[0] = anchorKm[0]; anchor[1] = anchorKm[1]; anchor[2] = anchorKm[2]; }
      group.position.set(anchor[0], anchor[1], anchor[2]);
      const moved = anchor[0] !== lastAnchor[0] || anchor[1] !== lastAnchor[1] || anchor[2] !== lastAnchor[2];
      if (moved) { lastAnchor = [anchor[0], anchor[1], anchor[2]]; uploadLines(); tubeDirty = true; }
      if (tubeDirty) { uploadTube(); tubeDirty = false; }
      for (let s = 0; s < S; s++) {
        const sample = Math.floor(s / perSample);
        const mine = cur >= 0 && sample === cur;
        const isNominal = sample === nominal;
        dimMask[s] = isNominal ? 2 : ((dimOthers && cur >= 0 && !mine) ? 1 : 0);
        if (n > 6) labelMask[s] = (isNominal || mine) ? 0 : 1;
        const line = lines[s];
        line.material.opacity = isNominal ? 1.0 : (cur < 0 ? 0.3 : (mine ? 0.9 : (dimOthers ? 0.06 : 0.3)));
        // The followed sample's own history is drawn by its trail at full precision.
        line.visible = spaghetti && !(follow && mine);
      }
      if (tube) { tube.mesh.visible = tubeVisible; tube.outline.visible = tubeVisible; }
      if (nominalLabel) {
        const b = nominal * perSample;
        const posN = new Float64Array(3);
        frames.positionAt(t, b, posN);
        nominalLabel.visible = Number.isFinite(posN[0]);
        if (nominalLabel.visible) {
          nominalLabel.position.set(posN[0] - anchor[0], posN[1] - anchor[1], posN[2] - anchor[2]);
          if (camera && groupMatrix) {
            const w = nominalLabel.position.clone().applyMatrix4(groupMatrix);
            const h = 0.032 * camera.position.distanceTo(w);
            nominalLabel.scale.set(h * nominalLabel.userData.aspect, h, 1);
          }
        }
      }
      const ring = ringAt(t);
      if (ring) {
        readout.sr.textContent = fmtKm(tubeScale * ring.sigmaR);
        readout.st.textContent = fmtKm(tubeScale * ring.sigmaT);
        readout.sn.textContent = fmtKm(tubeScale * ring.sigmaN);
        readout.members.textContent = `${ring.members} of ${n}`;
      }
    },
  };
}
