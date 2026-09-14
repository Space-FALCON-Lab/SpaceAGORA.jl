// Far-view spacecraft rendering: one point cloud for every marker (constant
// pixel size, one draw call for tens of thousands of members), per-spacecraft
// trails that fade with age for small runs, and name labels.
import * as THREE from 'three';
import { inferno } from 'viewer/colormaps.js';

export const TRAIL_MAX_SPACECRAFT = 64;
// Trail colour modes: 'age' fades with time; the others map a diagnostic
// through inferno over the run's range (log scale for the physical ones).
export const TRAIL_COLOR_MODES = {
  age: { label: 'age', log: false },
  heat_rate: { label: 'heat rate (W/m²)', log: true },
  dynamic_pressure: { label: 'dynamic pressure (Pa)', log: true },
  density: { label: 'density (kg/m³)', log: true },
  altitude: { label: 'altitude (km)', log: false },
  speed: { label: 'speed (km/s)', log: false },
};
export const LABEL_MAX_SPACECRAFT = 64;
const TRAIL_MAX_POINTS = 8192;

export function spacecraftColor(index, count) {
  const hue = count <= 1 ? 0.12 : (index / count) % 1;
  return new THREE.Color().setHSL(hue, 0.85, 0.6);
}

// Median spacing of radius minima along the saved trajectory: the orbital
// period without needing the gravitational parameter. Falls back to the run
// span when fewer than two periapsis passages were saved.
export function estimateOrbitPeriod(frames, sat = 0) {
  const n = frames.count, S = frames.sats;
  if (n < 5) return frames.tEnd - frames.tStart;
  const r = new Float64Array(n);
  for (let k = 0; k < n; k++) {
    const b = (k * S + sat) * 3;
    r[k] = Math.hypot(frames.pos[b], frames.pos[b + 1], frames.pos[b + 2]);
  }
  const minima = [];
  for (let k = 1; k < n - 1; k++) {
    if (r[k] <= r[k - 1] && r[k] < r[k + 1]) minima.push(frames.t[k]);
  }
  if (minima.length < 2) return frames.tEnd - frames.tStart;
  const gaps = [];
  for (let k = 1; k < minima.length; k++) gaps.push(minima[k] - minima[k - 1]);
  gaps.sort((a, b) => a - b);
  return gaps[gaps.length >> 1];
}

// `common` supplies EPSILON and isPerspectiveMatrix, which the log-depth
// chunks use; without it the shader fails to compile and every marker vanishes.
const POINT_VERTEX = `
  #include <common>
  attribute vec3 pointColor;
  attribute float state;
  uniform float size;
  uniform float pixelRatio;
  varying vec3 vColor;
  varying float vState;
  #include <logdepthbuf_pars_vertex>
  void main() {
    vColor = pointColor;
    vState = state;
    vec4 mvPosition = modelViewMatrix * vec4(position, 1.0);
    gl_Position = projectionMatrix * mvPosition;
    gl_PointSize = size * pixelRatio * (state > 1.5 && state < 2.5 ? 1.9 : 1.0);
    #include <logdepthbuf_vertex>
  }`;

const POINT_FRAGMENT = `
  #include <common>
  varying vec3 vColor;
  varying float vState;
  #include <logdepthbuf_pars_fragment>
  void main() {
    if (vState > 0.5 && vState < 1.5) discard;
    vec2 c = gl_PointCoord - vec2(0.5);
    float d = dot(c, c);
    if (d > 0.25) discard;
    #include <logdepthbuf_fragment>
    vec3 color = vColor;
    if (vState > 1.5 && vState < 2.5) {
      // selected: white ring around the coloured core
      if (d > 0.14) color = vec3(1.0);
    } else if (vState > 2.5) {
      // dimmed (ensemble: not the selected sample)
      color *= 0.3;
    }
    gl_FragColor = vec4(color, 1.0);
  }`;

export function makeLabelSprite(text, color) {
  const canvas = document.createElement('canvas');
  const ctx = canvas.getContext('2d');
  const font = '500 26px system-ui, sans-serif';
  ctx.font = font;
  const w = Math.ceil(ctx.measureText(text).width) + 24, h = 40;
  canvas.width = w; canvas.height = h;
  ctx.font = font;
  ctx.fillStyle = 'rgba(8,12,20,0.55)';
  ctx.fillRect(0, 0, w, h);
  ctx.fillStyle = color;
  ctx.textBaseline = 'middle';
  ctx.fillText(text, 12, h / 2);
  const texture = new THREE.CanvasTexture(canvas);
  texture.colorSpace = THREE.SRGBColorSpace;
  const sprite = new THREE.Sprite(new THREE.SpriteMaterial({ map: texture, transparent: true, depthTest: false }));
  sprite.userData.aspect = w / h;
  sprite.center.set(-0.15, 0.5);
  return sprite;
}

export function createSpacecraft(frames, sidecar, options = {}) {
  const S = frames.sats;
  const group = new THREE.Group();
  group.name = 'spacecraft';

  const positions = new Float32Array(S * 3);
  const colors = new Float32Array(S * 3);
  const state = new Float32Array(S);
  const colorObjects = [];
  for (let s = 0; s < S; s++) {
    const c = spacecraftColor(s, S);
    colorObjects.push(c);
    colors[3 * s] = c.r; colors[3 * s + 1] = c.g; colors[3 * s + 2] = c.b;
  }
  const markerGeometry = new THREE.BufferGeometry();
  markerGeometry.setAttribute('position', new THREE.BufferAttribute(positions, 3));
  markerGeometry.setAttribute('pointColor', new THREE.BufferAttribute(colors, 3));
  markerGeometry.setAttribute('state', new THREE.BufferAttribute(state, 1));
  const markerMaterial = new THREE.ShaderMaterial({
    uniforms: { size: { value: options.markerPixels ?? 7 }, pixelRatio: { value: options.pixelRatio ?? 1 } },
    vertexShader: POINT_VERTEX,
    fragmentShader: POINT_FRAGMENT,
    transparent: false,
  });
  const markers = new THREE.Points(markerGeometry, markerMaterial);
  markers.frustumCulled = false;
  group.add(markers);

  // Trails: window in seconds, chosen from an orbit count by main.js.
  let trailSeconds = options.trailSeconds ?? 0;
  const trailsEnabled = S <= TRAIL_MAX_SPACECRAFT;
  const trails = [];
  if (trailsEnabled) {
    for (let s = 0; s < S; s++) {
      const geometry = new THREE.BufferGeometry();
      geometry.setAttribute('position', new THREE.BufferAttribute(new Float32Array(TRAIL_MAX_POINTS * 3), 3));
      geometry.setAttribute('color', new THREE.BufferAttribute(new Float32Array(TRAIL_MAX_POINTS * 3), 3));
      geometry.setDrawRange(0, 0);
      const line = new THREE.Line(geometry, new THREE.LineBasicMaterial({ vertexColors: true, transparent: true, opacity: 0.9 }));
      line.frustumCulled = false;
      trails.push(line);
      group.add(line);
    }
  }

  const labelsEnabled = S <= LABEL_MAX_SPACECRAFT;
  const labels = [];
  if (labelsEnabled) {
    for (let s = 0; s < S; s++) {
      const name = (sidecar.spacecraft[s] && sidecar.spacecraft[s].name) || `sc${s + 1}`;
      const sprite = makeLabelSprite(name, `#${colorObjects[s].getHexString()}`);
      labels.push(sprite);
      group.add(sprite);
    }
  }

  const head = new Float64Array(3);
  const dark = new THREE.Color(0x0a0f18);
  const mixed = new THREE.Color();
  const Re = options.planetRadiusKm ?? 0;
  let colorMode = 'age';
  let colorRange = { lo: 0, hi: 1, log: false };

  // Range of a diagnostic over the whole run (positive values only on log scales).
  function scalarRange(kind) {
    const log = TRAIL_COLOR_MODES[kind].log;
    let lo = Infinity, hi = -Infinity;
    for (let k = 0; k < frames.count; k++) {
      for (let s = 0; s < S; s++) {
        const v = frames.scalarAtIndex(kind, k, s, Re);
        if (!Number.isFinite(v) || (log && v <= 0)) continue;
        const x = log ? Math.log10(v) : v;
        if (x < lo) lo = x; if (x > hi) hi = x;
      }
    }
    if (!(lo < hi)) { lo = 0; hi = 1; }
    if (log && hi - lo > 8) lo = hi - 8; // eight decades is plenty for a colour scale
    return { lo, hi, log };
  }

  function scalarColor(value, out) {
    const { lo, hi, log } = colorRange;
    const x = log ? (value > 0 ? Math.log10(value) : lo) : value;
    const c = inferno((x - lo) / (hi - lo));
    return out.setRGB(c[0], c[1], c[2]);
  }
  // Positions are uploaded relative to `anchor` (km, Float64), and `group` is
  // placed at the anchor. Near the followed spacecraft the buffers then hold
  // small numbers, so Float32 quantisation stays at the millimetre level
  // instead of the half-metre it would be at planetary distances.
  const anchor = new Float64Array(3);

  function updateTrails(t) {
    if (trailSeconds <= 0) { for (const l of trails) l.geometry.setDrawRange(0, 0); return; }
    const t0 = Math.max(frames.tStart, t - trailSeconds);
    const { i: iEnd } = frames.locate(t);
    const { i: iStart } = frames.locate(t0);
    const span = Math.max(0, iEnd - iStart + 1);
    const step = Math.max(1, Math.ceil(span / (TRAIL_MAX_POINTS - 1)));
    for (let s = 0; s < S; s++) {
      const posAttr = trails[s].geometry.getAttribute('position');
      const colAttr = trails[s].geometry.getAttribute('color');
      const arr = posAttr.array, col = colAttr.array;
      const c = colorObjects[s];
      let n = 0;
      for (let k = iStart; k <= iEnd && k < frames.count; k += step) {
        if (frames.t[k] > t) break;
        const base = (k * S + s) * 3;
        if (!Number.isFinite(frames.pos[base])) continue;
        arr[3 * n] = frames.pos[base] - anchor[0]; arr[3 * n + 1] = frames.pos[base + 1] - anchor[1]; arr[3 * n + 2] = frames.pos[base + 2] - anchor[2];
        if (colorMode === 'age') {
          // fade: oldest samples sink toward the background colour
          const age = (t - frames.t[k]) / trailSeconds;
          mixed.copy(c).lerp(dark, Math.min(1, Math.max(0, 0.15 + 0.8 * age)));
        } else {
          scalarColor(frames.scalarAtIndex(colorMode, k, s, Re), mixed);
        }
        col[3 * n] = mixed.r; col[3 * n + 1] = mixed.g; col[3 * n + 2] = mixed.b;
        n++;
      }
      frames.positionAt(t, s, head);
      if (Number.isFinite(head[0])) {
        arr[3 * n] = head[0] - anchor[0]; arr[3 * n + 1] = head[1] - anchor[1]; arr[3 * n + 2] = head[2] - anchor[2];
        col[3 * n] = c.r; col[3 * n + 1] = c.g; col[3 * n + 2] = c.b;
        n++;
      }
      trails[s].geometry.setDrawRange(0, n);
      posAttr.needsUpdate = true;
      colAttr.needsUpdate = true;
    }
  }

  const labelWorld = new THREE.Vector3();
  const scratch = new Float64Array(S * 3);
  const absent = new Uint8Array(S);
  let labelsWanted = true;

  return {
    group,
    markers,
    trails,
    labels,
    trailsEnabled,
    labelsEnabled,
    positions,
    colors: colorObjects,
    get trailSeconds() { return trailSeconds; },
    setTrailSeconds(v) { trailSeconds = Math.max(0, v); },
    get trailColorMode() { return colorMode; },
    // Returns the legend {label, lo, hi, log} (values in natural units), or null for 'age'.
    setTrailColorMode(mode) {
      if (!TRAIL_COLOR_MODES[mode] || (mode !== 'age' && !frames.hasScalar(mode))) mode = 'age';
      colorMode = mode;
      if (mode === 'age') return null;
      colorRange = scalarRange(mode);
      const { lo, hi, log } = colorRange;
      return { label: TRAIL_COLOR_MODES[mode].label, lo: log ? Math.pow(10, lo) : lo, hi: log ? Math.pow(10, hi) : hi, log };
    },
    availableColorModes() { return Object.keys(TRAIL_COLOR_MODES).filter((m) => m === 'age' || frames.hasScalar(m)); },
    // markerHidden: per-spacecraft flags from the LOD layer; selected: index or -1;
    // groupMatrix: this group's matrixWorld; anchorKm: Float64 [x, y, z] the group sits at;
    // dimMask: optional per-spacecraft flags (ensemble: not the selected sample)
    update(t, camera, markerHidden, selected, groupMatrix, anchorKm, dimMask) {
      anchor[0] = anchorKm[0]; anchor[1] = anchorKm[1]; anchor[2] = anchorKm[2];
      group.position.set(anchor[0], anchor[1], anchor[2]);
      frames.positionsAt(t, scratch);
      for (let k = 0; k < S; k++) {
        const x = scratch[3 * k], y = scratch[3 * k + 1], z = scratch[3 * k + 2];
        const present = Number.isFinite(x) && Number.isFinite(y) && Number.isFinite(z);
        positions[3 * k] = present ? x - anchor[0] : 0;
        positions[3 * k + 1] = present ? y - anchor[1] : 0;
        positions[3 * k + 2] = present ? z - anchor[2] : 0;
        absent[k] = present ? 0 : 1;
      }
      markerGeometry.getAttribute('position').needsUpdate = true;
      for (let s = 0; s < S; s++) {
        state[s] = absent[s] || (markerHidden && markerHidden[s]) ? 1 : (s === selected ? 2 : (dimMask && dimMask[s] ? 3 : 0));
      }
      markerGeometry.getAttribute('state').needsUpdate = true;
      if (trailsEnabled) updateTrails(t);
      if (labelsEnabled && labelsWanted) {
        for (let s = 0; s < S; s++) {
          const sprite = labels[s];
          sprite.visible = !absent[s];
          sprite.position.set(positions[3 * s], positions[3 * s + 1], positions[3 * s + 2]);
          labelWorld.copy(sprite.position).applyMatrix4(groupMatrix);
          const dist = camera.position.distanceTo(labelWorld);
          const h = 0.028 * dist; // constant on-screen height
          sprite.scale.set(h * sprite.userData.aspect, h, 1);
        }
      }
    },
    // Replace every spacecraft colour (ensemble mode); trails follow, labels keep their text.
    setColors(list) {
      for (let s = 0; s < S; s++) {
        const c = list[s] || colorObjects[s];
        colorObjects[s] = c;
        colors[3 * s] = c.r; colors[3 * s + 1] = c.g; colors[3 * s + 2] = c.b;
      }
      markerGeometry.getAttribute('pointColor').needsUpdate = true;
    },
    setMarkerPixels(px) { markerMaterial.uniforms.size.value = px; },
    setPixelRatio(r) { markerMaterial.uniforms.pixelRatio.value = r; },
    setTrailsVisible(v) { for (const l of trails) l.visible = v; },
    setLabelsVisible(v) { labelsWanted = v; for (const l of labels) l.visible = v; },
  };
}
