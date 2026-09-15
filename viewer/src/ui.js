// DOM overlay: playback controls, display toggles, the run info panel and
// the selected-spacecraft panel.
import { utcString, elapsedString } from 'viewer/timeline.js';
import { INFERNO_CSS } from 'viewer/colormaps.js';

const SPEEDS = [1, 10, 60, 300, 600, 3600, 21600, 86400];
const TRAIL_CHOICES = [
  { label: 'trail off', orbits: 0 },
  { label: '½ orbit', orbits: 0.5 },
  { label: '1 orbit', orbits: 1 },
  { label: '3 orbits', orbits: 3 },
  { label: '5 orbits', orbits: 5 },
  { label: '10 orbits', orbits: 10 },
  { label: 'all history', orbits: Infinity },
];

export function createUI(container, timeline, state, info) {
  const style = document.createElement('style');
  style.textContent = `
    .sa-ui { position: absolute; left: 0; right: 0; bottom: 0; padding: 10px 14px;
      font: 13px/1.4 system-ui, sans-serif; color: #e8ecf1;
      background: linear-gradient(to top, rgba(8,12,20,0.92), rgba(8,12,20,0.0)); pointer-events: none; }
    .sa-row { display: flex; gap: 10px; align-items: center; flex-wrap: wrap; pointer-events: auto; }
    .sa-ui button, .sa-ui select { font: inherit; color: #e8ecf1; background: #1d2633; border: 1px solid #3a4658;
      border-radius: 4px; padding: 3px 9px; cursor: pointer; }
    .sa-ui button.active { background: #2f5d9e; border-color: #4f86d6; }
    .sa-ui button:focus-visible, .sa-ui select:focus-visible, .sa-ui input:focus-visible { outline: 2px solid #f2b950; outline-offset: 1px; }
    .sa-ui input[type=range] { flex: 1 1 240px; min-width: 160px; }
    .sa-time { font-variant-numeric: tabular-nums; min-width: 22ch; }
    .sa-panel { position: absolute; padding: 8px 10px; max-width: 36ch;
      font: 12px/1.45 system-ui, sans-serif; color: #dbe2ea; background: rgba(8,12,20,0.72);
      border: 1px solid #2a3442; border-radius: 6px; pointer-events: auto; }
    .sa-info { left: 14px; top: 12px; }
    .sa-select { right: 14px; top: 12px; min-width: 26ch; }
    .sa-select[hidden] { display: none; }
    .sa-panel h1 { margin: 0 0 4px; font-size: 14px; font-weight: 600; color: #fff; }
    .sa-panel dl { margin: 0; display: grid; grid-template-columns: max-content 1fr; gap: 1px 8px; }
    .sa-panel dt { color: #8fa1b8; }
    .sa-panel dd { margin: 0; font-variant-numeric: tabular-nums; }
    .sa-panel .sa-hint { margin: 6px 0 0; color: #8fa1b8; }
    .sa-toggle label { margin-right: 6px; cursor: pointer; user-select: none; white-space: nowrap; }
    .sa-sep { width: 1px; height: 18px; background: #3a4658; }
    .sa-trail-legend { display: inline-flex; align-items: center; gap: 6px; font-variant-numeric: tabular-nums; }
    .sa-trail-legend .sa-bar-inferno { width: 90px; height: 9px; border-radius: 2px; background: ${INFERNO_CSS}; }
    .sa-trail-legend[hidden] { display: none; }
  `;
  document.head.appendChild(style);

  const infoBox = document.createElement('div');
  infoBox.className = 'sa-panel sa-info';
  const rows = Object.entries(info).filter(([k]) => k !== 'title').map(([k, v]) => `<dt>${k}</dt><dd>${v}</dd>`).join('');
  infoBox.innerHTML = `<h1>${info.title ?? 'SpaceAGORA run'}</h1><dl>${rows}</dl><p class="sa-hint">Click a spacecraft to inspect it.</p>`;
  container.appendChild(infoBox);

  const selectBox = document.createElement('div');
  selectBox.className = 'sa-panel sa-select';
  selectBox.hidden = true;
  container.appendChild(selectBox);

  const ui = document.createElement('div');
  ui.className = 'sa-ui';
  ui.innerHTML = `
    <div class="sa-row">
      <button data-role="play">Pause</button>
      <label>Speed <select data-role="speed"></select></label>
      <input type="range" data-role="scrub" min="0" max="1" step="0.0001" value="0">
      <span class="sa-time" data-role="time"></span>
    </div>
    <div class="sa-row sa-toggle" style="margin-top:6px">
      <span>Frame:</span>
      <button data-role="frame-inertial">Inertial</button>
      <button data-role="frame-fixed">Planet-fixed</button>
      <button data-role="follow" title="Keep the selected spacecraft centerd (F)">Follow</button>
      <span class="sa-sep"></span>
      <select data-role="trail" title="Trail length"></select>
      <select data-role="trailcolor" title="Trail color"></select>
      <span class="sa-trail-legend" data-role="trail-legend" hidden><span data-role="legend-lo"></span><span class="sa-bar-inferno"></span><span data-role="legend-hi"></span></span>
      <label><input type="checkbox" data-role="labels" checked> labels</label>
      <label data-role="paths-label"><input type="checkbox" data-role="paths" checked> planned paths</label>
      <label data-role="references-label"><input type="checkbox" data-role="references" checked> reference ghosts</label>
      <label><input type="checkbox" data-role="grid" checked> graticule</label>
      <span class="sa-sep"></span>
      <label><input type="checkbox" data-role="assemblies" checked> 3D models</label>
      <label><input type="checkbox" data-role="thrusters" checked> thrusters</label>
      <label><input type="checkbox" data-role="facets" checked> facets</label>
      <label><input type="checkbox" data-role="axes" checked> body axes</label>
      <label data-role="heating-label"><input type="checkbox" data-role="heating" checked> heating</label>
      <span class="sa-trail-legend" data-role="heat-legend" hidden><span data-role="heat-lo"></span><span class="sa-bar-inferno"></span><span data-role="heat-hi"></span></span>
      <span class="sa-sep" data-role="atmo-sep"></span>
      <label data-role="atmo-limb-label"><input type="checkbox" data-role="atmo-limb" checked> atmosphere</label>
      <label data-role="atmo-layers-label"><input type="checkbox" data-role="atmo-layers" checked> density shells</label>
      <label data-role="atmo-map-label"><input type="checkbox" data-role="atmo-map" checked> density map</label>
      <button data-role="reset">Reset view</button>
      <button data-role="video" title="Render the animation to an MP4 file">Save video…</button>
    </div>`;
  container.appendChild(ui);

  const q = (role) => ui.querySelector(`[data-role="${role}"]`);
  const play = q('play'), speed = q('speed'), scrub = q('scrub'), time = q('time');
  const frameInertial = q('frame-inertial'), frameFixed = q('frame-fixed'), follow = q('follow'), trail = q('trail');

  for (const s of SPEEDS) {
    const opt = document.createElement('option');
    opt.value = String(s);
    opt.textContent = s >= 3600 ? `${s / 3600}h/s` : s >= 60 ? `${s / 60}min/s` : `${s}x`;
    if (s === timeline.speed) opt.selected = true;
    speed.appendChild(opt);
  }
  for (const c of TRAIL_CHOICES) {
    const opt = document.createElement('option');
    opt.value = String(c.orbits);
    opt.textContent = c.label;
    if (c.orbits === state.trailOrbits) opt.selected = true;
    trail.appendChild(opt);
  }
  if (!state.trailsAvailable) trail.disabled = true;
  const trailColor = q('trailcolor');
  for (const mode of state.trailColorModes || ['age']) {
    const opt = document.createElement('option');
    opt.value = mode;
    opt.textContent = mode === 'age' ? 'color: age' : `color: ${state.trailColorLabels?.[mode] ?? mode}`;
    trailColor.appendChild(opt);
  }
  trailColor.addEventListener('change', () => state.setTrailColor(trailColor.value));
  if (!state.trailsAvailable) trailColor.disabled = true;
  const legend = q('trail-legend');
  const fmtLegend = (v) => (Math.abs(v) >= 1e4 || (Math.abs(v) < 1e-2 && v !== 0) ? v.toExponential(1) : v.toPrecision(3));
  function setTrailLegend(info) {
    legend.hidden = !info;
    if (!info) return;
    q('legend-lo').textContent = fmtLegend(info.lo);
    q('legend-hi').textContent = `${fmtLegend(info.hi)} ${info.log ? '(log)' : ''}`;
  }
  const atmoAvailable = !!state.hasAtmosphere;
  for (const role of ['atmo-sep', 'atmo-limb-label', 'atmo-layers-label', 'atmo-map-label']) q(role).hidden = !atmoAvailable;
  if (!state.hasDensityMap) q('atmo-map-label').hidden = true;
  q('atmo-limb').addEventListener('change', (e) => state.setAtmosphereLimb(e.target.checked));
  q('atmo-layers').addEventListener('change', (e) => state.setAtmosphereLayers(e.target.checked));
  q('atmo-map').addEventListener('change', (e) => state.setAtmosphereMap(e.target.checked));

  play.addEventListener('click', () => timeline.togglePlay());
  speed.addEventListener('change', () => timeline.setSpeed(Number(speed.value)));
  let scrubbing = false;
  scrub.addEventListener('input', () => { scrubbing = true; timeline.seekFraction(Number(scrub.value)); });
  scrub.addEventListener('change', () => { scrubbing = false; });
  frameInertial.addEventListener('click', () => state.setFrame('inertial'));
  frameFixed.addEventListener('click', () => state.setFrame('planet_fixed'));
  follow.addEventListener('click', () => state.setFollow(!state.follow));
  trail.addEventListener('change', () => state.setTrailOrbits(Number(trail.value)));
  q('labels').addEventListener('change', (e) => state.setLabels(e.target.checked));
  q('paths-label').hidden = !state.hasPaths;
  q('paths').addEventListener('change', (e) => state.setPaths(e.target.checked));
  q('references-label').hidden = !state.hasReferences;
  q('references').addEventListener('change', (e) => state.setReferences(e.target.checked));
  q('grid').addEventListener('change', (e) => state.setGraticule(e.target.checked));
  q('assemblies').addEventListener('change', (e) => state.setAssemblies(e.target.checked));
  q('thrusters').addEventListener('change', (e) => state.setThrusters(e.target.checked));
  q('facets').addEventListener('change', (e) => state.setFacets(e.target.checked));
  q('axes').addEventListener('change', (e) => state.setAxes(e.target.checked));
  q('heating-label').hidden = !state.hasHeating;
  q('heating').addEventListener('change', (e) => state.setHeating(e.target.checked));
  const heatLegend = q('heat-legend');
  function setHeatLegend(info) {
    heatLegend.hidden = !info;
    if (!info) return;
    q('heat-lo').textContent = `${fmtLegend(info.lo / 1e4)}`;
    q('heat-hi').textContent = `${fmtLegend(info.hi / 1e4)} W/cm² ½ρV³cosθ (log)`;
  }
  q('reset').addEventListener('click', () => state.resetView());
  q('video').addEventListener('click', () => state.openVideoDialog && state.openVideoDialog());

  window.addEventListener('keydown', (e) => {
    if (e.target && ['INPUT', 'SELECT', 'TEXTAREA'].includes(e.target.tagName)) return;
    if (e.code === 'Space') { e.preventDefault(); timeline.togglePlay(); }
    if (e.code === 'ArrowRight') timeline.seek(timeline.t + timeline.speed);
    if (e.code === 'ArrowLeft') timeline.seek(timeline.t - timeline.speed);
    if (e.code === 'KeyF') state.setFollow(!state.follow);
    if (e.code === 'Escape') state.select(-1);
  });

  function render() {
    play.textContent = timeline.playing ? 'Pause' : 'Play';
    if (!scrubbing) scrub.value = String(timeline.fraction());
    time.textContent = `${elapsedString(timeline.t - timeline.tStart)}  ${utcString(state.epochUtc, timeline.t)}`;
    frameInertial.classList.toggle('active', state.frame === 'inertial');
    frameFixed.classList.toggle('active', state.frame === 'planet_fixed');
    follow.classList.toggle('active', state.follow);
    follow.disabled = state.selected < 0 && !state.singleSpacecraft;
  }

  // Selected-spacecraft panel; `sel` is null or an object of label -> value strings.
  function setSelection(sel) {
    if (!sel) { selectBox.hidden = true; return; }
    selectBox.hidden = false;
    const body = Object.entries(sel).filter(([k]) => k !== 'title').map(([k, v]) => `<dt>${k}</dt><dd>${v}</dd>`).join('');
    selectBox.innerHTML = `<h1>${sel.title}</h1><dl>${body}</dl><p class="sa-hint">F follows · Esc deselects</p>`;
  }

  timeline.onChange(render);
  render();
  return { render, setSelection, setTrailLegend, setHeatLegend, infoBox, selectBox, root: ui };
}
