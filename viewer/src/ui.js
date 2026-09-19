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
    .sa-brand { position: absolute; right: 14px; top: 12px; display: flex; align-items: center;
      gap: 9px; height: 64px; color: #e8ecf1; font: 600 13px/1.3 system-ui, sans-serif;
      pointer-events: none; user-select: none; }
    .sa-brand img { display: block; width: 64px; height: 64px; object-fit: contain; }
    @media (max-width: 600px) {
      .sa-brand span { display: none; }
      .sa-info { max-width: calc(100% - 128px); }
    }
    .sa-select { right: 14px; top: var(--sa-panel-top, 88px); min-width: 26ch; max-height: calc(100% - var(--sa-toolbar, 100px) - var(--sa-panel-top, 88px) - 32px); overflow-y: auto; overscroll-behavior: contain; }
    .sa-select[hidden] { display: none; }
    .sa-panel h1 { margin: 0 0 4px; font-size: 14px; font-weight: 600; color: #fff; }
    .sa-panel dl { margin: 0; display: grid; grid-template-columns: max-content 1fr; gap: 1px 8px; }
    .sa-panel dt { color: #8fa1b8; }
    .sa-panel dd { margin: 0; font-variant-numeric: tabular-nums; }
    .sa-panel .sa-hint { margin: 6px 0 0; color: #8fa1b8; }
    .sa-panel h2 { margin: 6px 0 2px; font-size: 12px; font-weight: 600; color: #c9d4e0; }
    .sa-panel dt[data-key], .sa-panel dd[data-key] { cursor: pointer; }
    .sa-panel dt[data-key] { text-decoration: underline dotted #4f86d6; text-underline-offset: 2px; }
    .sa-panel dt[data-key]:hover, .sa-panel dd[data-key]:hover, .sa-panel dt.plotted { color: #f2b950; }
    .sa-face { right: 14px; top: 300px; min-width: 26ch; max-width: 40ch; }
    .sa-face[hidden] { display: none; }
    .sa-plot { left: 14px; bottom: 96px; width: 460px; max-width: calc(100% - 28px); }
    .sa-plot[hidden] { display: none; }
    .sa-plot-head { display: flex; gap: 10px; align-items: baseline; flex-wrap: wrap; }
    .sa-plot-head h1 { flex: 1 1 auto; margin: 0; }
    .sa-plot-head label { color: #8fa1b8; cursor: pointer; user-select: none; }
    .sa-plot-head button { font: inherit; color: #e8ecf1; background: #1d2633; border: 1px solid #3a4658; border-radius: 4px; padding: 0 7px; cursor: pointer; line-height: 1.4; }
    .sa-plot-value { font-variant-numeric: tabular-nums; color: #f2b950; }
    .sa-plot canvas { display: block; width: 100%; height: 220px; margin-top: 4px; }
    .sa-toggle label { margin-right: 6px; cursor: pointer; user-select: none; white-space: nowrap; }
    .sa-sep { width: 1px; height: 18px; background: #3a4658; }
    .sa-trail-legend { display: inline-flex; align-items: center; gap: 6px; font-variant-numeric: tabular-nums; }
    .sa-trail-legend .sa-bar-inferno { width: 90px; height: 9px; border-radius: 2px; background: ${INFERNO_CSS}; }
    .sa-trail-legend[hidden] { display: none; }
  `;
  document.head.appendChild(style);

  const brand = document.createElement('div');
  brand.className = 'sa-brand';
  const logo = document.createElement('img');
  logo.src = SPACE_FALCON_LAB_LOGO;
  logo.alt = 'Space-FALCON Lab';
  logo.width = 64;
  logo.height = 64;
  logo.draggable = false;
  const labName = document.createElement('span');
  labName.textContent = 'Space-FALCON Lab';
  labName.setAttribute('aria-hidden', 'true');
  brand.append(logo, labName);
  container.appendChild(brand);

  const infoBox = document.createElement('div');
  infoBox.className = 'sa-panel sa-info';
  // Labels can come from user titles, model names and saved scene metadata.
  // Keep them out of innerHTML even though the outer JSON is script-safe.
  infoBox.innerHTML = '<h1></h1><dl></dl><p class="sa-hint">Click a spacecraft to inspect it.</p>';
  infoBox.querySelector('h1').textContent = info.title ?? 'SpaceAGORA run';
  const rows = infoBox.querySelector('dl');
  for (const [key, value] of Object.entries(info)) {
    if (key === 'title') continue;
    const label = document.createElement('dt');
    label.textContent = key;
    const entry = document.createElement('dd');
    entry.dataset.info = key;
    entry.textContent = String(value);
    rows.append(label, entry);
  }
  container.appendChild(infoBox);

  const selectBox = document.createElement('div');
  selectBox.className = 'sa-panel sa-select';
  selectBox.hidden = true;
  container.appendChild(selectBox);

  // Picked-face panel, placed under the selection panel.
  const faceBox = document.createElement('div');
  faceBox.className = 'sa-panel sa-face';
  faceBox.hidden = true;
  container.appendChild(faceBox);
  // A click on a quantity (a row carrying data-key) opens its time history.
  for (const box of [selectBox, faceBox]) {
    box.addEventListener('click', (e) => {
      const el = e.target.closest('[data-key]');
      if (el && state.plot) state.plot(el.dataset.key);
    });
  }

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
      <label><input type="checkbox" data-role="groundtracks"> ground tracks</label>
      <label><input type="checkbox" data-role="grid" checked> graticule</label>
      <span class="sa-sep"></span>
      <label><input type="checkbox" data-role="assemblies" checked> 3D models</label>
      <label><input type="checkbox" data-role="thrusters" checked> thrusters</label>
      <label data-role="plumes-label"><input type="checkbox" data-role="plumes" checked> plumes</label>
      <label><input type="checkbox" data-role="facets" checked> facets</label>
      <label><input type="checkbox" data-role="axes" checked> body axes</label>
      <label data-role="heating-label"><input type="checkbox" data-role="heating" checked> heating</label>
      <select data-role="lighting" title="Sun lighting and renderer"></select>
      <span class="sa-trail-legend" data-role="heat-legend" hidden><span data-role="heat-lo"></span><span class="sa-bar-inferno"></span><span data-role="heat-hi"></span></span>
      <label data-role="dust-label"><input type="checkbox" data-role="dust" checked> dust</label>
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
  const groundTracks = q('groundtracks');
  groundTracks.checked = !!state.groundTracks;
  groundTracks.addEventListener('change', (e) => state.setGroundTracks(e.target.checked));
  q('grid').addEventListener('change', (e) => state.setGraticule(e.target.checked));
  q('assemblies').addEventListener('change', (e) => state.setAssemblies(e.target.checked));
  q('thrusters').addEventListener('change', (e) => state.setThrusters(e.target.checked));
  q('plumes-label').hidden = !state.hasPlumes;
  q('plumes').addEventListener('change', (e) => state.setPlumes(e.target.checked));
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
  q('dust-label').hidden = !state.hasDust;
  q('dust').addEventListener('change', (e) => state.setDust(e.target.checked));
  // Lighting: real-time, plus "path traced when paused" once the path tracer
  // has been imported (main.js refreshes the list when it knows).
  const lightingSelect = q('lighting');
  function setLightingModes(modes) {
    const list = modes && modes.length ? modes : [{ value: 'realtime', label: 'lighting: real-time' }];
    const current = state.lightingMode || list[0].value;
    lightingSelect.innerHTML = '';
    for (const m of list) {
      const opt = document.createElement('option');
      opt.value = m.value;
      opt.textContent = m.label;
      if (m.value === current) opt.selected = true;
      lightingSelect.appendChild(opt);
    }
  }
  setLightingModes(state.lightingModes);
  lightingSelect.addEventListener('change', () => state.setLighting && state.setLighting(lightingSelect.value));
  // The lighting row of the run panel, refreshed with the accumulated samples.
  const lightingStatus = infoBox.querySelector('[data-info="lighting"]');
  function setLightingStatus(text) {
    if (lightingStatus && text != null && lightingStatus.textContent !== text) lightingStatus.textContent = text;
  }

  q('reset').addEventListener('click', () => state.resetView());
  q('video').addEventListener('click', () => state.openVideoDialog && state.openVideoDialog());

  window.addEventListener('keydown', (e) => {
    if (e.target && ['INPUT', 'SELECT', 'TEXTAREA'].includes(e.target.tagName)) return;
    if (e.code === 'Space') { e.preventDefault(); timeline.togglePlay(); }
    if (e.code === 'ArrowRight') timeline.seek(timeline.t + timeline.speed);
    if (e.code === 'ArrowLeft') timeline.seek(timeline.t - timeline.speed);
    if (e.code === 'KeyF') state.setFollow(!state.follow);
    if (e.code === 'Escape') { if (state.face) state.setFace(null); else state.select(-1); }
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

  // Panel body from rows [{ label, text, key }]: rows with a key are clickable
  // quantities. The DOM is rebuilt only when the set of rows changes; between
  // rebuilds only the values are written, so a click lands on a stable element.
  const escapeHtml = (v) => String(v).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;').replace(/"/g, '&quot;');
  function renderRows(box, model, hint) {
    const sig = `${model.title}|${model.rows.map((r) => `${r.label}:${r.key || ''}`).join(';')}|${hint}`;
    if (box.dataset.sig !== sig) {
      box.dataset.sig = sig;
      const body = model.rows.map((r) => {
        const attr = r.key ? ` data-key="${escapeHtml(r.key)}" title="Click for the time history"` : '';
        return `<dt${attr}>${escapeHtml(r.label)}</dt><dd${attr}></dd>`;
      }).join('');
      box.innerHTML = `<h1>${escapeHtml(model.title)}</h1><dl>${body}</dl>${hint ? `<p class="sa-hint">${hint}</p>` : ''}`;
      box._dds = Array.from(box.querySelectorAll('dd'));
      box._dts = Array.from(box.querySelectorAll('dt'));
    }
    model.rows.forEach((r, i) => {
      const dd = box._dds[i];
      if (dd && dd.textContent !== r.text) dd.textContent = r.text;
      const dt = box._dts[i];
      if (dt) dt.classList.toggle('plotted', !!r.key && r.key === plottedKey);
    });
  }
  let plottedKey = null;
  function setPlotted(key) { plottedKey = key; }

  // Selected-spacecraft panel; `sel` is null or { title, rows: [{ label, text, key }] }.
  function setSelection(sel) {
    if (!sel) { selectBox.hidden = true; return; }
    selectBox.hidden = false;
    renderRows(selectBox, sel, 'Click a value for its history · click the spacecraft body for a face · F follows · Esc deselects');
  }
  // Picked-face panel; `face` is null or { title, rows }.
  function setFace(face) {
    if (!face) { faceBox.hidden = true; return; }
    faceBox.hidden = false;
    renderRows(faceBox, face, 'Esc clears the face');
    faceBox.style.top = `${selectBox.hidden ? parseFloat(getComputedStyle(selectBox).top) : selectBox.offsetTop + selectBox.offsetHeight + 8}px`;
  }

  timeline.onChange(render);
  render();
  return { render, setSelection, setFace, setPlotted, setTrailLegend, setHeatLegend, setLightingModes, setLightingStatus, infoBox, selectBox, faceBox, root: ui };
}

// Existing Space-FALCON Lab blue signature artwork, embedded unchanged for offline pages.
// PNG SHA-256: 8e961cfeb03c73671acba3c4ba420d99e8598dcfb8cabbeb8ec1e90feae7b79d
const SPACE_FALCON_LAB_LOGO = "data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAEsAAABKCAYAAADzEqlPAAAABGdBTUEAALGPC/xhBQAAACBjSFJNAAB6JgAAgIQAAPoAAACA6AAAdTAAAOpgAAA6mAAAF3CculE8AAAApmVYSWZNTQAqAAAACAAFARIAAwAAAAEAAQAAARoABQAAAAEAAABKARsABQAAAAEAAABSATEAAgAAACEAAABah2kABAAAAAEAAAB8AAAAAAAAASwAAAABAAABLAAAAAFBZG9iZSBQaG90b3Nob3AgMjUuMCAoTWFjaW50b3NoKQAAAAOgAQADAAAAAQABAACgAgAEAAAAAQAAAEugAwAEAAAAAQAAAEoAAAAA+bt6JQAAAAlwSFlzAAAuIwAALiMBeKU/dgAABlppVFh0WE1MOmNvbS5hZG9iZS54bXAAAAAAADx4OnhtcG1ldGEgeG1sbnM6eD0iYWRvYmU6bnM6bWV0YS8iIHg6eG1wdGs9IlhNUCBDb3JlIDYuMC4wIj4KICAgPHJkZjpSREYgeG1sbnM6cmRmPSJodHRwOi8vd3d3LnczLm9yZy8xOTk5LzAyLzIyLXJkZi1zeW50YXgtbnMjIj4KICAgICAgPHJkZjpEZXNjcmlwdGlvbiByZGY6YWJvdXQ9IiIKICAgICAgICAgICAgeG1sbnM6ZXhpZj0iaHR0cDovL25zLmFkb2JlLmNvbS9leGlmLzEuMC8iCiAgICAgICAgICAgIHhtbG5zOnhtcE1NPSJodHRwOi8vbnMuYWRvYmUuY29tL3hhcC8xLjAvbW0vIgogICAgICAgICAgICB4bWxuczpzdFJlZj0iaHR0cDovL25zLmFkb2JlLmNvbS94YXAvMS4wL3NUeXBlL1Jlc291cmNlUmVmIyIKICAgICAgICAgICAgeG1sbnM6eG1wPSJodHRwOi8vbnMuYWRvYmUuY29tL3hhcC8xLjAvIgogICAgICAgICAgICB4bWxuczpwaG90b3Nob3A9Imh0dHA6Ly9ucy5hZG9iZS5jb20vcGhvdG9zaG9wLzEuMC8iCiAgICAgICAgICAgIHhtbG5zOnRpZmY9Imh0dHA6Ly9ucy5hZG9iZS5jb20vdGlmZi8xLjAvIj4KICAgICAgICAgPGV4aWY6Q29sb3JTcGFjZT4xPC9leGlmOkNvbG9yU3BhY2U+CiAgICAgICAgIDxleGlmOlBpeGVsWERpbWVuc2lvbj4xMDA8L2V4aWY6UGl4ZWxYRGltZW5zaW9uPgogICAgICAgICA8ZXhpZjpQaXhlbFlEaW1lbnNpb24+OTg8L2V4aWY6UGl4ZWxZRGltZW5zaW9uPgogICAgICAgICA8eG1wTU06T3JpZ2luYWxEb2N1bWVudElEPnhtcC5kaWQ6ZDExZDM2OTctZjg4OC00NzUwLWIzZDgtNTcwOGM2MzZiN2ZiPC94bXBNTTpPcmlnaW5hbERvY3VtZW50SUQ+CiAgICAgICAgIDx4bXBNTTpEb2N1bWVudElEPnhtcC5kaWQ6RTlBNjg0RUY0Q0I5Q0Q0RkFDRUYxOTYxOTY4MUMzMEY8L3htcE1NOkRvY3VtZW50SUQ+CiAgICAgICAgIDx4bXBNTTpEZXJpdmVkRnJvbSByZGY6cGFyc2VUeXBlPSJSZXNvdXJjZSI+CiAgICAgICAgICAgIDxzdFJlZjppbnN0YW5jZUlEPnhtcC5paWQ6MzhkNmFiOGMtMzVhNi00NTkzLTgyOTUtMDFkMzY1YjA2NGI4PC9zdFJlZjppbnN0YW5jZUlEPgogICAgICAgICAgICA8c3RSZWY6ZG9jdW1lbnRJRD5hZG9iZTpkb2NpZDpwaG90b3Nob3A6MmZjMmM3OGYtNmQzYi0wMDQyLTlkZjEtM2ZhZWM4NmM2ZmUxPC9zdFJlZjpkb2N1bWVudElEPgogICAgICAgICA8L3htcE1NOkRlcml2ZWRGcm9tPgogICAgICAgICA8eG1wTU06SW5zdGFuY2VJRD54bXAuaWlkOkFFNDBBMzQzMTAxRUJENDFBNzE3MjIyMEVFRkIzNjdDPC94bXBNTTpJbnN0YW5jZUlEPgogICAgICAgICA8eG1wOkNyZWF0b3JUb29sPkFkb2JlIFBob3Rvc2hvcCAyNS4wIChNYWNpbnRvc2gpPC94bXA6Q3JlYXRvclRvb2w+CiAgICAgICAgIDxwaG90b3Nob3A6SUNDUHJvZmlsZT5zUkdCIElFQzYxOTY2LTIuMTwvcGhvdG9zaG9wOklDQ1Byb2ZpbGU+CiAgICAgICAgIDx0aWZmOk9yaWVudGF0aW9uPjE8L3RpZmY6T3JpZW50YXRpb24+CiAgICAgICAgIDx0aWZmOlhSZXNvbHV0aW9u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";
