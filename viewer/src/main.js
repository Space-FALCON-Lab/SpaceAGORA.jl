// Entry point. Expects `window.SPACEAGORA_VIEWER` (embedded by the Julia
// bundler, or loaded from dev_data.js in viewer/dev.html):
//   { scene, frames, textures, models, options }
import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';
import { FrameData, geodetic, rotateByConjugate } from 'viewer/data.js';
import { createGlobe } from 'viewer/globe.js';
import { createSpacecraft, estimateOrbitPeriod } from 'viewer/spacecraft.js';
import { createAssemblies } from 'viewer/lod.js';
import { createEnsemble } from 'viewer/ensemble.js';
import { createAtmosphere } from 'viewer/atmosphere.js';
import { createPaths } from 'viewer/paths.js';
import { createReferences } from 'viewer/references.js';
import { TRAIL_COLOR_MODES } from 'viewer/spacecraft.js';
import { Timeline } from 'viewer/timeline.js';
import { createUI } from 'viewer/ui.js';
import { createVideoDialog } from 'viewer/video.js';

const DEFAULT_TRAIL_ORBITS = 3;

export function start(payload, container = document.body) {
  const { scene: sidecar, frames: rawFrames, textures = {}, models = {}, options = {}, ensemble: ensembleSpec = null, paths: pathSpecs = [], references: referenceSpecs = [] } = payload;
  const frames = new FrameData(rawFrames);
  const planet = sidecar.planet;
  const Re = planet.equatorial_radius_m / 1000, Rp = planet.polar_radius_m / 1000;

  // Renderer. Scene units are kilometers; the logarithmic depth buffer keeps a
  // 6000 km globe and a 3 m spacecraft box in one view.
  const renderer = new THREE.WebGLRenderer({ antialias: true, logarithmicDepthBuffer: true });
  const pixelRatio = Math.min(window.devicePixelRatio || 1, 2);
  renderer.setPixelRatio(pixelRatio);
  renderer.setSize(container.clientWidth || window.innerWidth, container.clientHeight || window.innerHeight);
  renderer.setClearColor(0x05070c, 1);
  container.style.position = container.style.position || 'relative';
  container.appendChild(renderer.domElement);

  const scene = new THREE.Scene();
  const camera = new THREE.PerspectiveCamera(45, 1, 1e-5, 1e8);
  const controls = new OrbitControls(camera, renderer.domElement);
  controls.enableDamping = true;
  controls.dampingFactor = 0.08;
  controls.minDistance = 0;
  controls.maxDistance = 5e6;

  scene.add(new THREE.AmbientLight(0xffffff, 0.55));
  const sun = new THREE.DirectionalLight(0xffffff, 1.6);
  sun.position.set(1, 0.3, 0.2).multiplyScalar(1e6);
  scene.add(sun);

  // `world` holds everything expressed in the inertial frame. In planet-fixed
  // mode it is counter-rotated by the globe's rotation so the body stands
  // still; in follow mode it is translated so the followed spacecraft sits at
  // the origin (a floating origin keeps meter-scale geometry steady at
  // planetary distances).
  const world = new THREE.Group();
  scene.add(world);

  const textureKey = (planet.texture || planet.name || '').toLowerCase();
  const globe = createGlobe(planet, textures[textureKey] || null, {
    anisotropy: renderer.capabilities.getMaxAnisotropy(),
    maxTextureSize: renderer.capabilities.maxTextureSize,
  });
  world.add(globe.group);

  const craft = createSpacecraft(frames, sidecar, { markerPixels: options.marker_pixels ?? 7, pixelRatio, planetRadiusKm: Re });
  world.add(craft.group);

  // Atmosphere: limb glow, density shells and map, rotating with the body.
  const atmosphere = sidecar.atmosphere ? createAtmosphere(sidecar.atmosphere, planet, {}) : null;
  if (atmosphere) globe.group.add(atmosphere.group);

  const lod = createAssemblies(sidecar, frames, { models, enabled: options.assemblies ?? true, assemblyLimit: options.assembly_limit ?? 256 });
  world.add(lod.group);

  const refPaths = createPaths(pathSpecs, frames, {});
  if (refPaths.items.length) world.add(refPaths.group);

  // Reference ghosts (SPICE reconstructions, telemetry, plans) beside the flown spacecraft.
  const refs = createReferences(referenceSpecs, sidecar, frames, models, {});
  if (refs.items.length) world.add(refs.group);

  // The ensemble panel needs `state.select` before `state` is built; it reads
  // through this reference, which is filled in below.
  const stateRef = { selected: -1, select: (i) => state.select(i) };

  const ensemble = ensembleSpec ? createEnsemble(ensembleSpec, frames, craft, container, stateRef, {}) : null;
  if (ensemble) world.add(ensemble.group);

  const period = estimateOrbitPeriod(frames, 0);
  const span = frames.tEnd - frames.tStart;
  const timeline = new Timeline(frames.tStart, frames.tEnd, { speed: options.speed ?? defaultSpeed(frames), autoplay: true });

  const state = {
    frame: options.frame === 'planet_fixed' ? 'planet_fixed' : 'inertial',
    follow: false,
    selected: -1,
    trailOrbits: options.trail_orbits ?? DEFAULT_TRAIL_ORBITS,
    trailsAvailable: craft.trailsEnabled,
    singleSpacecraft: frames.sats === 1,
    trailColorModes: craft.availableColorModes(),
    trailColorLabels: Object.fromEntries(Object.entries(TRAIL_COLOR_MODES).map(([k, v]) => [k, v.label])),
    hasAtmosphere: !!atmosphere,
    hasPaths: refPaths.items.length > 0,
    setPaths(v) { refPaths.setVisible(v); },
    hasReferences: refs.items.length > 0,
    setReferences(v) { refs.setVisible(v); },
    hasDensityMap: !!(atmosphere && atmosphere.map),
    setTrailColor(mode) { ui.setTrailLegend(craft.setTrailColorMode(mode)); },
    setAtmosphereLimb(v) { atmosphere && atmosphere.setLimbVisible(v); },
    setAtmosphereLayers(v) { atmosphere && atmosphere.setLayersVisible(v); },
    setAtmosphereMap(v) { atmosphere && atmosphere.setMapVisible(v); },
    epochUtc: sidecar.epoch.utc,
    setFrame(f) { state.frame = f; ui.render(); },
    setFollow(v) {
      // A single-spacecraft run has an obvious target; otherwise a selection is needed.
      if (v && state.selected < 0) {
        if (frames.sats === 1) state.selected = 0; else return;
      }
      state.follow = v;
      if (v) approachSelected(); else { world.position.set(0, 0, 0); }
      ui.render();
    },
    select(i) {
      state.selected = i;
      stateRef.selected = i;
      if (i < 0 && state.follow) state.setFollow(false);
      ui.render();
    },
    setTrailOrbits(n) { state.trailOrbits = n; applyTrail(); },
    setLabels(v) { craft.setLabelsVisible(v); refs.setLabelsVisible(v); },
    setGraticule(v) { globe.grid.visible = v; globe.axis.visible = v; },
    setAssemblies(v) { lod.setEnabled(v); },
    setThrusters(v) { lod.setThrustersVisible(v); },
    setFacets(v) { lod.setFacetsVisible(v); },
    setAxes(v) { lod.setAxesVisible(v); },
    hasHeating: lod.heatingAvailable,
    setHeating(v) { lod.setHeatingVisible(v); ui.setHeatLegend(v && lod.heatingAvailable ? lod.heatRange() : null); },
    resetView() { state.setFollow(false); placeCamera(); },
    openVideoDialog() { videoDialog && videoDialog.open(); },
  };
  let videoDialog = null;

  // An explicit trail_s option wins until the user picks an orbit count.
  let trailTouched = false;
  function applyTrail() {
    if (options.trail_s != null && !trailTouched) craft.setTrailSeconds(options.trail_s);
    else craft.setTrailSeconds(state.trailOrbits === Infinity ? span : state.trailOrbits * period);
  }
  applyTrail();
  const setTrailOrbits = state.setTrailOrbits;
  state.setTrailOrbits = (n) => { trailTouched = true; setTrailOrbits(n); };

  const cadence = frames.count > 1 ? span / (frames.count - 1) : 0;
  const ui = createUI(container, timeline, state, {
    title: options.title ?? `SpaceAGORA · ${planet.name}`,
    body: `${planet.name}${textures[textureKey] ? ` (${textures[textureKey].resolution || 'texture'}, GPU max ${renderer.capabilities.maxTextureSize}px)` : ' (no texture)'}`,
    epoch: sidecar.epoch.utc,
    spacecraft: String(frames.sats),
    frames: `${frames.count} of ${frames.sourceRows} rows (every ${frames.strideRows})`,
    cadence: cadence > 0 ? `${cadence.toFixed(1)} s` : 'n/a',
    'orbit period': period > 0 && period < span ? `${(period / 60).toFixed(1)} min` : 'longer than run',
    attitude: frames.q ? 'saved quaternion' : 'velocity-aligned',
    'link poses': frames.linkPose ? 'recorded' : 'configured',
    ...(ensemble ? { samples: `${ensembleSpec.count} × ${ensemble.perSample} spacecraft` } : {}),
    ...(refPaths.items.length ? { paths: refPaths.items.map((it) => it.spec.name).join(', ') } : {}),
    ...(refs.items.length ? { references: refs.items.map((it) => it.spec.name).join(', ') } : {}),
    ...(atmosphere ? {
      atmosphere: `${atmosphere.info.model.replace('AtmosphereModel', '')}, EI ${atmosphere.info.ei_km.toFixed(0)} km`,
      ...(atmosphere.info.map ? { 'density map': `${atmosphere.info.map.altitude_km.toFixed(0)} km, ${atmosphere.info.map.min.toExponential(1)}–${atmosphere.info.map.max.toExponential(1)} kg/m³` } : {}),
    } : {}),
  });
  // Default trail coloring: heat rate when the run has it, else age.
  const initialColor = options.trail_color && TRAIL_COLOR_MODES[options.trail_color] ? options.trail_color : (frames.hasScalar('heat_rate') ? 'heat_rate' : 'age');
  ui.setTrailLegend(craft.setTrailColorMode(initialColor));
  if (lod.heatingAvailable && (options.heating ?? true)) { lod.setHeatingVisible(true); ui.setHeatLegend(lod.heatRange()); }
  const trailColorSelect = container.querySelector('[data-role="trailcolor"]');
  if (trailColorSelect) trailColorSelect.value = craft.trailColorMode;

  const tmpPos = new Float64Array(3), tmpVel = new Float64Array(3), qPi = new Float32Array(4), rBody = new Float64Array(3);
  const v3 = new THREE.Vector3();

  function placeCamera() {
    let far = Re;
    const p = new Float32Array(frames.sats * 3);
    frames.positionsAt(frames.tStart, p);
    for (let s = 0; s < frames.sats; s++) far = Math.max(far, Math.hypot(p[3 * s], p[3 * s + 1], p[3 * s + 2]));
    camera.position.set(-1.2 * far, -2.6 * far, 1.3 * far);
    controls.target.set(0, 0, 0);
    camera.lookAt(controls.target);
    controls.update();
  }

  function approachSelected() {
    const item = lod.items[state.selected];
    // Ten bounding radii: the assembly fills roughly a fifth of the view height.
    const dist = Math.max(0.02, 10 * item.radiusKm);
    const dir = camera.position.clone().sub(controls.target);
    if (dir.lengthSq() === 0) dir.set(-1, -1, 0.6);
    dir.normalize();
    controls.target.set(0, 0, 0);
    camera.position.copy(dir.multiplyScalar(dist));
    controls.update();
  }

  // Picking: nearest marker (or visible assembly center) within a few pixels.
  let downX = 0, downY = 0, downT = 0;
  renderer.domElement.addEventListener('pointerdown', (e) => { downX = e.clientX; downY = e.clientY; downT = performance.now(); });
  renderer.domElement.addEventListener('pointerup', (e) => {
    if (Math.hypot(e.clientX - downX, e.clientY - downY) > 5 || performance.now() - downT > 400) return;
    const rect = renderer.domElement.getBoundingClientRect();
    const px = e.clientX - rect.left, py = e.clientY - rect.top;
    let best = -1, bestScore = Infinity;
    world.updateMatrixWorld();
    for (let s = 0; s < frames.sats; s++) {
      if (lod.markerHidden[s] && !lod.visible[s]) continue; // absent at this time
      v3.set(craft.positions[3 * s], craft.positions[3 * s + 1], craft.positions[3 * s + 2]).applyMatrix4(craft.group.matrixWorld).project(camera);
      if (v3.z > 1) continue;
      const sx = (v3.x + 1) / 2 * rect.width, sy = (1 - v3.y) / 2 * rect.height;
      const d = Math.hypot(sx - px, sy - py);
      const tolerance = Math.max(14, lod.pxSize[s]);
      if (d <= tolerance && d / tolerance < bestScore) { best = s; bestScore = d / tolerance; }
    }
    state.select(best);
  });

  function selectionInfo(t) {
    const s = state.selected;
    if (s < 0) return null;
    const spec = sidecar.spacecraft[s] || sidecar.spacecraft[0];
    frames.positionAt(t, s, tmpPos);
    if (!Number.isFinite(tmpPos[0])) {
      return { title: spec.name || `sc${s + 1}`, id: String(spec.id), state: 'not present at this time' };
    }
    frames.velocityAt(t, s, tmpVel);
    globe.rotationAt(t, qPi);
    rotateByConjugate(qPi, tmpPos, rBody);
    const g = geodetic(rBody, Re, Rp);
    const mass = frames.massAt(t, s);
    const out = {
      title: spec.name || `sc${s + 1}`,
      id: String(spec.id),
      altitude: `${g.alt.toFixed(1)} km`,
      latitude: `${g.latDeg.toFixed(3)}°`,
      longitude: `${g.lonDeg.toFixed(3)}°`,
      radius: `${Math.hypot(tmpPos[0], tmpPos[1], tmpPos[2]).toFixed(1)} km`,
      speed: `${Math.hypot(tmpVel[0], tmpVel[1], tmpVel[2]).toFixed(3)} km/s`,
    };
    if (mass != null) out.mass = `${mass.toFixed(1)} kg`;
    if (frames.hasScalar('density')) {
      const rho = frames.scalarAt('density', t, s, Re);
      out.density = rho > 0 ? `${rho.toExponential(2)} kg/m³` : '0 (above atmosphere)';
      if (frames.hasScalar('dynamic_pressure')) out['dyn. pressure'] = `${frames.scalarAt('dynamic_pressure', t, s, Re).toFixed(3)} Pa`;
    }
    if (frames.hasScalar('heat_rate')) out['heat rate'] = `${(frames.scalarAt('heat_rate', t, s, Re) / 1e4).toFixed(4)} W/cm²`;
    if (frames.hasScalar('drag')) out.drag = `${frames.scalarAt('drag', t, s, Re).toFixed(3)} N`;
    if (frames.hasScalar('wind')) out.wind = `${frames.scalarAt('wind', t, s, Re).toFixed(1)} m/s`;
    out.links = String(spec.links.length);
    out.model = lod.visible[s] ? `3D (${lod.pxSize[s].toFixed(0)} px)` : `marker (${lod.pxSize[s].toFixed(1)} px)`;
    const status = lod.modelStatus(s);
    if (status) out['3D model'] = status;
    refs.items.forEach((it, k) => {
      if (it.target !== s) return;
      const sep = refs.separationKm(t, k);
      out[`vs ${it.spec.name}`] = Number.isFinite(sep) ? (sep >= 10 ? `${sep.toFixed(1)} km` : `${(1000 * sep).toFixed(1)} m`) : 'not covered';
      const ghost = refs.modelStatus(k);
      if (ghost && ghost !== 'boxes') out[`${it.spec.name} model`] = ghost;
    });
    return out;
  }

  function resize() {
    if (recording) return;
    const w = container.clientWidth || window.innerWidth, h = container.clientHeight || window.innerHeight;
    resizeTo(w, h);
  }
  function resizeTo(w, h) {
    camera.aspect = w / h;
    camera.updateProjectionMatrix();
    renderer.setSize(w, h, !recording);
  }
  let recording = false;
  window.addEventListener('resize', resize);
  resize();
  placeCamera();

  const qWorld = new THREE.Quaternion();
  const followPos = new THREE.Vector3();
  // Floating origin: the followed spacecraft's position (km, Float64). Every
  // buffer is uploaded relative to it and the world group is shifted by its
  // negative, so the GPU never subtracts two 7000 km numbers to place a 3 m box.
  const anchor = new Float64Array(3);
  let last = performance.now();
  let lastInfo = 0;
  function animate(now) {
    requestAnimationFrame(animate);
    const dt = Math.min(0.25, (now - last) / 1000);
    last = now;
    if (recording) return; // the recorder drives frame(t) itself
    timeline.tick(dt);
    frame(timeline.t);
    if (now - lastInfo > 100) { ui.setSelection(selectionInfo(timeline.t)); lastInfo = now; }
  }
  // One rendered frame at elapsed time t: every scene update and the draw.
  function frame(t) {
    globe.update(t);
    if (state.frame === 'planet_fixed') {
      globe.rotationAt(t, qWorld);
      world.quaternion.copy(qWorld).invert();
    } else {
      world.quaternion.identity();
    }
    if (state.follow && state.selected >= 0) {
      frames.positionAt(t, state.selected, tmpPos);
      if (Number.isFinite(tmpPos[0])) { anchor[0] = tmpPos[0]; anchor[1] = tmpPos[1]; anchor[2] = tmpPos[2]; }
      followPos.set(anchor[0], anchor[1], anchor[2]).applyQuaternion(world.quaternion).negate();
      world.position.copy(followPos);
      controls.target.set(0, 0, 0);
    } else {
      anchor[0] = 0; anchor[1] = 0; anchor[2] = 0;
    }
    craft.group.position.set(anchor[0], anchor[1], anchor[2]);
    lod.group.position.set(anchor[0], anchor[1], anchor[2]);
    refs.group.position.set(anchor[0], anchor[1], anchor[2]);
    if (ensemble) ensemble.group.position.set(anchor[0], anchor[1], anchor[2]);
    world.updateMatrixWorld();
    const viewportHeight = renderer.domElement.clientHeight || window.innerHeight;
    lod.update(t, camera, viewportHeight, lod.group.matrixWorld, anchor);
    if (refPaths.items.length) refPaths.update(t, anchor);
    if (refs.items.length) refs.update(t, camera, viewportHeight, refs.group.matrixWorld, anchor);
    if (ensemble) ensemble.update(state.follow, t, anchor, camera, ensemble.group.matrixWorld);
    craft.update(t, camera, lod.markerHidden, state.selected, craft.group.matrixWorld, anchor, ensemble ? ensemble.dimMask : null);
    controls.update();
    renderer.render(scene, camera);
  }
  requestAnimationFrame(animate);

  const viewer = {
    renderer, scene, camera, controls, world, globe, atmosphere, craft, lod, ensemble, references: refs, paths: refPaths, timeline, frames, state, period,
    // Deterministic rendering for exports: seek and draw one frame at t.
    renderAt(t) { timeline.seek(t); frame(t); },
    setRecording(v) { recording = v; if (!v) resize(); },
    resizeTo,
  };
  videoDialog = createVideoDialog(viewer, container);
  return viewer;
}

function defaultSpeed(frames) {
  const span = frames.tEnd - frames.tStart;
  if (span <= 0) return 1;
  // Aim for a ~60 s loop at default speed, snapped to the UI's speed ladder.
  const want = span / 60;
  const ladder = [1, 10, 60, 300, 600, 3600, 21600, 86400];
  return ladder.reduce((best, s) => (Math.abs(Math.log(s / want)) < Math.abs(Math.log(best / want)) ? s : best), ladder[0]);
}
