// Entry point. Expects `window.SPACEAGORA_VIEWER` (embedded by the Julia
// bundler, or loaded from dev_data.js in viewer/dev.html):
//   { scene, frames, textures, models, options }
import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';
import { FrameData, geodetic, rotateByConjugate } from 'viewer/data.js';
import { createGlobe } from 'viewer/globe.js';
import { createSpacecraft, estimateOrbitPeriod } from 'viewer/spacecraft.js';
import { createAssemblies } from 'viewer/lod.js';
import { createPlumes } from 'viewer/plumes.js';
import { createEnsemble } from 'viewer/ensemble.js';
import { createAtmosphere } from 'viewer/atmosphere.js';
import { createPaths } from 'viewer/paths.js';
import { createReferences } from 'viewer/references.js';
import { createGroundTracks } from 'viewer/groundtrack.js';
import { TRAIL_COLOR_MODES } from 'viewer/spacecraft.js';
import { Timeline } from 'viewer/timeline.js';
import { createUI } from 'viewer/ui.js';
import { createVideoDialog } from 'viewer/video.js';
import { createPlotPanel } from 'viewer/plots.js';
import { createTerrain } from 'viewer/terrain.js';
import { createDust } from 'viewer/dust.js';
import { createLighting } from 'viewer/lighting.js';

const DEFAULT_TRAIL_ORBITS = 3;

export function start(payload, container = document.body) {
  const { scene: sidecar, frames: rawFrames, textures = {}, models = {}, options = {}, ensemble: ensembleSpec = null, paths: pathSpecs = [], references: referenceSpecs = [], terrain: terrainSpec = null } = payload;
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

  // Lighting (viewer/lighting.js) is created below, once the bodies it has to
  // light and shadow exist.

  // `world` holds everything expressed in the inertial frame. In planet-fixed
  // mode it is counter-rotated by the globe's rotation so the body stands
  // still; in follow mode it is translated so the followed spacecraft sits at
  // the origin (a floating origin keeps meter-scale geometry steady at
  // planetary distances).
  const world = new THREE.Group();
  scene.add(world);

  const textureKey = (planet.texture || planet.name || '').toLowerCase();
  // Site terrain (a view-dependent quadtree of DEM patches with draped imagery)
  // sits in the globe group and cuts a hole in the sphere under the region it covers.
  const terrain = createTerrain(terrainSpec, planet, {
    anisotropy: renderer.capabilities.getMaxAnisotropy(),
    globeLonLeft: textures[textureKey] ? textures[textureKey].lon_left_deg : -180,
    splitPixels: options.terrain_split_pixels,
  });
  const globe = createGlobe(planet, textures[textureKey] || null, {
    anisotropy: renderer.capabilities.getMaxAnisotropy(),
    maxTextureSize: renderer.capabilities.maxTextureSize,
    hole: terrain.hole,
  });
  world.add(globe.group);
  if (terrain.levels.length) {
    globe.group.add(terrain.group);
    // The outer ring of the covered region blends into the globe's own map, so
    // the boundary of the terrain is not an edge.
    terrain.setGlobeTexture(globe.mesh.material.map || null);
  }

  const craft = createSpacecraft(frames, sidecar, { markerPixels: options.marker_pixels ?? 7, pixelRatio, planetRadiusKm: Re });
  world.add(craft.group);

  // Atmosphere: limb glow, density shells and map, rotating with the body.
  const atmosphere = sidecar.atmosphere ? createAtmosphere(sidecar.atmosphere, planet, {}) : null;
  if (atmosphere) globe.group.add(atmosphere.group);

  const lod = createAssemblies(sidecar, frames, { models, enabled: options.assemblies ?? true, assemblyLimit: options.assembly_limit ?? 256 });
  world.add(lod.group);

  // Regolith blown off the surface by a landing vehicle's descent engine; it
  // adds its own group to `world`, so it rides the inertial frame like the rest.
  const dust = createDust(world, frames, terrain, lod, { rotationAt: (t, out) => globe.rotationAt(t, out) });
  // Thruster plumes: one per thruster glyph, driven by the recorded firing levels.
  const plumes = createPlumes(sidecar, frames, lod, { raw: rawFrames, enabled: options.plumes ?? true });

  const refPaths = createPaths(pathSpecs, frames, {});
  if (refPaths.items.length) world.add(refPaths.group);

  // Reference ghosts (SPICE reconstructions, telemetry, plans) beside the flown spacecraft.
  const refs = createReferences(referenceSpecs, sidecar, frames, models, {});
  if (refs.items.length) world.add(refs.group);

  // Sub-satellite tracks on the body's surface. They live in the globe group,
  // so they ride the body's rotation in both the inertial and planet-fixed views.
  const groundTracks = createGroundTracks(globe, frames, sidecar, {
    equatorialRadiusKm: Re, polarRadiusKm: Rp, pixelRatio,
    enabled: options.ground_tracks ?? false,
  });
  globe.group.add(groundTracks.group);

  // The ensemble panel needs `state.select` before `state` is built; it reads
  // through this reference, which is filled in below.
  const stateRef = { selected: -1, select: (i) => state.select(i), plotSeries: (spec) => state.plotSeries(spec) };

  const ensemble = ensembleSpec ? createEnsemble(ensembleSpec, frames, craft, container, stateRef, {}) : null;
  if (ensemble) world.add(ensemble.group);

  // Sun lighting: the directional sun follows frames.sun_dir, shadows are cast
  // around the followed vehicle, and 'path traced when paused' hands the scene
  // to three-gpu-pathtracer once the timeline and the camera are still.
  const lighting = createLighting(scene, renderer, frames, planet, {
    camera, world, globe, terrain, lod,
    helpers: [craft.group, atmosphere && atmosphere.group, refPaths.group, refs.group, ensemble && ensemble.group],
  });

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
    groundTracks: groundTracks.visible,
    setGroundTracks(v) { groundTracks.setVisible(v); state.groundTracks = groundTracks.visible; },
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
      if (i !== state.selected) {
        state.setFace(null);
        if (plots.key && plots.key.startsWith('sc')) { plots.hide(); ui.setPlotted(null); }
      }
      state.selected = i;
      stateRef.selected = i;
      if (i < 0 && state.follow) state.setFollow(false);
      ui.render();
    },
    // Picked face: { s, link, normal, point, source } from lod.pickFace, or null.
    face: null,
    setFace(pick) {
      state.face = pick;
      faceQuantities = pick ? buildFaceQuantities(pick) : [];
      lod.setFaceMarker(pick);
      if (!pick && plots.key && plots.key.includes(':face:')) { plots.hide(); ui.setPlotted(null); }
      ui.setFace(faceInfo(timeline.t));
    },
    // Open (or close) the time history of a panel quantity by its row key.
    plot(key) {
      const q = spacecraftQuantities.find((x) => x.key === key) || faceQuantities.find((x) => x.key === key);
      if (!q) return;
      plots.toggle(seriesOf(q));
      ui.setPlotted(plots.key);
    },
    plotSeries(spec) { plots.toggle(spec); ui.setPlotted(plots.key); },
    setTrailOrbits(n) { state.trailOrbits = n; applyTrail(); },
    setLabels(v) { craft.setLabelsVisible(v); refs.setLabelsVisible(v); },
    setGraticule(v) { globe.grid.visible = v; globe.axis.visible = v; },
    setAssemblies(v) { lod.setEnabled(v); },
    setThrusters(v) { lod.setThrustersVisible(v); },
    hasPlumes: plumes.available,
    setPlumes(v) { plumes.setVisible(v); },
    setFacets(v) { lod.setFacetsVisible(v); },
    setAxes(v) { lod.setAxesVisible(v); },
    hasDust: frames.hasPlume(),
    setDust(v) { dust.setVisible(v); },
    hasHeating: lod.heatingAvailable,
    setHeating(v) { lod.setHeatingVisible(v); ui.setHeatLegend(v && lod.heatingAvailable ? lod.heatRange() : null); },
    resetView() { state.setFollow(false); placeCamera(); },
    lightingModes: lighting.modes,
    lightingMode: lighting.mode,
    setLighting(m) { state.lightingMode = lighting.setMode(m); },
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
    lighting: lighting.status,
    'link poses': frames.linkPose ? 'recorded' : 'configured',
    ...(ensemble ? { samples: `${ensembleSpec.count} × ${ensemble.perSample} spacecraft` } : {}),
    ...(refPaths.items.length ? { paths: refPaths.items.map((it) => it.spec.name).join(', ') } : {}),
    ...(refs.items.length ? { references: refs.items.map((it) => it.spec.name).join(', ') } : {}),
    ...(terrain.levels.length ? { terrain: `${terrain.site ? terrain.site.name + ', ' : ''}${terrain.modelStatus}` } : {}),
    ...(terrain.levels.length && terrain.attribution ? { imagery: terrain.attribution } : {}),
    ...(atmosphere ? {
      atmosphere: `${atmosphere.info.model.replace('AtmosphereModel', '')}, EI ${atmosphere.info.ei_km.toFixed(0)} km`,
      ...(atmosphere.info.map ? { 'density map': `${atmosphere.info.map.altitude_km.toFixed(0)} km, ${atmosphere.info.map.min.toExponential(1)}–${atmosphere.info.map.max.toExponential(1)} kg/m³` } : {}),
    } : {}),
  });
  // Default trail coloring: heat rate when the run has it, else age.
  const initialColor = options.trail_color && TRAIL_COLOR_MODES[options.trail_color] ? options.trail_color : (frames.hasScalar('heat_rate') ? 'heat_rate' : 'age');
  ui.setTrailLegend(craft.setTrailColorMode(initialColor));
  // The graticule is drawn on the reference sphere. A landing site whose ground
  // lies below that radius (Tranquility Base is 1.5 km under it) would wear the
  // lines as bright streaks across the landscape, so a run with terrain starts
  // with the graticule off; the toolbar still turns it back on.
  if (terrain.levels.length) {
    const graticuleBox = ui.root.querySelector('[data-role="grid"]');
    if (graticuleBox) graticuleBox.checked = false;
    state.setGraticule(false);
  }
  if (lod.heatingAvailable && (options.heating ?? true)) { lod.setHeatingVisible(true); ui.setHeatLegend(lod.heatRange()); }
  lighting.ready.then(() => ui.setLightingModes(lighting.modes));
  const trailColorSelect = container.querySelector('[data-role="trailcolor"]');
  if (trailColorSelect) trailColorSelect.value = craft.trailColorMode;

  const tmpPos = new Float64Array(3), tmpVel = new Float64Array(3), qPi = new Float32Array(4), rBody = new Float64Array(3);
  const v3 = new THREE.Vector3();
  const plots = createPlotPanel(container, timeline, frames);

  // ---- Quantities. Every number in the selection and face panels is a
  // quantity: a sampler at(t, out) filling `dims` values, a row formatter,
  // and a unit, so the same definition feeds the panel text and the time
  // history the plot panel draws when the row is clicked.
  const fmt = (v, d = 1) => (Number.isFinite(v) ? v.toFixed(d) : '–');
  const vecText = (o, d, n) => { const parts = []; for (let i = 0; i < n; i++) parts.push(fmt(o[i], d)); return parts.join(', '); };
  const geo3 = new Float64Array(3), qAtt = new Float32Array(4), qLnk = new Float32Array(4), pLnk = new Float64Array(3), air = new Float64Array(3);
  const QA = new THREE.Quaternion(), QI = new THREE.Quaternion(), QL = new THREE.Quaternion(), NV = new THREE.Vector3(), FV = new THREE.Vector3();
  let spacecraftQuantities = [], spacecraftQuantitiesFor = -1, faceQuantities = [];

  function geodeticAt(t, s, o) {
    frames.positionAt(t, s, tmpPos);
    globe.rotationAt(t, qPi);
    rotateByConjugate(qPi, tmpPos, rBody);
    const g = geodetic(rBody, Re, Rp);
    o[0] = g.alt; o[1] = g.latDeg; o[2] = g.lonDeg;
    return o;
  }

  function buildSpacecraftQuantities(s) {
    const list = [];
    const add = (key, label, unit, dims, at, text, opts = {}) => list.push({ key: `sc${s}:${key}`, label, unit, dims, at, text, names: opts.names || null, log: opts.log ?? false });
    add('altitude', 'altitude', 'km', 1, (t, o) => { geodeticAt(t, s, geo3); o[0] = geo3[0]; }, (o) => `${fmt(o[0], 1)} km`);
    add('latitude', 'latitude', '°', 1, (t, o) => { geodeticAt(t, s, geo3); o[0] = geo3[1]; }, (o) => `${fmt(o[0], 3)}°`);
    add('longitude', 'longitude', '°', 1, (t, o) => { geodeticAt(t, s, geo3); o[0] = geo3[2]; }, (o) => `${fmt(o[0], 3)}°`);
    add('radius', 'radius', 'km', 1, (t, o) => { frames.positionAt(t, s, tmpPos); o[0] = Math.hypot(tmpPos[0], tmpPos[1], tmpPos[2]); }, (o) => `${fmt(o[0], 1)} km`);
    if (terrain.hasHeights) {
      // height above the terrain grids (the radar altitude), from the body-fixed position
      add('terrain_altitude', 'height above terrain', 'm', 1, (t, o) => {
        frames.positionAt(t, s, tmpPos); globe.rotationAt(t, qPi); rotateByConjugate(qPi, tmpPos, rBody);
        const r = Math.hypot(rBody[0], rBody[1], rBody[2]);
        const lat = THREE.MathUtils.radToDeg(Math.asin(rBody[2] / r)), lon = THREE.MathUtils.radToDeg(Math.atan2(rBody[1], rBody[0]));
        const h = terrain.heightAt(lat, lon);
        o[0] = Number.isFinite(h) ? 1000 * (r - terrain.referenceRadiusKm) - h : NaN;
      }, (o) => (Number.isFinite(o[0]) ? `${fmt(o[0], 1)} m` : 'off the terrain grids'), { log: 'auto' });
    }
    add('speed', 'speed', 'km/s', 1, (t, o) => { frames.velocityAt(t, s, tmpVel); o[0] = Math.hypot(tmpVel[0], tmpVel[1], tmpVel[2]); }, (o) => `${fmt(o[0], 3)} km/s`);
    add('airspeed', 'airspeed', 'km/s', 1, (t, o) => { lod.airspeedAt(t, s, air); o[0] = Math.hypot(air[0], air[1], air[2]); }, (o) => `${fmt(o[0], 3)} km/s`);
    if (frames.mass) add('mass', 'mass', 'kg', 1, (t, o) => { o[0] = frames.massAt(t, s); }, (o) => `${fmt(o[0], 1)} kg`);
    if (frames.hasScalar('density')) {
      add('density', 'density', 'kg/m³', 1, (t, o) => { o[0] = frames.scalarAt('density', t, s, Re); }, (o) => (o[0] > 0 ? `${o[0].toExponential(2)} kg/m³` : '0 (above atmosphere)'), { log: 'auto' });
      if (frames.hasScalar('dynamic_pressure')) add('dynamic_pressure', 'dyn. pressure', 'Pa', 1, (t, o) => { o[0] = frames.scalarAt('dynamic_pressure', t, s, Re); }, (o) => `${fmt(o[0], 3)} Pa`, { log: 'auto' });
    }
    if (frames.hasScalar('heat_rate')) add('heat_rate', 'heat rate', 'W/cm²', 1, (t, o) => { o[0] = frames.scalarAt('heat_rate', t, s, Re) / 1e4; }, (o) => `${fmt(o[0], 4)} W/cm²`, { log: 'auto' });
    if (frames.hasScalar('drag')) add('drag', 'drag', 'N', 1, (t, o) => { o[0] = frames.scalarAt('drag', t, s, Re); }, (o) => `${fmt(o[0], 3)} N`, { log: 'auto' });
    if (frames.hasScalar('wind')) add('wind', 'wind', 'm/s', 1, (t, o) => { o[0] = frames.scalarAt('wind', t, s, Re); }, (o) => `${fmt(o[0], 1)} m/s`);
    if (frames.hasPlume()) {
      // Plume-surface interaction of the descent engine with the regolith.
      const plume = (key, label, unit, name, digits, opts = {}) =>
        add(key, label, unit, 1, (t, o) => { o[0] = frames.plumeAt(name, t, s); }, (o) => `${fmt(o[0], digits)} ${unit}`, opts);
      plume('plume_height', 'engine height', 'm', 'height_m', 1, { log: 'auto' });
      plume('plume_shear', 'plume shear', 'Pa', 'shear_pa', 3, { log: 'auto' });
      plume('plume_pressure', 'plume pressure', 'Pa', 'pressure_pa', 1, { log: 'auto' });
      plume('plume_erosion', 'erosion rate', 'kg/s', 'erosion_kg_s', 2);
      plume('plume_eroded', 'eroded mass', 'kg', 'eroded_kg', 1);
      plume('plume_ejecta', 'ejecta speed', 'm/s', 'ejecta_mps', 1);
      plume('plume_ground_effect', 'ground effect', 'N', 'ground_effect_n', 1);
    }
    // One row per thruster: its firing level over the run, 0 (idle) to 1 (full).
    const nThrusters = plumes.counts(s);
    if (nThrusters > 0) {
      const levelBuf = new Float32Array(nThrusters);
      for (let k = 0; k < nThrusters; k++) {
        add(`thruster_${k + 1}`, `thruster ${k + 1} level`, '', 1,
          (t, o) => { const L = plumes.levelsAt(t, s, levelBuf); o[0] = L ? L[k] : NaN; },
          (o) => fmt(o[0], 3));
      }
    }
    // Extra channels the run asked the page to carry (export_visualization's
    // `channels`): one panel row and one plot each, in the order given.
    frames.channels.forEach((c, k) => {
      add(`ch_${c.name}`, c.label, c.unit, 1, (t, o) => { o[0] = frames.channelAt(k, t, s); },
        (o) => `${fmt(o[0], c.digits)}${c.unit ? ' ' + c.unit : ''}`, { log: c.log });
    });
    refs.items.forEach((it, k) => {
      if (it.target !== s) return;
      add(`ref${k}`, `vs ${it.spec.name}`, 'km', 1, (t, o) => { o[0] = refs.separationKm(t, k); },
        (o) => (Number.isFinite(o[0]) ? (o[0] >= 10 ? `${o[0].toFixed(1)} km` : `${(1000 * o[0]).toFixed(1)} m`) : 'not covered'), { log: 'auto' });
    });
    return list;
  }

  // Face quantities: the picked face's normal follows its link's pose, the
  // flow is the airspeed in body axes, and the local heating is the same
  // ½ρV³cosθ the overlay shades (full accommodation, no shadowing).
  function buildFaceQuantities(face) {
    const s = face.s, k = face.link;
    const nLink = new THREE.Vector3(face.normal[0], face.normal[1], face.normal[2]);
    const list = [];
    const add = (key, label, unit, dims, at, text, opts = {}) => list.push({ key: `sc${s}:face:${key}`, label, unit, dims, at, text, names: opts.names || null, log: opts.log ?? false });
    const normalBody = (t) => { lod.linkPoseAt(t, s, k, pLnk, qLnk); QL.set(qLnk[0], qLnk[1], qLnk[2], qLnk[3]); return NV.copy(nLink).applyQuaternion(QL); };
    const flowBody = (t) => {
      lod.airspeedAt(t, s, air); lod.attitudeAt(t, s, qAtt);
      QI.set(qAtt[0], qAtt[1], qAtt[2], qAtt[3]).invert();
      return FV.set(air[0], air[1], air[2]).normalize().applyQuaternion(QI);
    };
    const cosTheta = (t) => { const f = flowBody(t); return normalBody(t).dot(f); };
    add('incidence', 'incidence θ', '°', 1, (t, o) => { o[0] = THREE.MathUtils.radToDeg(Math.acos(Math.min(1, Math.max(-1, cosTheta(t))))); }, (o) => `${fmt(o[0], 1)}° ${o[0] < 90 ? '(windward)' : '(lee)'}`);
    if (frames.hasScalar('density')) {
      add('flux', 'local heat flux ½ρV³cosθ', 'W/cm²', 1, (t, o) => {
        const rho = frames.scalarAt('density', t, s, Re); const c = Math.max(0, cosTheta(t)); const v = Math.hypot(air[0], air[1], air[2]) * 1000;
        o[0] = rho > 0 ? 0.5 * rho * v * v * v * c / 1e4 : 0;
      }, (o) => `${o[0] > 0 ? o[0].toExponential(3) : '0'} W/cm²`, { log: 'auto' });
      add('ram', 'ram pressure ρV²cos²θ', 'Pa', 1, (t, o) => {
        const rho = frames.scalarAt('density', t, s, Re); const c = Math.max(0, cosTheta(t)); const v = Math.hypot(air[0], air[1], air[2]) * 1000;
        o[0] = rho > 0 ? rho * v * v * c * c : 0;
      }, (o) => `${o[0] > 0 ? o[0].toExponential(3) : '0'} Pa`, { log: 'auto' });
    }
    add('airspeed', 'airspeed', 'km/s', 1, (t, o) => { lod.airspeedAt(t, s, air); o[0] = Math.hypot(air[0], air[1], air[2]); }, (o) => `${fmt(o[0], 3)} km/s`);
    add('flow', 'flow in body axes', '', 3, (t, o) => { const f = flowBody(t); o[0] = f.x; o[1] = f.y; o[2] = f.z; }, (o) => vecText(o, 3, 3), { names: ['x', 'y', 'z'] });
    add('normal', 'face normal (body)', '', 3, (t, o) => { const n = normalBody(t); o[0] = n.x; o[1] = n.y; o[2] = n.z; }, (o) => vecText(o, 3, 3), { names: ['x', 'y', 'z'] });
    add('position', 'position (inertial)', 'km', 3, (t, o) => { frames.positionAt(t, s, tmpPos); o[0] = tmpPos[0]; o[1] = tmpPos[1]; o[2] = tmpPos[2]; }, (o) => vecText(o, 1, 3), { names: ['x', 'y', 'z'] });
    add('velocity', 'velocity (inertial)', 'km/s', 3, (t, o) => { frames.velocityAt(t, s, tmpVel); o[0] = tmpVel[0]; o[1] = tmpVel[1]; o[2] = tmpVel[2]; }, (o) => vecText(o, 3, 3), { names: ['x', 'y', 'z'] });
    add('attitude', `attitude q${frames.q ? '' : ' (velocity-aligned)'}`, '', 4, (t, o) => { lod.attitudeAt(t, s, qAtt); o[0] = qAtt[0]; o[1] = qAtt[1]; o[2] = qAtt[2]; o[3] = qAtt[3]; }, (o) => vecText(o, 3, 4), { names: ['x', 'y', 'z', 'w'] });
    if (k > 0) {
      add('linkpos', `link ${k + 1} position (body)`, 'm', 3, (t, o) => { lod.linkPoseAt(t, s, k, pLnk, qLnk); o[0] = pLnk[0]; o[1] = pLnk[1]; o[2] = pLnk[2]; }, (o) => vecText(o, 2, 3), { names: ['x', 'y', 'z'] });
      add('linkq', `link ${k + 1} q`, '', 4, (t, o) => { lod.linkPoseAt(t, s, k, pLnk, qLnk); o[0] = qLnk[0]; o[1] = qLnk[1]; o[2] = qLnk[2]; o[3] = qLnk[3]; }, (o) => vecText(o, 3, 4), { names: ['x', 'y', 'z', 'w'] });
    }
    return list;
  }

  // Full history of a quantity over the kept frames, as the plot panel wants it.
  function seriesOf(q) {
    const n = frames.count;
    const t = new Float64Array(n);
    for (let k = 0; k < n; k++) t[k] = frames.t[k];
    const ys = []; for (let d = 0; d < q.dims; d++) ys.push(new Float64Array(n));
    const o = new Float64Array(q.dims);
    for (let k = 0; k < n; k++) { o.fill(NaN); q.at(t[k], o); for (let d = 0; d < q.dims; d++) ys[d][k] = o[d]; }
    const s = state.selected;
    const owner = s >= 0 ? ((sidecar.spacecraft[s] || sidecar.spacecraft[0]).name || `sc${s + 1}`) : '';
    return { key: q.key, title: `${q.label}${owner ? ' · ' + owner : ''}`, unit: q.unit, t, series: ys.map((y, d) => ({ name: q.names ? q.names[d] : '', y })), log: q.log };
  }
  const qOut = new Float64Array(4);
  function quantityRows(list, t) {
    return list.map((q) => { qOut.fill(NaN); q.at(t, qOut); return { label: q.label, text: q.text(qOut), key: q.key }; });
  }

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
    // Near the ground the camera must start above it: tilt the view down along local up.
    if (terrain.levels.length && state.selected >= 0) {
      frames.positionAt(timeline.t, state.selected, tmpPos);
      const up = new THREE.Vector3(tmpPos[0], tmpPos[1], tmpPos[2]).applyQuaternion(world.quaternion).normalize();
      const r = Math.hypot(tmpPos[0], tmpPos[1], tmpPos[2]);
      if (r - terrain.referenceRadiusKm < 20) {
        const horizontal = dir.clone().sub(up.clone().multiplyScalar(dir.dot(up)));
        if (horizontal.lengthSq() < 1e-6) horizontal.set(1, 0, 0).sub(up.clone().multiplyScalar(up.x));
        dir.copy(horizontal.normalize().multiplyScalar(0.8).add(up.multiplyScalar(0.6))).normalize();
      }
    }
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
    if (pickFaceAt(px, py)) return;
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

  // Face pick at canvas pixel (px, py): selects that spacecraft and its face. Returns the pick or null.
  const raycaster = new THREE.Raycaster();
  const ndc = new THREE.Vector2();
  function pickFaceAt(px, py) {
    if (!lod.enabled) return null;
    const rect = renderer.domElement.getBoundingClientRect();
    ndc.set((px / rect.width) * 2 - 1, -(py / rect.height) * 2 + 1);
    world.updateMatrixWorld();
    raycaster.setFromCamera(ndc, camera);
    const pick = lod.pickFace(raycaster);
    if (!pick) return null;
    state.select(pick.s);
    state.setFace(pick);
    return pick;
  }

  function selectionInfo(t) {
    const s = state.selected;
    if (s < 0) return null;
    const spec = sidecar.spacecraft[s] || sidecar.spacecraft[0];
    const title = spec.name || `sc${s + 1}`;
    frames.positionAt(t, s, tmpPos);
    if (!Number.isFinite(tmpPos[0])) {
      return { title, rows: [{ label: 'id', text: String(spec.id) }, { label: 'state', text: 'not present at this time' }] };
    }
    if (spacecraftQuantitiesFor !== s) { spacecraftQuantities = buildSpacecraftQuantities(s); spacecraftQuantitiesFor = s; }
    const rows = [{ label: 'id', text: String(spec.id) }, ...quantityRows(spacecraftQuantities, t)];
    rows.push({ label: 'links', text: String(spec.links.length) });
    rows.push({ label: 'model', text: lod.visible[s] ? `3D (${lod.pxSize[s].toFixed(0)} px)` : `marker (${lod.pxSize[s].toFixed(1)} px)` });
    const status = lod.modelStatus(s);
    if (status) rows.push({ label: '3D model', text: status });
    refs.items.forEach((it, k) => {
      if (it.target !== s) return;
      const ghost = refs.modelStatus(k);
      if (ghost && ghost !== 'boxes') rows.push({ label: `${it.spec.name} model`, text: ghost });
    });
    return { title, rows };
  }

  function faceInfo(t) {
    const face = state.face;
    if (!face || face.s !== state.selected) return null;
    const spec = sidecar.spacecraft[face.s] || sidecar.spacecraft[0];
    const link = spec.links[face.link] || {};
    const rows = [
      { label: 'face', text: `${face.source === 'model' ? 'model surface' : 'box face'} on link ${face.link + 1}${link.name ? ` (${link.name})` : ''}` },
      { label: 'point (link frame)', text: `${vecText(face.point, 2, 3)} m` },
      ...quantityRows(faceQuantities, t),
    ];
    return { title: `face · ${spec.name || `sc${face.s + 1}`}`, rows };
  }

  function resize() {
    if (recording) return;
    const w = container.clientWidth || window.innerWidth, h = container.clientHeight || window.innerHeight;
    resizeTo(w, h);
    // The plot panel sits just above the toolbar, whose height depends on how its rows wrap.
    plots.panel.style.bottom = `${(ui.root.offsetHeight || 90) + 10}px`;
    // The selection panel scrolls instead of running into the toolbar (a lander lists every thruster).
    container.style.setProperty('--sa-toolbar', `${ui.root.offsetHeight || 90}px`);
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
  const shadowFocus = new THREE.Vector3();
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
    if (now - lastInfo > 100) { ui.setSelection(selectionInfo(timeline.t)); ui.setFace(faceInfo(timeline.t)); ui.setLightingStatus(lighting.status); lastInfo = now; }
    plots.update(timeline.t);
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
    if (terrain.levels.length) terrain.update(camera, viewportHeight);
    lod.update(t, camera, viewportHeight, lod.group.matrixWorld, anchor);
    dust.update(t);
    plumes.update(t);
    groundTracks.update(t);
    if (refPaths.items.length) refPaths.update(t, anchor);
    if (refs.items.length) refs.update(t, camera, viewportHeight, refs.group.matrixWorld, anchor);
    if (ensemble) ensemble.update(state.follow, t, anchor, camera, ensemble.group.matrixWorld);
    craft.update(t, camera, lod.markerHidden, state.selected, craft.group.matrixWorld, anchor, ensemble ? ensemble.dimMask : null);
    controls.update();
    // Shadow focus: the selected assembly in scene space (the origin in follow
    // mode), or the scene origin when nothing is selected.
    const selectedItem = state.selected >= 0 ? lod.items[state.selected] : null;
    lighting.update(t, selectedItem && selectedItem.group ? selectedItem.group.getWorldPosition(shadowFocus) : null);
    // A video export draws every frame itself, so it always takes the real-time
    // path: a path-traced frame would take seconds and never accumulate.
    if (recording || !lighting.render(scene, camera)) renderer.render(scene, camera);
  }
  requestAnimationFrame(animate);

  const viewer = {
    renderer, scene, camera, controls, world, globe, atmosphere, craft, lod, plumes, dust, ensemble, references: refs, paths: refPaths, groundTracks, timeline, frames, state, period, lighting,
    // Deterministic rendering for exports: seek and draw one frame at t.
    renderAt(t) { timeline.seek(t); frame(t); },
    setRecording(v) { recording = v; if (!v) resize(); },
    resizeTo,
    plots,
    terrain,
    pickFaceAt,
    quantities() { return { spacecraft: spacecraftQuantities.map((q) => q.key), face: faceQuantities.map((q) => q.key) }; },
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
