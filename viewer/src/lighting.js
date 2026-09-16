// Sun and earthshine lighting, the physical camera, and the optional
// path-traced renderer.
//
// The scene is lit from where the Sun actually was: `frames.sun_dir` is the
// unit vector from the planet's center to the Sun in the same inertial frame
// as `frames.pos_km` (written by the `sun_dir` save field, see
// src/simulation/callbacks/save_fields.jl), interpolated on the timeline and
// rotated into the scene exactly as positions are -- by the world group's
// quaternion, which is identity in the inertial frame and the inverse planet
// rotation in planet-fixed mode. The world group's translation (the floating
// origin of follow mode) moves positions, not directions, so it plays no part
// here; the shadow camera is instead re-centered on the followed vehicle every
// frame, because a directional shadow map has to be a few tens of meters wide
// to resolve a lander's legs.
//
// `frames.earth_dir` is the same vector for Earth, written away from Earth, and
// it carries the second light source a vacuum scene actually has: earthshine.
//
// Radiometry. The render buffer is kept in units of "a white Lambertian surface
// facing the Sun reads 1.0": the directional sun's intensity is pi, so a
// diffuse surface of albedo a at incidence theta reads a*cos(theta) (three.js
// shades diffuse as albedo * dotNL * intensity / pi). One scale factor,
// `radianceScale = sunIrradiance / pi`, turns a buffer value into an absolute
// radiance in W/m^2/sr, which is what the exposure is set against. Every other
// source is placed on the same scale, so their ratios are the physical ones:
//
//   sun         1361 / d_au^2 W/m^2, d from the body's mean Sun distance
//   earthshine  0.15 W/m^2 at full Earth, scaled by the Lambert-sphere phase
//               law for the Sun-Earth-body angle; about 1e-4 of sunlight, so
//               it is thirteen stops down and only shows at a long exposure
//   ground      the lit surface bounced back up: irradiance albedo*E*cos(theta)
//               on a downward-facing surface, carried by a hemisphere light
//               with a black sky (vacuum has none) and a regolith-colored
//               ground. This replaces the fill term that used to be a guess.
//
// Exposure is a stated camera setting, not a measurement of the scene. `ev` is
// the photographic exposure value at ISO 100 (EV = log2(L * S / K) with
// S = 100, K = 12.5, radiance converted to luminance at 98 lm/W, the luminous
// efficacy of unattenuated sunlight); the exposure factor is whatever maps the
// radiance that EV meters onto the tone curve's middle gray. The default EV is
// the one metered off a 0.12-albedo surface facing the Sun, which comes out at
// 15.3 at 1 au -- the "sunny 16" exposure, as it should. `[` and `]` step it by
// a third of a stop.
//
// Two modes:
//   'realtime'    the directional sun with a PCF-soft shadow map, earthshine,
//                 and the ground-bounce term.
//   'pathtraced'  while the timeline runs or the camera moves the real-time
//                 renderer draws; once nothing has changed for a moment the
//                 scene is handed to three-gpu-pathtracer, which accumulates
//                 samples into the same canvas until the next change. It gets
//                 the same irradiances (the sun as a disc of the true angular
//                 diameter, same total irradiance) and the same exposure. The
//                 library is loaded from the CDN on demand and the mode is
//                 simply not offered when it cannot be imported (the offline
//                 standalone page).
import * as THREE from 'three';

const LIGHTING_SUN_COLOR = 0xfff6ec;
const LIGHTING_FALLBACK_DIRECTION = [1, 0.3, 0.2];
// Intensities. The fallback (no sun_dir) keeps the flat lighting earlier pages
// were built with; a physically placed sun carries the irradiance instead.
const LIGHTING_FALLBACK_SUN_INTENSITY = 1.6;
const LIGHTING_FALLBACK_AMBIENT_INTENSITY = 0.55;
const LIGHTING_FALLBACK_HEMISPHERE_INTENSITY = 0.08;
// The directional intensity that puts the render buffer in "white Lambertian
// facing the Sun = 1.0" units, since three shades diffuse as a*dotNL*I/pi.
const LIGHTING_SUN_UNIT_INTENSITY = Math.PI;

// --- radiometry ------------------------------------------------------------
// Solar irradiance at 1 au (W/m^2), and the mean Sun distance of the bodies the
// viewer draws. The planet spec carries no Sun distance, so it is a table; an
// unlisted body is assumed to be at 1 au and `options.sunDistanceAu` overrides.
const LIGHTING_SOLAR_CONSTANT_W_M2 = 1361.0;
const LIGHTING_SUN_DISTANCE_AU = {
  mercury: 0.387, venus: 0.723, earth: 1.0, moon: 1.0, mars: 1.523,
  jupiter: 5.203, saturn: 9.537, titan: 9.537, uranus: 19.19, neptune: 30.07,
};
// Earthshine at the Moon with a full Earth (W/m^2). The Earth is about 3.7
// times the Moon's diameter and 2.5 times as bright per unit area, so a full
// Earth is some fifty times a full Moon; 0.15 W/m^2 is the round number that
// gives. It is calibrated at the Moon's 384,400 km, and the columns carry only
// a direction, so it is applied only to a body in the Earth's neighborhood --
// from Mars the Earth returns something like 1e-5 W/m^2, which is nothing.
const LIGHTING_EARTHSHINE_FULL_W_M2 = 0.15;
const LIGHTING_EARTHSHINE_MAX_AU_OFFSET = 0.05;
const LIGHTING_EARTH_COLOR = 0xc6d8ff;      // ocean and cloud, slightly blue
// The ground bounce: the albedo the surface is assumed to have when its own lit
// radiance is fed back as fill, and the color that fill carries.
const LIGHTING_GROUND_ALBEDO = 0.12;
const LIGHTING_BOUNCE_COLOR = 0xffe9d5;
const LIGHTING_SKY_COLOR = 0x000000;        // vacuum: no sky term at all

// --- the camera ------------------------------------------------------------
// Photographic exposure at ISO 100: EV = log2(L * S / K), with the reflected
// light meter constant K = 12.5 and the radiance converted to luminance at the
// luminous efficacy of sunlight above the atmosphere (1361 W/m^2 is about
// 133 klx, so 98 lm/W).
const LIGHTING_LUMINOUS_EFFICACY = 98.0;
const LIGHTING_METER_CONSTANT = 12.5;
const LIGHTING_METER_ISO = 100.0;
// The reference the default EV is metered off, and where it should land on the
// display: a 0.12-albedo surface facing the Sun renders at 0.18 after tone
// mapping. `renderer.toneMappingExposure` is then whatever puts the metered
// radiance at the tone curve's middle gray.
const LIGHTING_REFERENCE_ALBEDO = 0.12;
const LIGHTING_DISPLAY_TARGET = 0.18;
const LIGHTING_EV_STEP = 1 / 3;
const LIGHTING_EV_MIN = -8;
const LIGHTING_EV_MAX = 24;

// Shadow camera: an orthographic box `radius` wide around the followed vehicle,
// with the light `distance` up-sun of it. Kilometers, the scene's unit.
const LIGHTING_SHADOW_RADIUS_KM = 0.05;
const LIGHTING_SHADOW_DISTANCE_KM = 0.5;
const LIGHTING_SHADOW_MAP_SIZE = 2048;
// Biases sized for a grazing sun, which is the hard case and the interesting
// one: at the 10.6 degrees Apollo 11 landed under, a shadow texel (5 cm here)
// spans 26 cm of depth, so the default biases leave the whole terrain
// self-shadowed and black. These clear that and still keep the lander's
// shadow attached to its legs.
const LIGHTING_SHADOW_BIAS = -5e-4;
const LIGHTING_SHADOW_NORMAL_BIAS = 2e-4;   // 20 cm along the receiver's normal
// The Sun subtends 0.5334 degrees from 1 au; the penumbra that gives is what
// separates a rendered shadow from a stencil. The path tracer gets it from a
// circular area light of the matching solid angle (its directional light is an
// ideal parallel source), placed far enough up-sun that no terrain in between
// can occlude it.
const LIGHTING_SUN_ANGULAR_DIAMETER_DEG = 0.5334;
const LIGHTING_SUN_DISC_DISTANCE_KM = 5.0;
const LIGHTING_PATHTRACER_SPECIFIER = 'three-gpu-pathtracer';
const LIGHTING_SETTLE_MS = 300;
const LIGHTING_MODE_LABELS = { realtime: 'lighting: real-time', pathtraced: 'lighting: path traced when paused' };

// Environment probe for the metal foils: a small cube render of the lit ground
// and the black sky around the vehicle, prefiltered into a PMREM. It is rebuilt
// when the Sun has turned by more than a few degrees or the vehicle has moved
// far enough for the ground to fill a different part of its sky.
const LIGHTING_ENV_SIZE = 64;
const LIGHTING_ENV_SUN_STEP_DEG = 3.0;
const LIGHTING_ENV_MOVE_KM = 0.5;

// Sun direction at elapsed time t, linearly interpolated between the kept
// frames and renormalized. Returns false when the run carries none.
function lightingDirAt(table, frames, t, out) {
  const d = table;
  if (!d) return false;
  const { i, f } = frames.locate(t);
  const a = 3 * i, b = frames.count < 2 ? a : 3 * (i + 1);
  let x = d[a] + f * (d[b] - d[a]);
  let y = d[a + 1] + f * (d[b + 1] - d[a + 1]);
  let z = d[a + 2] + f * (d[b + 2] - d[a + 2]);
  const n = Math.hypot(x, y, z);
  if (!(n > 0)) return false;
  out.set(x / n, y / n, z / n);
  return true;
}

function lightingSunDirAt(frames, t, out) { return lightingDirAt(frames.sunDir, frames, t, out); }
function lightingEarthDirAt(frames, t, out) { return lightingDirAt(frames.earthDir, frames, t, out); }

// Phase law of a Lambertian sphere: the fraction of its full-phase brightness
// at phase angle alpha (0 = fully lit face toward the observer).
function lightingLambertPhase(alpha) {
  const a = Math.min(Math.PI, Math.max(0, alpha));
  return (Math.sin(a) + (Math.PI - a) * Math.cos(a)) / Math.PI;
}

// three.js's ACES filmic curve on the neutral axis. Its input and output
// matrices both have unit row sums, so a gray stays gray and the whole curve
// reduces to the RRT/ODT rational fit; inverting it numerically is how the
// exposure is turned into "this radiance lands at this display value".
function lightingAcesNeutral(x) {
  const v = Math.max(0, x);
  const y = (v * (v + 0.0245786) - 0.000090537) / (v * (0.983729 * v + 0.4329510) + 0.238081);
  return Math.min(1, Math.max(0, y));
}

function lightingAcesInverse(target) {
  let lo = 0, hi = 16;
  for (let k = 0; k < 60; k++) {
    const mid = 0.5 * (lo + hi);
    if (lightingAcesNeutral(mid) < target) lo = mid; else hi = mid;
  }
  return 0.5 * (lo + hi);
}

// Objects the path tracer must not see. Everything that is not a Mesh is
// ignored by the scene generator already (trails, graticules, axes and arrow
// helpers are lines; labels are sprites), so this covers the meshes that are
// user interface rather than vehicle or ground: markers, thruster and facet
// glyphs, the site ring, the density shells.
function lightingIsHelper(object) {
  return object.userData && object.userData.uiOnly === true;
}

export function createLighting(scene, renderer, frames, planet, options = {}) {
  const world = options.world || null;
  const camera = options.camera || null;
  const globe = options.globe || null;
  const terrain = options.terrain || null;
  const lod = options.lod || null;
  const helpers = (options.helpers || []).filter(Boolean);
  const hasSunDir = !!frames.sunDir;
  const hasEarthDir = !!frames.earthDir;

  renderer.shadowMap.enabled = true;
  renderer.shadowMap.type = THREE.PCFSoftShadowMap;
  renderer.shadowMap.autoUpdate = true;

  // ---- radiometry --------------------------------------------------------
  const bodyKeys = [planet.texture, planet.name].filter(Boolean).map((k) => String(k).toLowerCase());
  const sunDistanceAu = options.sunDistanceAu
    ?? bodyKeys.map((k) => LIGHTING_SUN_DISTANCE_AU[k]).find((v) => v != null)
    ?? 1.0;
  const sunIrradiance = LIGHTING_SOLAR_CONSTANT_W_M2 / (sunDistanceAu * sunDistanceAu);
  let earthshineIrradiance = 0;

  const sun = new THREE.DirectionalLight(LIGHTING_SUN_COLOR, hasSunDir ? LIGHTING_SUN_UNIT_INTENSITY : LIGHTING_FALLBACK_SUN_INTENSITY);
  sun.castShadow = hasSunDir;
  sun.shadow.mapSize.set(LIGHTING_SHADOW_MAP_SIZE, LIGHTING_SHADOW_MAP_SIZE);
  sun.shadow.bias = options.shadowBias ?? LIGHTING_SHADOW_BIAS;
  sun.shadow.normalBias = options.shadowNormalBias ?? LIGHTING_SHADOW_NORMAL_BIAS;
  const shadowRadiusKm = options.shadowRadiusKm ?? LIGHTING_SHADOW_RADIUS_KM;
  const shadowDistanceKm = options.shadowDistanceKm ?? LIGHTING_SHADOW_DISTANCE_KM;
  {
    const c = sun.shadow.camera;
    c.left = -shadowRadiusKm; c.right = shadowRadiusKm;
    c.top = shadowRadiusKm; c.bottom = -shadowRadiusKm;
    c.near = Math.max(1e-4, shadowDistanceKm - 6 * shadowRadiusKm);
    c.far = shadowDistanceKm + 6 * shadowRadiusKm;
    c.updateProjectionMatrix();
  }
  scene.add(sun);
  scene.add(sun.target);

  // Earthshine: a second directional light toward Earth, its irradiance the
  // sunlight the Earth reflects onto this body at the phase it is in. It does
  // not cast shadows -- at 1e-4 of the Sun there is nothing to see in them.
  // Earthshine rides on the same radiometric scale as the sun, so it is offered
  // only when the sun is physically placed as well.
  const earthshineApplies = hasSunDir && hasEarthDir && Math.abs(sunDistanceAu - 1.0) < LIGHTING_EARTHSHINE_MAX_AU_OFFSET;
  const earth = new THREE.DirectionalLight(LIGHTING_EARTH_COLOR, 0);
  earth.castShadow = false;
  earth.visible = earthshineApplies;
  scene.add(earth);
  scene.add(earth.target);

  // Fill. Without sun_dir the page keeps the flat ambient it always had. With
  // it there is no ambient term at all: what lights the shadows is the ground
  // bounce below, and (far below that) earthshine.
  const ambient = new THREE.AmbientLight(0xffffff, hasSunDir ? 0 : LIGHTING_FALLBACK_AMBIENT_INTENSITY);
  scene.add(ambient);
  const hemisphere = new THREE.HemisphereLight(LIGHTING_SKY_COLOR, LIGHTING_BOUNCE_COLOR, hasSunDir ? 0 : LIGHTING_FALLBACK_HEMISPHERE_INTENSITY);
  scene.add(hemisphere);
  // The bounce color is a tint, not a brightness: divide it out so the light's
  // intensity is the irradiance the ground actually returns.
  const lightingBounceLuminance = (() => {
    const c = new THREE.Color(LIGHTING_BOUNCE_COLOR);
    return Math.max(1e-3, 0.2126 * c.r + 0.7152 * c.g + 0.0722 * c.b);
  })();

  // Scene-space directions, and the focus the shadow box is centered on.
  const sunInertial = new THREE.Vector3(...LIGHTING_FALLBACK_DIRECTION).normalize();
  const sunScene = new THREE.Vector3().copy(sunInertial);
  const earthInertial = new THREE.Vector3();
  const earthScene = new THREE.Vector3();
  const focus = new THREE.Vector3();
  const upScene = new THREE.Vector3(0, 0, 1);

  let mode = 'realtime';
  const modes = [{ value: 'realtime', label: LIGHTING_MODE_LABELS.realtime }];

  // ---- exposure ----------------------------------------------------------
  // A buffer value of 1 is the radiance of a white Lambertian surface facing
  // the Sun; this turns buffer values into W/m^2/sr.
  const radianceScale = sunIrradiance / LIGHTING_SUN_UNIT_INTENSITY;
  // The pre-tone-map value that lands on LIGHTING_DISPLAY_TARGET.
  const middleGray = lightingAcesInverse(LIGHTING_DISPLAY_TARGET);
  const lightingEvOfBuffer = (b) => Math.log2(Math.max(1e-12, b) * radianceScale * LIGHTING_LUMINOUS_EFFICACY * LIGHTING_METER_ISO / LIGHTING_METER_CONSTANT);
  const lightingBufferOfEv = (v) => Math.pow(2, v) * LIGHTING_METER_CONSTANT / (LIGHTING_METER_ISO * LIGHTING_LUMINOUS_EFFICACY * radianceScale);
  const defaultEv = lightingEvOfBuffer(LIGHTING_REFERENCE_ALBEDO);
  let ev = options.ev != null ? Number(options.ev) : defaultEv;
  let exposure = 1;

  function lightingApplyExposure() {
    if (!hasSunDir) { exposure = 1; return; }
    exposure = middleGray / lightingBufferOfEv(ev);
    renderer.toneMapping = THREE.ACESFilmicToneMapping;
    renderer.toneMappingExposure = exposure;
  }
  lightingApplyExposure();

  const lightingEvText = () => `EV ${ev >= 0 ? '+' : '−'}${Math.abs(ev).toFixed(1)}`;

  function lightingStatusText() {
    if (!hasSunDir) return 'real-time, fixed light (no sun_dir in this run)';
    const parts = [`real-time, sun from the run`, lightingEvText()];
    if (earthshineApplies) parts.push(`earthshine ${earthshineIrradiance.toExponential(1)} W/m²`);
    return parts.join(', ');
  }
  let status = lightingStatusText();

  // ---- path tracer -------------------------------------------------------
  let tracerLib = null;
  let tracerError = null;
  let tracer = null;
  let tracerScene = null;      // the scene handed over, so a change re-hands it
  let restore = [];            // [obj, key, value] to put back after a trace
  let lastChangeMs = (typeof performance !== 'undefined' ? performance.now() : Date.now());
  let signature = '';
  let tracing = false;

  const ready = (async () => {
    if (options.pathTracing === false) return false;
    try {
      tracerLib = await import(LIGHTING_PATHTRACER_SPECIFIER);
      if (!tracerLib || !tracerLib.WebGLPathTracer) throw new Error('WebGLPathTracer is not exported');
      modes.push({ value: 'pathtraced', label: LIGHTING_MODE_LABELS.pathtraced });
      return true;
    } catch (err) {
      tracerLib = null;
      tracerError = err && err.message ? err.message : String(err);
      console.info(`path tracing unavailable (${tracerError}); real-time lighting only`);
      return false;
    }
  })();

  // A circular light of the Sun's angular diameter, radiance calibrated so its
  // irradiance matches the directional sun it replaces: L = I / omega, with
  // omega the solid angle of the disc. The tracer therefore works in the same
  // buffer units as the raster path, and the same exposure applies to both.
  function lightingMakeSunDisc() {
    const Shaped = tracerLib && tracerLib.ShapedAreaLight;
    const halfAngle = 0.5 * THREE.MathUtils.degToRad(LIGHTING_SUN_ANGULAR_DIAMETER_DEG);
    const radius = LIGHTING_SUN_DISC_DISTANCE_KM * Math.tan(halfAngle);
    const solidAngle = Math.PI * halfAngle * halfAngle;
    if (!Shaped) return null;
    const disc = new Shaped(LIGHTING_SUN_COLOR, sun.intensity / solidAngle, 2 * radius, 2 * radius);
    disc.isCircular = true;
    return disc;
  }

  function lightingHide(object) {
    restore.push([object, 'visible', object.visible]);
    object.visible = false;
  }

  // Hand-off preparation: hide the user interface meshes, give the materials
  // that predate physically based shading a sane roughness (the path tracer
  // defaults a missing `roughness` to 0, i.e. a mirror), and swap any custom
  // shader overlay back to the material it replaced, since the material
  // uploader reads `color` off every material it is given.
  function lightingPrepareScene() {
    restore = [];
    const helperSet = new Set(helpers);
    // The globe sphere is cut open under the terrain by a shader `discard` the
    // path tracer never runs, so with terrain in the scene it would be an
    // unbroken sphere through the ground: leave it out.
    if (globe && globe.mesh && terrain && terrain.hole) helperSet.add(globe.mesh);
    scene.traverse((object) => {
      if (!object.visible) return;
      if (helperSet.has(object) || lightingIsHelper(object)) { lightingHide(object); return; }
      if (!object.isMesh) return;
      const materials = Array.isArray(object.material) ? object.material : [object.material];
      if (materials.some((m) => m && m.color === undefined)) {
        if (object.userData.baseMaterial) {
          restore.push([object, 'material', object.material]);
          object.material = object.userData.baseMaterial;
        } else {
          lightingHide(object);
          return;
        }
      }
      for (const m of (Array.isArray(object.material) ? object.material : [object.material])) {
        if (m && m.color !== undefined && m.roughness === undefined) m.roughness = 1;
      }
    });
    restore.push([scene, 'background', scene.background]);
    scene.background = new THREE.Color(0x000000);   // vacuum: no sky to light the shadows
  }

  function lightingRestoreScene() {
    for (let k = restore.length - 1; k >= 0; k--) restore[k][0][restore[k][1]] = restore[k][2];
    restore = [];
  }

  function lightingHandOver(sceneToTrace, cam) {
    if (!tracer) {
      // three r160's Scene has no `backgroundRotation` / `environmentRotation`
      // (r165 added them) and the path tracer reads both without guarding --
      // including on the empty scene its own constructor hands itself. Two
      // identity Eulers on the prototype are all it wants, and nothing in the
      // raster renderer looks at them.
      if (THREE.Scene.prototype.backgroundRotation === undefined) {
        THREE.Scene.prototype.backgroundRotation = new THREE.Euler();
        THREE.Scene.prototype.environmentRotation = new THREE.Euler();
      }
      tracer = new tracerLib.WebGLPathTracer(renderer);
      tracer.renderScale = options.pathTracingScale ?? 0.5;
      tracer.bounces = options.pathTracingBounces ?? 3;
      tracer.tiles.set(options.pathTracingTiles ?? 3, options.pathTracingTiles ?? 3);
      tracer.rasterizeScene = false;      // no raster pass, so the area sun never reaches WebGLRenderer
      tracer.renderToCanvas = true;
      tracer.renderDelay = 0;
      tracer.minSamples = 1;
      tracer.filterGlossyFactor = 0.5;
      tracer.multipleImportanceSampling = true;
    }
    lightingPrepareScene();
    // The lights are read into a texture during setScene, so the sun disc only
    // has to be in the scene for that one call; the raster sun steps aside for
    // it so the two are not both counted.
    const disc = lightingMakeSunDisc();
    const sunWasVisible = sun.visible;
    if (disc) {
      disc.position.copy(focus).addScaledVector(sunScene, LIGHTING_SUN_DISC_DISTANCE_KM);
      disc.lookAt(focus);
      scene.add(disc);
      sun.visible = false;
    }
    tracer.setScene(sceneToTrace, cam);
    if (disc) { scene.remove(disc); disc.dispose && disc.dispose(); sun.visible = sunWasVisible; }
    lightingRestoreScene();
    tracerScene = sceneToTrace;
  }

  function lightingStopTracing() {
    tracing = false;
    tracerScene = null;
  }

  // Everything that invalidates an accumulated image: the timeline moved, the
  // camera moved, or the canvas was resized.
  function lightingSignature(t) {
    if (!camera) return String(t);
    const e = camera.matrixWorld.elements;
    let s = `${t.toFixed(6)}|${renderer.domElement.width}x${renderer.domElement.height}`;
    for (let k = 0; k < 16; k++) s += `|${e[k].toFixed(9)}`;
    return s;
  }

  // ---- environment probe --------------------------------------------------
  // What a mirror at the vehicle sees: the lit ground below, black sky above.
  // The foil materials on the 3D model reflect it, which is the whole reason a
  // metal reads as metal instead of as a flat gray. One 64-pixel cube render
  // prefiltered into a PMREM, rebuilt only when the Sun or the vehicle has
  // moved enough to change what it sees.
  let envTarget = null;
  let envCamera = null;
  let pmrem = null;
  let envGenerated = null;
  let envTexture = null;
  const envSun = new THREE.Vector3();
  const envAt = new THREE.Vector3();
  const envPoint = new THREE.Vector3();
  const envPos = new Float64Array(3);
  let envBuilt = false;
  const envCos = Math.cos(THREE.MathUtils.degToRad(LIGHTING_ENV_SUN_STEP_DEG));

  // Where the probe stands: the followed vehicle when there is one, else the
  // first spacecraft placed into the scene the way positions are. Never the
  // scene origin -- that is the planet's center, and a probe rendered from
  // inside the body sees nothing but its own shell.
  function lightingProbePoint(t, cameraTarget) {
    if (cameraTarget) return envPoint.copy(focus);
    frames.positionAt(t, 0, envPos);
    if (!Number.isFinite(envPos[0])) return envPoint.copy(focus);
    envPoint.set(envPos[0], envPos[1], envPos[2]);
    if (world) envPoint.applyQuaternion(world.quaternion).add(world.position);
    return envPoint;
  }

  function lightingEnvironmentStale(at) {
    if (!hasSunDir || !lod || typeof lod.setEnvironmentMap !== 'function') return false;
    if (!envBuilt) return true;
    if (envSun.dot(sunScene) < envCos) return true;
    return envAt.distanceTo(at) > LIGHTING_ENV_MOVE_KM;
  }

  function lightingBuildEnvironment(at) {
    if (!envTarget) {
      envTarget = new THREE.WebGLCubeRenderTarget(LIGHTING_ENV_SIZE, { type: THREE.HalfFloatType });
      envCamera = new THREE.CubeCamera(1e-5, 1e5, envTarget);
      pmrem = new THREE.PMREMGenerator(renderer);
      pmrem.compileCubemapShader();
    }
    // The vehicles step aside: a probe inside the model would see its own
    // interior, and the interface meshes are not part of the environment.
    const hidden = [];
    const hide = (o) => { if (o && o.visible) { hidden.push(o); o.visible = false; } };
    if (lod.group) hide(lod.group);
    for (const h of helpers) hide(h);
    scene.traverse((o) => { if (o.isMesh && lightingIsHelper(o)) hide(o); });
    const background = scene.background;
    scene.background = null;
    const toneMapping = renderer.toneMapping, toneExposure = renderer.toneMappingExposure;
    renderer.toneMapping = THREE.NoToneMapping;          // the probe holds radiance, not display values
    renderer.toneMappingExposure = 1;
    envCamera.position.copy(at);
    try {
      envCamera.update(renderer, scene);
      const generated = pmrem.fromCubemap(envTarget.texture);
      if (envGenerated) envGenerated.dispose();
      envGenerated = generated;
      envTexture = generated.texture;
      lod.setEnvironmentMap(envTexture);
      envBuilt = true;
      envSun.copy(sunScene);
      envAt.copy(at);
    } catch (err) {
      console.warn('the environment probe could not be rendered; foils keep their flat look', err);
      envBuilt = true;    // do not retry every frame
    }
    renderer.toneMapping = toneMapping;
    renderer.toneMappingExposure = toneExposure;
    scene.background = background;
    for (const o of hidden) o.visible = true;
  }

  // ---- earthshine ---------------------------------------------------------
  function lightingUpdateEarthshine(t) {
    if (!earthshineApplies || !lightingEarthDirAt(frames, t, earthInertial)) {
      earthshineIrradiance = 0;
      earth.visible = false;
      return;
    }
    // Phase angle Sun-Earth-body: the body sees a full Earth when Earth lies
    // opposite the Sun in its sky, so cos(alpha) = -dot(sunDir, earthDir) (the
    // Sun is far enough that its direction from the Earth and from the body are
    // the same to a thousandth of a radian).
    const cosAlpha = Math.min(1, Math.max(-1, -sunInertial.dot(earthInertial)));
    const phase = lightingLambertPhase(Math.acos(cosAlpha));
    earthshineIrradiance = LIGHTING_EARTHSHINE_FULL_W_M2 * phase;
    earth.visible = earthshineIrradiance > 0;
    // Against the unit intensity rather than the sun light's own, so that
    // turning the sun off (to look at the earthshine alone) leaves this term
    // where it was.
    earth.intensity = LIGHTING_SUN_UNIT_INTENSITY * earthshineIrradiance / sunIrradiance;
    earthScene.copy(earthInertial);
    if (world) earthScene.applyQuaternion(world.quaternion);
    earth.target.position.copy(focus);
    earth.position.copy(focus).addScaledVector(earthScene, shadowDistanceKm);
    earth.target.updateMatrixWorld();
    earth.updateMatrixWorld();
  }

  // The sunlight the ground returns: a downward-facing surface just above a
  // Lambertian half-space of albedo a lit at incidence theta receives
  // a*E*cos(theta), which in buffer units is a*cos(theta)*sun.intensity.
  function lightingUpdateBounce() {
    if (!hasSunDir) return;
    const cosine = Math.max(0, upScene.dot(sunScene));
    hemisphere.intensity = LIGHTING_SUN_UNIT_INTENSITY * LIGHTING_GROUND_ALBEDO * cosine / lightingBounceLuminance;
  }

  const api = {
    sun,
    earth,
    ambient,
    hemisphere,
    ready,
    modes,
    get mode() { return mode; },
    get status() { return status; },
    get pathTracingAvailable() { return !!tracerLib; },
    get exposure() { return exposure; },
    get samples() { return tracing && tracer ? tracer.samples : 0; },
    // Physical quantities the other viewer modules light themselves by.
    get sunIrradiance() { return sunIrradiance; },
    get sunDistanceAu() { return sunDistanceAu; },
    get earthshineIrradiance() { return earthshineIrradiance; },
    get ev() { return ev; },
    get defaultEv() { return defaultEv; },
    get evText() { return lightingEvText(); },
    get environmentMap() { return envTexture; },

    // Camera exposure. `setEv` is the absolute photographic value; `stepEv`
    // moves it by thirds of a stop, which is what `[` and `]` do.
    setEv(next) {
      if (!Number.isFinite(next)) return ev;
      ev = Math.min(LIGHTING_EV_MAX, Math.max(LIGHTING_EV_MIN, next));
      lightingApplyExposure();
      status = mode === 'realtime' ? lightingStatusText() : status;
      lastChangeMs = performance.now();
      if (tracing) lightingStopTracing();
      return ev;
    },
    stepEv(steps = 1) { return api.setEv(ev + steps * LIGHTING_EV_STEP); },

    setMode(next) {
      const want = next === 'pathtraced' && tracerLib ? 'pathtraced' : 'realtime';
      if (want === mode) return mode;
      mode = want;
      if (mode === 'realtime') lightingStopTracing();
      lastChangeMs = performance.now();
      return mode;
    },

    // Place the sun for elapsed time `t` and aim the shadow box at
    // `cameraTarget` (the followed or selected vehicle in scene space, or null
    // for the scene origin).
    update(t, cameraTarget) {
      if (lightingSunDirAt(frames, t, sunInertial)) {
        sunScene.copy(sunInertial);
        if (world) sunScene.applyQuaternion(world.quaternion);
      } else {
        sunScene.copy(sunInertial);
      }
      if (cameraTarget) focus.copy(cameraTarget); else focus.set(0, 0, 0);
      sun.target.position.copy(focus);
      sun.position.copy(focus).addScaledVector(sunScene, shadowDistanceKm);
      sun.target.updateMatrixWorld();
      sun.updateMatrixWorld();
      // Local vertical of the focus, for the ground-bounce fill: the planet
      // center sits at the world group's origin.
      if (world) {
        upScene.copy(focus).sub(world.position);
        if (upScene.lengthSq() < 1e-12) upScene.copy(sunScene);
        upScene.normalize();
      }
      hemisphere.position.copy(upScene);
      hemisphere.updateMatrixWorld();
      lightingUpdateBounce();
      lightingUpdateEarthshine(t);

      const sig = lightingSignature(t);
      if (sig !== signature) {
        signature = sig;
        lastChangeMs = performance.now();
        if (tracing) lightingStopTracing();
      }
      const probeAt = lightingProbePoint(t, cameraTarget);
      if (lightingEnvironmentStale(probeAt)) lightingBuildEnvironment(probeAt);
      if (mode !== 'pathtraced') status = lightingStatusText();
    },

    // Draw the frame when this module owns the canvas. Returns false when the
    // caller should render normally, which is every frame in real-time mode and
    // every unsettled frame in path-traced mode.
    render(sceneToDraw, cam) {
      if (mode !== 'pathtraced' || !tracerLib) return false;
      const waited = performance.now() - lastChangeMs;
      if (waited < LIGHTING_SETTLE_MS) {
        status = `path traced, waiting for the view to settle, ${lightingEvText()}`;
        return false;
      }
      try {
        if (!tracing || tracerScene !== sceneToDraw) {
          lightingHandOver(sceneToDraw, cam);
          tracing = true;
        }
        tracer.renderSample();
        status = tracer.isCompiling
          ? `path traced, compiling the shader, ${lightingEvText()}`
          : `path traced, ${tracer.samples.toFixed(0)} samples, ${lightingEvText()}`;
        return true;
      } catch (err) {
        tracerError = err && err.message ? err.message : String(err);
        console.warn('path tracing failed; falling back to real-time', err);
        lightingRestoreScene();
        lightingStopTracing();
        tracerLib = null;
        mode = 'realtime';
        status = `path tracing failed: ${tracerError}`;
        return false;
      }
    },

    // Scene-space unit vector toward the Sun, for callers that want it.
    direction(out) { return out.copy(sunScene); },
    // The same toward Earth, or null when the run carries no Earth direction.
    earthDirection(out) { return earthshineApplies && earth.visible ? out.copy(earthScene) : null; },
  };
  return api;
}
