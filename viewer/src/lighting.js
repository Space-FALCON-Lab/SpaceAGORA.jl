// Sun lighting and the optional path-traced renderer.
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
// Two modes:
//   'realtime'    one directional sun with a PCF-soft shadow map, plus a faint
//                 hemisphere term standing in for earthshine (or sky glow).
//   'pathtraced'  while the timeline runs or the camera moves the real-time
//                 renderer draws; once nothing has changed for a moment the
//                 scene is handed to three-gpu-pathtracer, which accumulates
//                 samples into the same canvas until the next change. The
//                 library is loaded from the CDN on demand and the mode is
//                 simply not offered when it cannot be imported (the offline
//                 standalone page).
import * as THREE from 'three';

const LIGHTING_SUN_COLOR = 0xfff6ec;
const LIGHTING_FALLBACK_DIRECTION = [1, 0.3, 0.2];
// Intensities. The fallback (no sun_dir) keeps the flat lighting earlier pages
// were built with; a physically placed sun is brighter and its fill is faint,
// because vacuum has no sky to scatter light into the shadows.
const LIGHTING_FALLBACK_SUN_INTENSITY = 1.6;
const LIGHTING_FALLBACK_AMBIENT_INTENSITY = 0.55;
const LIGHTING_SUN_INTENSITY = 2.4;
const LIGHTING_HEMISPHERE_INTENSITY = 0.08;
const LIGHTING_AMBIENT_INTENSITY = 0.03;
// Exposure. A real surface lit at a grazing angle is genuinely dark: the LROC
// mosaic of Tranquility Base averages 0.021 in linear light, and at the 10.6
// degrees of the landing only a fifth of that reaches the eye, which leaves the
// ground -- and the shadows falling on it -- inside a handful of display levels.
// So the exposure is measured rather than assumed: the mean albedo of the ground
// texture times the sun's elevation at the site gives what a lit surface will
// render at, and the exposure is whatever puts that at LIGHTING_EXPOSURE_TARGET.
// The constant folds in the low-end gain of the ACES curve (0.30 as measured
// against the untone-mapped render), so the target is the linear value the
// display sees. A scene bright enough not to need the lift keeps the linear
// mapping, so pages of Earth and Mars look as they always did.
const LIGHTING_EXPOSURE_TARGET = 0.15;
const LIGHTING_EXPOSURE_MAX = 20.0;
const LIGHTING_EXPOSURE_MIN_USEFUL = 2.0;
const LIGHTING_SKY_COLOR = 0x2a3a52;
const LIGHTING_GROUND_COLOR = 0x6b6258;
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

// Sun direction at elapsed time t, linearly interpolated between the kept
// frames and renormalized. Returns false when the run carries none.
function lightingSunDirAt(frames, t, out) {
  const d = frames.sunDir;
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
  const helpers = (options.helpers || []).filter(Boolean);
  const hasSunDir = !!frames.sunDir;

  renderer.shadowMap.enabled = true;
  renderer.shadowMap.type = THREE.PCFSoftShadowMap;
  renderer.shadowMap.autoUpdate = true;
  let exposure = 1;
  if (hasSunDir && options.exposure != null) {
    exposure = options.exposure;
    renderer.toneMapping = THREE.ACESFilmicToneMapping;
    renderer.toneMappingExposure = exposure;
  }

  const sun = new THREE.DirectionalLight(LIGHTING_SUN_COLOR, hasSunDir ? LIGHTING_SUN_INTENSITY : LIGHTING_FALLBACK_SUN_INTENSITY);
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

  // Fill. Without sun_dir the page keeps the flat ambient it always had; with
  // it, a hemisphere term oriented along the local vertical stands in for
  // earthshine and the light the ground bounces back up into the vehicle.
  const ambient = new THREE.AmbientLight(0xffffff, hasSunDir ? LIGHTING_AMBIENT_INTENSITY : LIGHTING_FALLBACK_AMBIENT_INTENSITY);
  scene.add(ambient);
  const hemisphere = new THREE.HemisphereLight(LIGHTING_SKY_COLOR, LIGHTING_GROUND_COLOR, hasSunDir ? LIGHTING_HEMISPHERE_INTENSITY : 0);
  scene.add(hemisphere);

  // Scene-space sun direction, and the focus the shadow box is centered on.
  const sunInertial = new THREE.Vector3(...LIGHTING_FALLBACK_DIRECTION).normalize();
  const sunScene = new THREE.Vector3().copy(sunInertial);
  const focus = new THREE.Vector3();
  const upScene = new THREE.Vector3(0, 0, 1);

  let mode = 'realtime';
  let status = hasSunDir ? 'real-time, sun from the run' : 'real-time, fixed light (no sun_dir in this run)';
  const modes = [{ value: 'realtime', label: LIGHTING_MODE_LABELS.realtime }];

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
  // omega the solid angle of the disc.
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

  // ---- auto exposure ------------------------------------------------------
  let exposureMeasured = options.exposure != null || !hasSunDir;

  const hasTerrain = !!(terrain && terrain.levels && terrain.levels.length);

  // The surface the exposure is set for: the widest terrain patch when the run
  // has one (the ground the vehicle is on), else the globe's own map. A texture
  // that has not finished decoding yet reports nothing and is retried.
  function lightingGroundTexture() {
    if (hasTerrain) {
      const map = terrain.levels[0].mesh.material.map;
      return map && map.image && map.image.width ? map : null;
    }
    if (globe && globe.mesh && globe.mesh.material && globe.mesh.material.map) {
      const map = globe.mesh.material.map;
      if (map && map.image && map.image.width) return map;
    }
    return null;
  }

  // Cosine of the sun's incidence on the ground the run ends on -- for a landing
  // that is the site at touchdown, which is the moment the exposure has to serve.
  // Both vectors are inertial, so the camera and the frame toggle play no part
  // and the exposure is the same every time the page is opened.
  const lightingRefPos = new Float64Array(3);
  const lightingRefSun = new THREE.Vector3();
  function lightingReferenceCosine() {
    if (!hasTerrain) return 0.5;
    frames.positionAt(frames.tEnd, 0, lightingRefPos);
    if (!lightingSunDirAt(frames, frames.tEnd, lightingRefSun)) return 0.5;
    const radius = Math.hypot(lightingRefPos[0], lightingRefPos[1], lightingRefPos[2]);
    if (!(radius > 0)) return 0.5;
    const cosine = (lightingRefPos[0] * lightingRefSun.x + lightingRefPos[1] * lightingRefSun.y + lightingRefPos[2] * lightingRefSun.z) / radius;
    return Math.max(0.15, Math.abs(cosine));
  }

  // Mean linear luminance of a texture, from a 32x32 downsample.
  function lightingMeanAlbedo(texture) {
    try {
      const n = 32;
      const canvas = document.createElement('canvas');
      canvas.width = n; canvas.height = n;
      const ctx = canvas.getContext('2d', { willReadFrequently: true });
      ctx.drawImage(texture.image, 0, 0, n, n);
      const data = ctx.getImageData(0, 0, n, n).data;
      let sum = 0;
      for (let k = 0; k < data.length; k += 4) {
        const v = (0.2126 * data[k] + 0.7152 * data[k + 1] + 0.0722 * data[k + 2]) / 255;
        sum += v <= 0.04045 ? v / 12.92 : Math.pow((v + 0.055) / 1.055, 2.4);
      }
      return sum / (data.length / 4);
    } catch (err) {
      return NaN;   // a texture the page cannot read back leaves the exposure alone
    }
  }

  function lightingAutoExposure() {
    if (exposureMeasured) return;
    const texture = lightingGroundTexture();
    if (!texture) return;               // still loading; try again next frame
    exposureMeasured = true;
    const albedo = lightingMeanAlbedo(texture);
    if (!(albedo > 0)) return;
    const wanted = LIGHTING_EXPOSURE_TARGET / (albedo * sun.intensity * lightingReferenceCosine());
    if (!(wanted > LIGHTING_EXPOSURE_MIN_USEFUL)) return;
    exposure = Math.min(LIGHTING_EXPOSURE_MAX, wanted);
    renderer.toneMapping = THREE.ACESFilmicToneMapping;
    renderer.toneMappingExposure = exposure;
  }

  const api = {
    sun,
    ambient,
    hemisphere,
    ready,
    modes,
    get mode() { return mode; },
    get status() { return status; },
    get pathTracingAvailable() { return !!tracerLib; },
    get exposure() { return exposure; },
    get samples() { return tracing && tracer ? tracer.samples : 0; },

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
      // Local vertical of the focus, for the hemisphere fill: the planet center
      // sits at the world group's origin.
      if (world) {
        upScene.copy(focus).sub(world.position);
        if (upScene.lengthSq() < 1e-12) upScene.copy(sunScene);
        upScene.normalize();
      }
      hemisphere.position.copy(upScene);
      hemisphere.updateMatrixWorld();

      const sig = lightingSignature(t);
      if (sig !== signature) {
        signature = sig;
        lastChangeMs = performance.now();
        if (tracing) lightingStopTracing();
      }
      lightingAutoExposure();
      if (mode !== 'pathtraced') {
        status = hasSunDir
          ? `real-time, sun from the run${exposure > 1 ? `, exposure x${exposure.toFixed(1)}` : ''}`
          : 'real-time, fixed light (no sun_dir in this run)';
      }
    },

    // Draw the frame when this module owns the canvas. Returns false when the
    // caller should render normally, which is every frame in real-time mode and
    // every unsettled frame in path-traced mode.
    render(sceneToDraw, cam) {
      if (mode !== 'pathtraced' || !tracerLib) return false;
      const waited = performance.now() - lastChangeMs;
      if (waited < LIGHTING_SETTLE_MS) {
        status = 'path traced, waiting for the view to settle';
        return false;
      }
      try {
        if (!tracing || tracerScene !== sceneToDraw) {
          lightingHandOver(sceneToDraw, cam);
          tracing = true;
        }
        tracer.renderSample();
        status = tracer.isCompiling ? 'path traced, compiling the shader' : `path traced, ${tracer.samples.toFixed(0)} samples`;
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
  };
  return api;
}
