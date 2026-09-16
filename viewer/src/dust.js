// Blowing regolith under a landing vehicle: the ground-hugging sheet of gas and
// dust the descent engine's plume drives off the surface, the haze it fills the
// scene with, the ground it obscures, and the scour crater it digs. Driven
// entirely by the `frames.plume` block the plume-surface effector writes (see
// src/dynamics/coupled/force_torque_models/plume_surface_interaction.jl), so
// nothing appears where the physics says the soil is not moving.
//
// Payload contract (viewer_bundle.jl `build_viewer_frames`):
//   frames.plume.height_m        distance from the vehicle to the ground along the engine axis
//   frames.plume.shear_pa        peak wall shear stress under the plume
//   frames.plume.pressure_pa     peak surface pressure under the plume
//   frames.plume.erosion_kg_s    mass erosion rate
//   frames.plume.eroded_kg       its time integral
//   frames.plume.ejecta_mps      speed the grains leave at
//   frames.plume.ground_effect_n thrust augmentation in ground effect
// Each is a Float32 array, frame-major then spacecraft.
//
// What it is supposed to look like: NASA Langley's plume-surface interaction
// tests, where an engine fires into a bin of lunar simulant inside a vacuum
// chamber. Overhead, the first seconds show a continuous translucent veil
// streaming radially out of the impingement point in billowing filaments that
// fade with radius, not a spray of separate grains; from the side the sheet
// stays flat against the surface, only a couple of meters thick, while a fine
// haze builds up over the following seconds until the hardware behind it is
// barely readable. The page draws that as
//
//   * a stack of very flat disks whose fragment shader is fractional Brownian
//     motion advected radially at the ejecta speed, with angular streaks (the
//     sheet -- the main body of the effect);
//   * a thin, sparse set of points for the fast grains that outrun it;
//   * a scour crater: a swept-clear floor inside a darker rim of piled soil;
//   * a bank of view-aligned haze puffs whose opacity is a decaying integral
//     of the erosion rate, so it builds while the engine blows and thins
//     afterwards, and which veils the scene from a camera standing in it as
//     well as from one above;
//   * a soft elongated shadow the sheet casts on the ground away from the sun.
//
// Everything the sun touches is scaled by a Henyey-Greenstein phase function:
// fine dust scatters strongly forward, which is why the Apollo films are
// brightest looking down-sun through the sheet and nearly clear looking the
// other way. The sun direction, its irradiance and the camera's exposure come
// from the lighting handle (`options.lighting`); the shaders tone map and
// encode their output exactly as every lit material does, so one exposure
// applies to the whole frame, and without a handle -- or on a page that never
// raised its exposure -- the module keeps the fixed brightness it had before
// that module existed.
//
// Positions inside the sheet, the puffs and the points are derived in the
// shaders from the clock, so a frame update writes uniforms and nothing else.
//
// Every top-level name here is prefixed `dust`/`DUST_`: the CDN page builder
// concatenates all viewer modules into one script scope.
import * as THREE from 'three';
import { rotateByConjugate } from 'viewer/data.js';

const DUST_M_TO_KM = 1e-3;
const DUST_PARTICLE_COUNT = 4096;
const DUST_EJECTA_FRACTION = 0.12;                     // the points are the fast tail, not the sheet
const DUST_MAX_EMITTERS = 8;
const DUST_LIFETIME_S = 2.4;
const DUST_ELEVATION_MIN_RAD = 1.0 * Math.PI / 180;   // the sheet hugs the ground: a
const DUST_ELEVATION_MAX_RAD = 3.0 * Math.PI / 180;   // few degrees above the local slope
const DUST_GRAVITY_M_S2 = 1.62;                        // lunar, unless the payload says otherwise
const DUST_SIZE_KM = 4.5e-4;                           // a puff about half a meter across
const DUST_PROJECTION_SCALE = 820;                     // pixels per unit angle, ~800 px tall view at 45 deg
const DUST_POINT_MIN_PX = 1.2;
const DUST_POINT_MAX_PX = 16.0;   // a near grain is a speck, not a bubble: the sheet is the body of the effect
const DUST_COLOR = 0xc9bda6;
const DUST_SCOUR_FLOOR_COLOR = 0xa79e8d;               // soil swept down to a smoother, brighter floor
const DUST_SCOUR_RIM_COLOR = 0x574f43;                 // the ridge of piled-up material around it
const DUST_SHADOW_COLOR = 0x0b0d10;
const DUST_MAX_OPACITY = 0.11;                         // the points are a sparkle over the sheet, not a spray
const DUST_SCOUR_RIM_OPACITY = 0.3;
const DUST_SCOUR_FLOOR_OPACITY = 0.12;                 // added light, scaled by the display gain
const DUST_SCOUR_RIM_RADIUS = 0.78;                    // of the crater radius
const DUST_SCOUR_MAX_RADIUS_M = 14.0;
const DUST_SCOUR_RIM_WIDTH = 0.34;                     // as a fraction of the crater radius
const DUST_SHEET_LIFT_M = 0.7;                         // the saltating layer has thickness, and this
                                                       // keeps the sheet out of the terrain patches' facets
const DUST_DECAL_LIFT_M = 1.6;                         // clears the terrain patches' triangle facets
const DUST_PLUME_TAN_HALF_ANGLE = 0.4663;              // tan(25 deg), the plume core the effector models
const DUST_UP_AXIS = new THREE.Vector3(0, 0, 1);

// ---- the sheet -----------------------------------------------------------
// Three disks, a few tens of centimeters apart, each drawing the same advected
// noise field with its own seed and drift. Stacking them is what makes the
// veil read as a volume of dust from a low camera instead of a decal, and
// three is enough at the Apollo scale: the layer is 2.4 m thick against a
// 45 m radius, the aspect ratio the test films show.
const DUST_SHEET_LAYERS = [
  { height_m: 0.35, weight: 1.00, seed: 0.0, scroll: 1.00, swirl: 0.95 },
  { height_m: 1.10, weight: 0.78, seed: 4.7, scroll: 0.86, swirl: -0.65 },
  { height_m: 2.40, weight: 0.52, seed: 9.1, scroll: 0.72, swirl: 0.38 },
];
const DUST_SHEET_MAX_RADIUS_M = 45.0;
const DUST_SHEET_FOOTPRINTS = 2.2;                     // radius at the moment erosion starts
const DUST_SHEET_SPREAD = 0.18;                        // the visible front advances at this fraction of the ejecta speed
const DUST_SHEET_TAU = 2.8;                           // optical depth of a layer, straight through
const DUST_SHEET_EXTINCTION = 0.8;                     // how much of the sheet's alpha actually hides the ground
// Streaks: the noise is sampled on a circle of this radius, so it has about
// 2*pi*uSpokes cells around the sheet, against uRadialFreq cells from the
// center to the rim. Five to one is what draws filaments rather than blobs.
const DUST_SHEET_SPOKES = 5.0;
const DUST_SHEET_RADIAL_FREQ = 1.6;
const DUST_SHEET_SHARPNESS = 1.9;                      // exponent on the noise: lanes between the billows
const DUST_SHEET_DECAY = 1.9;                          // the veil thins outward, e-folding in sheet radii
const DUST_SHEET_FLOOR = 0.22;                         // ... but never quite to nothing inside the rim
const DUST_SHEET_SCROLL_FRACTION = 0.12;               // of the ejecta speed, in sheet radii per second
const DUST_SHEET_SEGMENTS = 96;
const DUST_SHEET_RINGS = 8;

// ---- the haze ------------------------------------------------------------
// The fines that never settle. They are drawn as a bank of view-aligned
// puffs -- flattened ellipsoids of dust, wide and low, each rendered as a quad
// turned to face the camera about the local vertical -- rather than as flat
// disks, because a camera at head height looks along the layer and a stack of
// horizontal sheets has nothing to show it, which is exactly the view the side
// footage is shot from. Opacity follows an exponentially decaying integral of
// the erosion rate, so the haze builds up over seconds and thins out after the
// engine stops. The chamber films hold theirs far longer than the Moon would
// -- a chamber has residual gas to keep the finest grains aloft -- so
// `hazeDecayS` and `hazeOpacity` are knobs, and their defaults are the vacuum
// case: thinner, and gone within a few tens of seconds.
const DUST_HAZE_PUFFS = 22;
const DUST_HAZE_DECAY_S = 9.0;
const DUST_HAZE_OPACITY = 0.30;
const DUST_HAZE_TAU = 1.3;                             // per puff, through its middle
const DUST_HAZE_RADIUS_M = 90.0;
const DUST_HAZE_HEIGHT_M = 16.0;                       // the bank's top, over the impingement point
const DUST_HAZE_PUFF_HEIGHT_M = 7.0;                   // half-height of one puff: wide and low
const DUST_HAZE_PUFF_MIN = 0.18;                       // puff radius, as a fraction of the bank's
const DUST_HAZE_PUFF_MAX = 0.40;
const DUST_HAZE_DRIFT_HZ = 0.055;                      // how fast a puff drifts out and is replaced
const DUST_HAZE_COLOR = 0xbdb3a1;
const DUST_HAZE_PHASE_G = 0.35;                        // multiply-scattered: much less directional than the sheet

// ---- the sun --------------------------------------------------------------
// Henyey-Greenstein asymmetry for lunar fines. g in 0.6 to 0.7 puts most of
// the scattered light within a few tens of degrees of the forward direction,
// which is the Apollo 12 and 17 observation: the sheet is blinding down-sun
// and almost invisible up-sun.
const DUST_PHASE_G = 0.65;
const DUST_PHASE_MIN = 0.25;                           // clamped for display: the true ratio is ~100:1
const DUST_PHASE_MAX = 3.0;
const DUST_REFERENCE_IRRADIANCE_W_M2 = 1361.0;
// The dust's own radiance, as a fraction of the light falling on it: a cloud
// of fines scatters a good deal of what reaches it, and unlike the ground it
// has no incidence cosine to cut it down, which is why the Apollo sheets are
// so much brighter than the surface they blow off. Multiplied by the scene's
// sun intensity and by the phase function, and then tone mapped by the same
// exposure as every other material.
const DUST_SCATTER_ALBEDO = 0.02;
// What a page with no physical exposure uses instead: the brightness the
// module had before there was a lighting module to ask.
const DUST_FALLBACK_BRIGHTNESS = 0.55;
const DUST_SCALE_MIN = 0.004;
const DUST_SCALE_MAX = 1.5;

// ---- the sheet's shadow ---------------------------------------------------
const DUST_SHADOW_OPACITY = 0.3;
const DUST_SHADOW_MAX_STRETCH = 6.0;                   // at a grazing sun the sheet's shadow runs a long way
const DUST_SHADOW_MIN_ELEVATION_RAD = 3.0 * Math.PI / 180;

// Value noise and the phase function, shared by every shader below.
const DUST_NOISE_GLSL = `
  float dustHash(vec3 p) {
    p = fract(p * 0.3183099 + vec3(0.1, 0.2, 0.3));
    p *= 17.0;
    return fract(p.x * p.y * p.z * (p.x + p.y + p.z));
  }
  float dustNoise(vec3 x) {
    vec3 i = floor(x), f = fract(x);
    f = f * f * (3.0 - 2.0 * f);
    return mix(mix(mix(dustHash(i), dustHash(i + vec3(1.0, 0.0, 0.0)), f.x),
                   mix(dustHash(i + vec3(0.0, 1.0, 0.0)), dustHash(i + vec3(1.0, 1.0, 0.0)), f.x), f.y),
               mix(mix(dustHash(i + vec3(0.0, 0.0, 1.0)), dustHash(i + vec3(1.0, 0.0, 1.0)), f.x),
                   mix(dustHash(i + vec3(0.0, 1.0, 1.0)), dustHash(i + vec3(1.0, 1.0, 1.0)), f.x), f.y), f.z);
  }
  #ifndef DUST_OCTAVES
  #define DUST_OCTAVES 3
  #endif
  float dustFbm(vec3 p) {
    float a = 0.5, v = 0.0, norm = 0.0;
    for (int k = 0; k < DUST_OCTAVES; k++) { v += a * dustNoise(p); norm += a; p *= 2.03; a *= 0.5; }
    return v / norm;
  }
`;

// Henyey-Greenstein, normalized so that isotropic scattering (g = 0) is 1 and
// `uPhaseNorm` (the reciprocal of the 90-degree value) puts a side-lit sheet at
// the module's nominal brightness. `viewDir` points from the fragment to the
// camera; `uSunWorld` points from the scene toward the sun, so the light
// travels along -uSunWorld and the forward peak is where the two agree.
const DUST_PHASE_GLSL = `
  uniform vec3 uSunWorld;
  uniform float uPhaseG;
  uniform float uPhaseNorm;
  uniform float uBrightness;
  float dustPhase(vec3 viewDir) {
    float g = uPhaseG;
    float cosTheta = dot(-uSunWorld, normalize(viewDir));
    float d = 1.0 + g * g - 2.0 * g * cosTheta;
    float p = (1.0 - g * g) / max(1e-4, d * sqrt(d));
    return clamp(p * uPhaseNorm, ${DUST_PHASE_MIN.toFixed(2)}, ${DUST_PHASE_MAX.toFixed(2)});
  }
`;

// The 90-degree scattering value, used to normalize the phase function so a
// side-lit sheet renders at the nominal brightness (and so a page with no
// lighting handle, whose uSunWorld stays zero, renders at exactly that).
function dustPhaseNorm(g) {
  const d = 1 + g * g;
  return Math.pow(d, 1.5) / Math.max(1e-6, 1 - g * g);
}

const DUST_VERTEX_SHADER = `
  #include <common>
  #include <logdepthbuf_pars_vertex>
  attribute vec4 aSeed;
  uniform float uTime;
  uniform float uLifetime;
  uniform float uSpeed;      // km/s
  uniform float uGravity;    // km/s^2
  uniform float uSizeKm;
  uniform float uProjScale;
  uniform float uFadeKm;
  uniform float uElevMin;
  uniform float uElevMax;
  uniform float uLiftKm;
  varying float vAlpha;
  varying vec3 vView;
  void main() {
    float age = fract(uTime / uLifetime + aSeed.w) * uLifetime;
    float ang = aSeed.x * 6.28318530718;
    float elev = mix(uElevMin, uElevMax, aSeed.z);
    // Most of the eroded mass creeps out slowly and a thin tail streaks away at
    // the full ejecta speed, so the sheet stays dense where the soil leaves it.
    float v = uSpeed * (0.05 + 0.95 * aSeed.y * aSeed.y);
    vec3 dir = vec3(cos(ang) * cos(elev), sin(ang) * cos(elev), sin(elev));
    vec3 p = dir * (v * age);
    p.z = max(p.z - 0.5 * uGravity * age * age, 0.0) + uLiftKm;
    float radial = length(p.xy);
    vAlpha = (1.0 - age / uLifetime) * exp(-radial / uFadeKm);
    vec4 world = modelMatrix * vec4(p, 1.0);
    vView = cameraPosition - world.xyz;
    vec4 mv = modelViewMatrix * vec4(p, 1.0);
    gl_Position = projectionMatrix * mv;
    #include <logdepthbuf_vertex>
    gl_PointSize = clamp(uSizeKm * uProjScale / max(1e-7, -mv.z), ${DUST_POINT_MIN_PX.toFixed(1)}, ${DUST_POINT_MAX_PX.toFixed(1)});
  }
`;

const DUST_FRAGMENT_SHADER = `
  #include <common>
  #include <logdepthbuf_pars_fragment>
  // (the tone mapping and color space helpers are already in three's fragment prefix)
  ${DUST_PHASE_GLSL}
  uniform vec3 uColor;
  uniform float uOpacity;
  varying float vAlpha;
  varying vec3 vView;
  void main() {
    #include <logdepthbuf_fragment>
    vec2 d = gl_PointCoord - vec2(0.5);
    float r2 = dot(d, d);
    if (r2 > 0.25) discard;
    float soft = 1.0 - smoothstep(0.0, 0.25, r2);
    gl_FragColor = vec4(uColor * uBrightness * dustPhase(vView), 1.0);
    #include <tonemapping_fragment>
    #include <colorspace_fragment>
    gl_FragColor.a = vAlpha * uOpacity * soft;
  }
`;

// Per-particle seeds: azimuth, speed fraction, elevation fraction, emission phase.
function dustSeedGeometry(count) {
  const seeds = new Float32Array(4 * count);
  let seed = 0x9e3779b9;
  const rand = () => {
    seed = (seed * 1664525 + 1013904223) >>> 0;
    return seed / 4294967296;
  };
  for (let i = 0; i < count; i++) {
    seeds[4 * i] = rand();
    seeds[4 * i + 1] = rand();
    seeds[4 * i + 2] = rand();
    seeds[4 * i + 3] = rand();
  }
  const geometry = new THREE.BufferGeometry();
  geometry.setAttribute('position', new THREE.BufferAttribute(new Float32Array(3 * count), 3));
  geometry.setAttribute('aSeed', new THREE.BufferAttribute(seeds, 4));
  geometry.boundingSphere = new THREE.Sphere(new THREE.Vector3(0, 0, 0), 1e3);
  return geometry;
}

// Sun uniforms every dust material carries; `dustSetSun` fills them in.
function dustSunUniforms(g) {
  return {
    uSunWorld: { value: new THREE.Vector3(0, 0, 0) },
    uPhaseG: { value: g },
    uPhaseNorm: { value: dustPhaseNorm(g) },
    uBrightness: { value: 1 },
  };
}

function dustMaterial(gravityKmS2, projScale) {
  return new THREE.ShaderMaterial({
    uniforms: {
      uTime: { value: 0 },
      uLifetime: { value: DUST_LIFETIME_S },
      uSpeed: { value: 0 },
      uGravity: { value: gravityKmS2 },
      uSizeKm: { value: DUST_SIZE_KM },
      uProjScale: { value: projScale },
      uFadeKm: { value: 0.05 },
      uElevMin: { value: DUST_ELEVATION_MIN_RAD },
      uElevMax: { value: DUST_ELEVATION_MAX_RAD },
      uLiftKm: { value: DUST_SHEET_LIFT_M * DUST_M_TO_KM },
      uColor: { value: new THREE.Color(DUST_COLOR) },
      uOpacity: { value: 0 },
      ...dustSunUniforms(DUST_PHASE_G),
    },
    vertexShader: DUST_VERTEX_SHADER,
    fragmentShader: DUST_FRAGMENT_SHADER,
    transparent: true,
    depthWrite: false,
    blending: THREE.AdditiveBlending,
  });
}

// The same particles, written into the shadow map. three r160 does render
// Points into a shadow pass (its shadow map walks isMesh, isLine and isPoints
// alike) and it does take an object's `customDepthMaterial`, which is the only
// way to get these particles in: the built-in depth material knows nothing of
// the vertex shader that places them, so it would pack the depth of the
// emitter's origin for every one of them. Sharing the display material's
// uniform object keeps the two passes on the same particle positions.
//
// It is off by default all the same. Measured on the Apollo page it changes
// almost nothing: the vertex shader sizes a particle for the perspective
// camera, so in the shadow camera's orthographic pass each one lands on a
// texel or two, and turning it on moves a couple of hundred pixels of the
// frame. Sized up it would be worse, not better -- a few hundred opaque disks
// speckle the ground where a dust cloud casts one soft darkening, which is
// what the decal below draws. `pointShadows: true` turns them on anyway.
const DUST_DEPTH_FRAGMENT_SHADER = `
  #include <common>
  #include <packing>
  #include <logdepthbuf_pars_fragment>
  varying vec2 vHighPrecisionZW;
  void main() {
    #include <logdepthbuf_fragment>
    vec2 d = gl_PointCoord - vec2(0.5);
    if (dot(d, d) > 0.25) discard;
    float fragCoordZ = 0.5 * vHighPrecisionZW[0] / vHighPrecisionZW[1] + 0.5;
    gl_FragColor = packDepthToRGBA(fragCoordZ);
  }
`;

function dustDepthMaterial(uniforms) {
  return new THREE.ShaderMaterial({
    uniforms,
    vertexShader: DUST_VERTEX_SHADER.replace('varying vec3 vView;', 'varying vec3 vView;\n  varying vec2 vHighPrecisionZW;')
      .replace('#include <logdepthbuf_vertex>', '#include <logdepthbuf_vertex>\n    vHighPrecisionZW = gl_Position.zw;'),
    fragmentShader: DUST_DEPTH_FRAGMENT_SHADER,
  });
}

// ---- sheet ----------------------------------------------------------------

const DUST_SHEET_VERTEX_SHADER = `
  #include <common>
  #include <logdepthbuf_pars_vertex>
  varying vec2 vSheetUv;
  varying vec3 vSheetView;
  varying vec3 vSheetNormal;
  void main() {
    vSheetUv = (uv - 0.5) * 2.0;
    vSheetNormal = normalize(normalMatrix * normal);
    vec4 world = modelMatrix * vec4(position, 1.0);
    vSheetView = cameraPosition - world.xyz;
    gl_Position = projectionMatrix * viewMatrix * world;
    #include <logdepthbuf_vertex>
  }
`;

// Advected turbulence in polar coordinates: the noise is sampled on a circle
// of radius uSpokes (so it is seamless across +/- pi and has roughly
// 2*pi*uSpokes cells around, which is what draws the radial streaks) crossed
// with a radial coordinate that scrolls outward at the ejecta speed. The layer
// is treated as a slab: a thin sheet seen edge-on is a long way through, which
// is why it reads as a low wall from a grazing camera and as a veil from
// overhead, and why the disks have no visible rim at a shallow angle. The
// output is premultiplied -- the color carries the light the dust scatters and
// the alpha carries how much of the scene behind it is hidden.
const DUST_SHEET_FRAGMENT_SHADER = `
  #include <common>
  #include <logdepthbuf_pars_fragment>
  // (the tone mapping and color space helpers are already in three's fragment prefix)
  ${DUST_NOISE_GLSL}
  ${DUST_PHASE_GLSL}
  uniform vec3 uColor;
  uniform float uTau;
  uniform float uExtinction;
  uniform float uTime;
  uniform float uSeed;
  uniform float uScroll;
  uniform float uSwirl;
  uniform float uSpokes;
  uniform float uRadialFreq;
  uniform float uSharpness;
  uniform float uDecay;
  uniform float uFloor;
  uniform float uInnerHole;
  varying vec2 vSheetUv;
  varying vec3 vSheetView;
  varying vec3 vSheetNormal;
  void main() {
    #include <logdepthbuf_fragment>
    float r = length(vSheetUv);
    if (r > 1.0) discard;
    float ang = atan(vSheetUv.y, vSheetUv.x) + uSwirl * r;
    vec3 q = vec3(cos(ang), sin(ang), 0.0) * uSpokes
           + vec3(0.0, 0.0, r * uRadialFreq - uTime * uScroll)
           + uSeed;
    float n = dustFbm(q);
    #ifdef DUST_FINE
    // A finer pass, mostly across the streaks, breaks the base cells into the
    // billowing filaments the overhead test footage shows.
    n *= 0.55 + 0.75 * dustFbm(q * vec3(2.7, 2.7, 1.9) + 11.0);
    #endif
    float turb = pow(clamp(n, 0.0, 1.0), uSharpness);
    float rim = smoothstep(1.0, 0.85, r);
    float decay = uFloor + (1.0 - uFloor) * exp(-uDecay * r);
    float inner = smoothstep(0.0, uInnerHole, r);
    float facing = max(abs(dot(normalize(vSheetNormal), normalize(vSheetView))), 0.10);
    float tau = uTau * turb * rim * decay * inner / facing;
    float a = 1.0 - exp(-tau);
    if (a <= 0.004) discard;
    gl_FragColor = vec4(uColor * (uBrightness * dustPhase(vSheetView)), 1.0);
    #include <tonemapping_fragment>
    #include <colorspace_fragment>
    gl_FragColor = vec4(gl_FragColor.rgb * a, a * uExtinction);
  }
`;

// Premultiplied "over": the sheet adds the light it scatters and hides what is
// behind it in the same pass, which is what a translucent cloud does.
function dustCloudBlending(material) {
  material.transparent = true;
  material.depthWrite = false;
  material.blending = THREE.CustomBlending;
  material.blendEquation = THREE.AddEquation;
  material.blendSrc = THREE.OneFactor;
  material.blendDst = THREE.OneMinusSrcAlphaFactor;
  material.blendSrcAlpha = THREE.OneFactor;
  material.blendDstAlpha = THREE.OneMinusSrcAlphaFactor;
  return material;
}

// One disk of a stack. `style` carries everything the two stacks differ in.
function dustStackLayer(spec, index, style) {
  const material = dustCloudBlending(new THREE.ShaderMaterial({
    uniforms: {
      uColor: { value: new THREE.Color(style.color) },
      uTau: { value: 0 },
      uExtinction: { value: style.extinction },
      uTime: { value: 0 },
      uSeed: { value: spec.seed },
      uScroll: { value: 0 },
      uSwirl: { value: spec.swirl },
      uSpokes: { value: style.spokes },
      uRadialFreq: { value: style.radialFreq },
      uSharpness: { value: style.sharpness },
      uDecay: { value: style.decay },
      uFloor: { value: style.floor },
      uInnerHole: { value: 0.08 },
      ...dustSunUniforms(style.phaseG),
    },
    defines: style.fine ? { DUST_FINE: '' } : { DUST_OCTAVES: 2 },
    vertexShader: DUST_SHEET_VERTEX_SHADER,
    fragmentShader: DUST_SHEET_FRAGMENT_SHADER,
    side: THREE.DoubleSide,
  }));
  const mesh = new THREE.Mesh(new THREE.RingGeometry(0, 1, DUST_SHEET_SEGMENTS, DUST_SHEET_RINGS), material);
  mesh.renderOrder = style.renderOrder + index;
  mesh.frustumCulled = false;
  mesh.visible = false;
  mesh.userData.weight = spec.weight;
  mesh.userData.heightM = spec.height_m;
  mesh.userData.scroll = spec.scroll;
  return mesh;
}

// The ground-hugging sheet, and the haze that builds up above it.
const DUST_SHEET_STYLE = {
  color: DUST_COLOR, extinction: DUST_SHEET_EXTINCTION, spokes: DUST_SHEET_SPOKES,
  radialFreq: DUST_SHEET_RADIAL_FREQ, sharpness: DUST_SHEET_SHARPNESS, decay: DUST_SHEET_DECAY,
  floor: DUST_SHEET_FLOOR, phaseG: DUST_PHASE_G, fine: true, renderOrder: 26,
};

// ---- haze puffs -----------------------------------------------------------

// One quad per puff, turned to face the camera about the local vertical, so a
// puff is wide and low from the side and a broad ellipse from overhead -- the
// projection of the flattened ellipsoid of dust it stands for. `aCenter` is in
// units of the bank's radius (z in units of its height) and `aParam` carries
// the puff's own radius, seed and drift phase.
const DUST_PUFF_VERTEX_SHADER = `
  #include <common>
  #include <logdepthbuf_pars_vertex>
  attribute vec3 aCenter;
  attribute vec3 aParam;      // radius fraction, seed, drift phase
  uniform vec3 uScale;        // bank radius, radius, height (km)
  uniform float uPuffHeight;  // km
  uniform float uTime;
  uniform float uDrift;
  varying vec2 vPuffUv;
  varying vec3 vPuffView;
  varying float vPuffSeed;
  varying float vPuffFade;
  void main() {
    // Each puff drifts out from the middle of the bank and is replaced, which
    // keeps the haze churning instead of hanging as a fixed blob.
    float phase = fract(uTime * uDrift + aParam.z);
    float spread = mix(0.30, 1.0, phase);
    vPuffFade = sin(3.14159265 * phase);
    vec3 local = vec3(aCenter.xy * spread, aCenter.z) * uScale;
    vec4 world = modelMatrix * vec4(local, 1.0);
    vPuffView = cameraPosition - world.xyz;
    vec3 upWorld = normalize((modelMatrix * vec4(0.0, 0.0, 1.0, 0.0)).xyz);
    vec4 c = viewMatrix * world;
    vec3 upView = normalize((viewMatrix * vec4(upWorld, 0.0)).xyz);
    vec3 toCam = normalize(-c.xyz);
    vec3 right = cross(upView, toCam);
    float axial = abs(dot(upView, toCam));
    right = length(right) > 1e-6 ? normalize(right) : vec3(1.0, 0.0, 0.0);
    vec3 up = normalize(cross(toCam, right));
    float wide = aParam.x * uScale.x;
    // Seen along its own axis the flattened puff shows its full width; seen
    // from the side it shows its height.
    float tall = mix(uPuffHeight, wide, axial);
    vec3 p = c.xyz + right * (position.x * wide) + up * (position.y * tall);
    // A puff the camera is nearly inside would read as one soft ball of fog
    // filling the frame, so it fades out as the camera closes on it and the
    // ones behind it carry the haze instead.
    vPuffFade *= smoothstep(0.7, 2.2, length(c.xyz) / max(1e-9, wide));
    vPuffUv = position.xy;
    vPuffSeed = aParam.y;
    gl_Position = projectionMatrix * vec4(p, 1.0);
    #include <logdepthbuf_vertex>
  }
`;

const DUST_PUFF_FRAGMENT_SHADER = `
  #include <common>
  #include <logdepthbuf_pars_fragment>
  // (the tone mapping and color space helpers are already in three's fragment prefix)
  ${DUST_NOISE_GLSL}
  ${DUST_PHASE_GLSL}
  uniform vec3 uColor;
  uniform float uTau;
  uniform float uTime;
  varying vec2 vPuffUv;
  varying vec3 vPuffView;
  varying float vPuffSeed;
  varying float vPuffFade;
  void main() {
    #include <logdepthbuf_fragment>
    float r = length(vPuffUv);
    if (r > 1.0) discard;
    float falloff = 1.0 - smoothstep(0.0, 1.0, r);
    // Torn edges and lanes inside: a puff close to the camera has to read as
    // part of a cloud, not as a sphere of fog with a rim.
    float n = dustFbm(vec3(vPuffUv * 2.6, uTime * 0.05) + vPuffSeed);
    n *= 0.6 + 0.8 * dustFbm(vec3(vPuffUv * 6.1, uTime * 0.09) + vPuffSeed * 1.7);
    n = clamp(n * 2.1, 0.0, 1.4);
    float a = 1.0 - exp(-uTau * falloff * n * vPuffFade);
    if (a <= 0.004) discard;
    gl_FragColor = vec4(uColor * (uBrightness * dustPhase(vPuffView)), 1.0);
    #include <tonemapping_fragment>
    #include <colorspace_fragment>
    gl_FragColor = vec4(gl_FragColor.rgb * a, a);
  }
`;

// The bank: `count` quads in one buffer, laid out in a disk around the
// impingement point and jittered in height, size and phase.
function dustPuffGeometry(count) {
  const corners = [[-1, -1], [1, -1], [1, 1], [-1, 1]];
  const position = new Float32Array(12 * count);
  const center = new Float32Array(12 * count);
  const param = new Float32Array(12 * count);
  const index = new Uint16Array(6 * count);
  let seed = 0x2545f491;
  const rand = () => { seed = (seed * 1664525 + 1013904223) >>> 0; return seed / 4294967296; };
  for (let i = 0; i < count; i++) {
    const theta = 2 * Math.PI * rand();
    const radius = 0.72 * Math.sqrt(rand());
    const cx = radius * Math.cos(theta), cy = radius * Math.sin(theta), cz = 0.15 + 0.85 * rand() * rand();
    const size = DUST_HAZE_PUFF_MIN + (DUST_HAZE_PUFF_MAX - DUST_HAZE_PUFF_MIN) * rand();
    const noiseSeed = 37.0 * rand(), phase = rand();
    for (let k = 0; k < 4; k++) {
      const o = 12 * i + 3 * k;
      position[o] = corners[k][0]; position[o + 1] = corners[k][1]; position[o + 2] = 0;
      center[o] = cx; center[o + 1] = cy; center[o + 2] = cz;
      param[o] = size; param[o + 1] = noiseSeed; param[o + 2] = phase;
    }
    const v = 4 * i, e = 6 * i;
    index[e] = v; index[e + 1] = v + 1; index[e + 2] = v + 2;
    index[e + 3] = v; index[e + 4] = v + 2; index[e + 5] = v + 3;
  }
  const geometry = new THREE.BufferGeometry();
  geometry.setAttribute('position', new THREE.BufferAttribute(position, 3));
  geometry.setAttribute('aCenter', new THREE.BufferAttribute(center, 3));
  geometry.setAttribute('aParam', new THREE.BufferAttribute(param, 3));
  geometry.setIndex(new THREE.BufferAttribute(index, 1));
  geometry.boundingSphere = new THREE.Sphere(new THREE.Vector3(0, 0, 0), 1e3);
  return geometry;
}

function dustHazeBank() {
  const material = dustCloudBlending(new THREE.ShaderMaterial({
    uniforms: {
      uColor: { value: new THREE.Color(DUST_HAZE_COLOR) },
      uTau: { value: 0 },
      uTime: { value: 0 },
      uScale: { value: new THREE.Vector3(0, 0, 0) },
      uPuffHeight: { value: DUST_HAZE_PUFF_HEIGHT_M * DUST_M_TO_KM },
      uDrift: { value: DUST_HAZE_DRIFT_HZ },
      ...dustSunUniforms(DUST_HAZE_PHASE_G),
    },
    defines: { DUST_OCTAVES: 2 },
    vertexShader: DUST_PUFF_VERTEX_SHADER,
    fragmentShader: DUST_PUFF_FRAGMENT_SHADER,
    side: THREE.DoubleSide,
  }));
  const mesh = new THREE.Mesh(dustPuffGeometry(DUST_HAZE_PUFFS), material);
  mesh.renderOrder = 32;
  mesh.frustumCulled = false;
  mesh.visible = false;
  return mesh;
}

// ---- decals ---------------------------------------------------------------

const DUST_DECAL_VERTEX_SHADER = `
  #include <common>
  #include <logdepthbuf_pars_vertex>
  varying vec2 vDecalUv;
  void main() {
    vDecalUv = uv;
    gl_Position = projectionMatrix * modelViewMatrix * vec4(position, 1.0);
    #include <logdepthbuf_vertex>
  }
`;

// A flat ellipse lying on the ground with a soft rim: the sheet's own shadow,
// stretched away from the sun. Drawn as ordinary alpha over a near-black
// color, which is a multiply: it darkens the ground it lies on whatever the
// ground's own brightness is, instead of painting a fixed gray over it.
const DUST_DECAL_FRAGMENT_SHADER = `
  #include <common>
  #include <logdepthbuf_pars_fragment>
  uniform vec3 uColor;
  uniform float uOpacity;
  uniform float uCore;
  varying vec2 vDecalUv;
  void main() {
    #include <logdepthbuf_fragment>
    float r = length(vDecalUv - vec2(0.5)) * 2.0;
    float a = uOpacity * (1.0 - smoothstep(uCore, 1.0, r));
    if (a <= 0.003) discard;
    gl_FragColor = vec4(uColor, a);
  }
`;

// The scour crater: a floor the plume has swept smooth and bright inside a
// ring of the darker soil pushed out of it. The two are separate meshes
// because they compose differently -- the rim darkens what is under it, the
// floor adds light -- and because an added brightness has to be scaled by the
// display gain, or it burns to white on a page exposed for a grazing sun.
const DUST_SCOUR_FRAGMENT_SHADER = `
  #include <common>
  #include <logdepthbuf_pars_fragment>
  uniform vec3 uColor;
  uniform float uOpacity;
  uniform float uRim;
  uniform float uRimWidth;
  uniform float uPickRim;
  varying vec2 vDecalUv;
  void main() {
    #include <logdepthbuf_fragment>
    float r = length(vDecalUv - vec2(0.5)) * 2.0;
    if (r > 1.0) discard;
    float floorMask = 1.0 - smoothstep(uRim - 0.30, uRim, r);
    float d = (r - uRim) / uRimWidth;
    float rimMask = exp(-d * d);
    float a = uOpacity * mix(floorMask, rimMask, uPickRim) * (1.0 - smoothstep(0.92, 1.0, r));
    if (a <= 0.003) discard;
    gl_FragColor = vec4(uColor, a);
  }
`;

function dustDecal(color, opacity, renderOrder, core) {
  const material = new THREE.ShaderMaterial({
    uniforms: {
      uColor: { value: new THREE.Color(color) },
      uOpacity: { value: 0 },
      uCore: { value: core },
    },
    vertexShader: DUST_DECAL_VERTEX_SHADER,
    fragmentShader: DUST_DECAL_FRAGMENT_SHADER,
    transparent: true,
    depthWrite: false,
    side: THREE.DoubleSide,
    blending: THREE.NormalBlending,
  });
  const mesh = new THREE.Mesh(new THREE.CircleGeometry(1, 64), material);
  mesh.renderOrder = renderOrder;
  mesh.visible = false;
  mesh.userData.maxOpacity = opacity;
  return mesh;
}

function dustScourPart(color, opacity, renderOrder, pickRim, additive) {
  const material = new THREE.ShaderMaterial({
    uniforms: {
      uColor: { value: new THREE.Color(color) },
      uOpacity: { value: 0 },
      uRim: { value: DUST_SCOUR_RIM_RADIUS },
      uRimWidth: { value: DUST_SCOUR_RIM_WIDTH },
      uPickRim: { value: pickRim },
    },
    vertexShader: DUST_DECAL_VERTEX_SHADER,
    fragmentShader: DUST_SCOUR_FRAGMENT_SHADER,
    transparent: true,
    depthWrite: false,
    side: THREE.DoubleSide,
    blending: additive ? THREE.AdditiveBlending : THREE.NormalBlending,
  });
  const mesh = new THREE.Mesh(new THREE.CircleGeometry(1, 96), material);
  mesh.renderOrder = renderOrder;
  mesh.visible = false;
  mesh.userData.maxOpacity = opacity;
  return mesh;
}

function dustMax(array) {
  let m = 0;
  if (!array) return m;
  for (let i = 0; i < array.length; i++) if (array[i] > m) m = array[i];
  return m;
}

// The lighting handle. main.js passes a function because the lighting module is
// built after this one (it has to light the bodies that already exist); a plain
// handle, or nothing at all, works just as well.
function dustLighting(options) {
  const l = options.lighting;
  return typeof l === 'function' ? l() : (l || null);
}

// How bright the scattered light is, in the units the scene is lit in. The
// dust shaders tone map and encode their output exactly as every other
// material does, so the camera's exposure is applied once, downstream of this,
// and what is wanted here is a radiance: the dust's scattering albedo times
// the sun's strength in the scene times its irradiance at the body.
//
// A page whose lighting module has not raised the exposure is not working in
// radiance at all -- it renders display values directly, the way the viewer
// did before there was a lighting module -- and there the sheet keeps the
// fixed brightness it always had.
function dustDisplayScale(lighting) {
  if (!lighting) return DUST_FALLBACK_BRIGHTNESS;
  const exposure = Number.isFinite(lighting.exposure) && lighting.exposure > 0 ? lighting.exposure : 1;
  if (!(exposure > 1.01)) return DUST_FALLBACK_BRIGHTNESS;
  const irradiance = Number.isFinite(lighting.sunIrradiance) && lighting.sunIrradiance > 0
    ? lighting.sunIrradiance : DUST_REFERENCE_IRRADIANCE_W_M2;
  const intensity = lighting.sun && Number.isFinite(lighting.sun.intensity) ? lighting.sun.intensity : 1;
  const scale = DUST_SCATTER_ALBEDO * intensity * (irradiance / DUST_REFERENCE_IRRADIANCE_W_M2);
  return Math.min(DUST_SCALE_MAX, Math.max(DUST_SCALE_MIN, scale));
}

// An exponentially decaying integral of the erosion rate, one value per frame
// and spacecraft, normalized to its own peak: how much fine dust is still in
// the air. Precomputed over the whole timeline so scrubbing backwards shows
// the same haze as playing forwards.
function dustHazeTrack(frames, tauS) {
  const S = frames.sats, N = frames.count;
  const rate = frames.plume.erosion_kg_s;
  const track = new Float32Array(N * S);
  let peak = 0;
  for (let s = 0; s < S; s++) {
    let acc = 0;
    for (let k = 0; k < N; k++) {
      const dt = k === 0 ? 0 : Math.max(0, frames.t[k] - frames.t[k - 1]);
      acc = acc * Math.exp(-dt / tauS) + (rate[k * S + s] || 0) * dt;
      track[k * S + s] = acc;
      if (acc > peak) peak = acc;
    }
  }
  if (peak > 0) for (let i = 0; i < track.length; i++) track[i] /= peak;
  return track;
}

// Elapsed time of the first frame with any erosion, per spacecraft (Infinity
// when a spacecraft never disturbs the ground): the sheet's radius grows from
// there.
function dustOnsetTimes(frames) {
  const S = frames.sats, N = frames.count;
  const rate = frames.plume.erosion_kg_s;
  const onset = new Float64Array(S).fill(Infinity);
  for (let s = 0; s < S; s++) {
    for (let k = 0; k < N; k++) {
      if (rate[k * S + s] > 0) { onset[s] = frames.t[k]; break; }
    }
  }
  return onset;
}

/**
 * Dust blown off the surface by a descent engine.
 *
 * `parent` is the object the dust hangs from — the viewer passes its inertial
 * `world` group, so the sheet rides the same frame as the spacecraft positions.
 * `terrain` is the site terrain (used through `heightAt` to lay the decals on
 * the DEM surface) and `lod` supplies the vehicle attitude that orients the
 * engine axis. Returns `{ update(t), setVisible(v), visible, group }`.
 *
 * Options: `rotationAt(t, out)` (the planet rotation, for the terrain lookup),
 * `lighting` (the lighting handle, or a function returning it — read every
 * update for the sun direction, its irradiance and the exposure; without it
 * the sheet keeps a fixed, unlit brightness), `gravity_m_s2`,
 * `projectionScale`, `maxEmitters`, `hazeOpacity`, `hazeDecayS`,
 * `sheetMaxRadiusM`, `pointShadows`.
 */
export function createDust(parent, frames, terrain, lod, options = {}) {
  const group = new THREE.Group();
  group.name = 'dust';
  const empty = { group, update() {}, setVisible() {}, get visible() { return false; } };
  if (!frames || !frames.hasPlume || !frames.hasPlume()) return empty;

  const gravityKmS2 = (options.gravity_m_s2 ?? DUST_GRAVITY_M_S2) * DUST_M_TO_KM;
  const projScale = options.projectionScale ?? DUST_PROJECTION_SCALE;
  const rotationAt = typeof options.rotationAt === 'function' ? options.rotationAt : null;
  const referenceRadiusKm = terrain && terrain.referenceRadiusKm ? terrain.referenceRadiusKm : 0;
  const heightAt = terrain && typeof terrain.heightAt === 'function' ? terrain.heightAt : null;
  const erosionMax = dustMax(frames.plume.erosion_kg_s) || 1;
  const erodedMax = dustMax(frames.plume.eroded_kg) || 1;
  const count = Math.min(frames.sats, options.maxEmitters ?? DUST_MAX_EMITTERS);
  const hazeOpacity = options.hazeOpacity ?? DUST_HAZE_OPACITY;
  const hazeTrack = dustHazeTrack(frames, Math.max(0.5, options.hazeDecayS ?? DUST_HAZE_DECAY_S));
  const onsetTimes = dustOnsetTimes(frames);
  const sheetMaxKm = (options.sheetMaxRadiusM ?? DUST_SHEET_MAX_RADIUS_M) * DUST_M_TO_KM;

  const emitters = [];
  for (let s = 0; s < count; s++) {
    const pivot = new THREE.Group();
    const material = dustMaterial(gravityKmS2, projScale);
    const points = new THREE.Points(dustSeedGeometry(DUST_PARTICLE_COUNT), material);
    points.frustumCulled = false;
    points.renderOrder = 30;
    points.customDepthMaterial = dustDepthMaterial(material.uniforms);
    points.castShadow = !!options.pointShadows;
    const sheet = DUST_SHEET_LAYERS.map((spec, k) => dustStackLayer(spec, k, DUST_SHEET_STYLE));
    const haze = dustHazeBank();
    const scourFloor = dustScourPart(DUST_SCOUR_FLOOR_COLOR, DUST_SCOUR_FLOOR_OPACITY, 21, 0, true);
    const scourRim = dustScourPart(DUST_SCOUR_RIM_COLOR, DUST_SCOUR_RIM_OPACITY, 20, 1, false);
    const shadow = dustDecal(DUST_SHADOW_COLOR, DUST_SHADOW_OPACITY, 22, 0.25);
    pivot.add(scourRim, scourFloor, shadow, ...sheet, points, haze);
    pivot.visible = false;
    group.add(pivot);
    emitters.push({ s, pivot, points, material, sheet, haze, scourFloor, scourRim, shadow, placed: false });
  }
  if (parent && typeof parent.add === 'function') parent.add(group);

  const posKm = new Float64Array(3), qAtt = new Float32Array(4), qPlanet = new Float32Array(4);
  const bodyKm = new Float64Array(3);
  const axis = new THREE.Vector3(), impact = new THREE.Vector3(), up = new THREE.Vector3();
  const quat = new THREE.Quaternion();
  const sunScene = new THREE.Vector3(), sunLocal = new THREE.Vector3(), pivotQuat = new THREE.Quaternion();
  let visible = true;

  // Ground point under the engine: march from the vehicle along the engine axis
  // by the height the effector measured, then snap onto the DEM when the page
  // carries one (the sheet is planar above that point).
  function dustImpactPoint(t, s, height) {
    frames.positionAt(t, s, posKm);
    if (!Number.isFinite(posKm[0])) return false;
    if (lod && typeof lod.attitudeAt === 'function') {
      lod.attitudeAt(t, s, qAtt);
      quat.set(qAtt[0], qAtt[1], qAtt[2], qAtt[3]);
      axis.set(0, 0, 1).applyQuaternion(quat);           // the engine fires along body +z
    } else {
      axis.set(posKm[0], posKm[1], posKm[2]).normalize().negate();
    }
    impact.set(posKm[0], posKm[1], posKm[2]).addScaledVector(axis, height * DUST_M_TO_KM);
    up.copy(impact).normalize();
    if (heightAt && rotationAt && referenceRadiusKm > 0) {
      rotationAt(t, qPlanet);
      bodyKm[0] = impact.x; bodyKm[1] = impact.y; bodyKm[2] = impact.z;
      rotateByConjugate(qPlanet, bodyKm, bodyKm);
      const r = Math.hypot(bodyKm[0], bodyKm[1], bodyKm[2]);
      const lat = THREE.MathUtils.radToDeg(Math.asin(bodyKm[2] / r));
      const lon = THREE.MathUtils.radToDeg(Math.atan2(bodyKm[1], bodyKm[0]));
      const h = heightAt(lat, lon);
      if (Number.isFinite(h)) impact.copy(up).multiplyScalar(referenceRadiusKm + h * DUST_M_TO_KM);
    }
    return true;
  }

  // Sun direction for this frame, in scene space, and whether there is one.
  function dustSunDirection(lighting) {
    if (!lighting || typeof lighting.direction !== 'function') return false;
    lighting.direction(sunScene);
    return sunScene.lengthSq() > 1e-12;
  }

  function dustApplySun(material, hasSun, brightness) {
    const u = material.uniforms;
    if (hasSun) u.uSunWorld.value.copy(sunScene); else u.uSunWorld.value.set(0, 0, 0);
    u.uBrightness.value = brightness;
  }

  // The sheet's shadow: an ellipse on the ground, offset and stretched away
  // from the sun by the sheet's own height over the surface. The sheet is only
  // meters thick and tens of meters wide, so at a grazing sun the shadow is a
  // long soft smear, and at a high sun it collapses onto the sheet itself.
  function dustUpdateShadow(item, hasSun, sheetRadiusKm, thicknessKm, strength) {
    if (!hasSun) { item.shadow.visible = false; return; }
    item.pivot.getWorldQuaternion(pivotQuat).invert();
    sunLocal.copy(sunScene).applyQuaternion(pivotQuat);
    const horizontal = Math.hypot(sunLocal.x, sunLocal.y);
    const elevation = Math.atan2(sunLocal.z, horizontal);
    if (!(elevation > DUST_SHADOW_MIN_ELEVATION_RAD) || horizontal < 1e-6) { item.shadow.visible = false; return; }
    const runKm = Math.min(DUST_SHADOW_MAX_STRETCH * sheetRadiusKm, thicknessKm / Math.tan(elevation));
    const halfLength = sheetRadiusKm + 0.5 * runKm;
    const azimuth = Math.atan2(-sunLocal.y, -sunLocal.x);   // away from the sun
    item.shadow.visible = true;
    item.shadow.rotation.z = azimuth;
    item.shadow.position.set(
      0.5 * runKm * Math.cos(azimuth),
      0.5 * runKm * Math.sin(azimuth),
      0.45 * DUST_DECAL_LIFT_M * DUST_M_TO_KM,
    );
    item.shadow.scale.set(halfLength, sheetRadiusKm, 1);
    item.shadow.material.uniforms.uOpacity.value = item.shadow.userData.maxOpacity * strength;
  }

  function dustUpdateEmitter(item, t, hasSun, displayScale) {
    const s = item.s;
    const erosion = frames.plumeAt('erosion_kg_s', t, s);
    const haze = dustHazeAt(t, s);
    const active = erosion > 0;
    if (!active && !(haze > 0.005)) { item.pivot.visible = false; return; }
    const height = frames.plumeAt('height_m', t, s);
    // Once the engine is off there may be no usable height any more, but the
    // haze is still hanging over the spot: keep the last impingement point.
    if (Number.isFinite(height) && dustImpactPoint(t, s, height)) {
      item.pivot.position.copy(impact);
      item.pivot.quaternion.setFromUnitVectors(DUST_UP_AXIS, up);
      item.placed = true;
    } else if (!item.placed) {
      item.pivot.visible = false;
      return;
    }
    item.pivot.visible = true;
    item.pivot.updateMatrixWorld();

    const ejecta = frames.plumeAt('ejecta_mps', t, s);
    const eroded = frames.plumeAt('eroded_kg', t, s);
    const strength = active ? Math.min(1, erosion / erosionMax) : 0;
    const speedKmS = Math.max(1e-4, (Number.isFinite(ejecta) ? ejecta : 0) * DUST_M_TO_KM);
    const reachKm = speedKmS * DUST_LIFETIME_S;
    const footprintKm = Math.max(1e-5, (Number.isFinite(height) ? height : 0) * DUST_M_TO_KM * DUST_PLUME_TAN_HALF_ANGLE);
    // The visible front of the sheet leaves the impingement point when erosion
    // starts and runs outward at a fraction of the ejecta speed.
    const elapsed = Math.max(0, t - onsetTimes[s]);
    const sheetRadiusKm = Math.min(sheetMaxKm,
      DUST_SHEET_FOOTPRINTS * footprintKm + DUST_SHEET_SPREAD * speedKmS * elapsed);

    // Fast grains that outrun the sheet: a thin sparkle over it, not the body
    // of the effect.
    const u = item.material.uniforms;
    u.uTime.value = t;
    u.uSpeed.value = speedKmS;
    u.uFadeKm.value = 0.22 * reachKm;
    u.uOpacity.value = DUST_MAX_OPACITY * Math.sqrt(strength);
    dustApplySun(item.material, hasSun, displayScale);
    item.points.visible = strength > 0;
    item.points.geometry.setDrawRange(0, Math.max(16, Math.round(DUST_PARTICLE_COUNT * DUST_EJECTA_FRACTION * strength)));

    // The sheet: a thin stack of advected noise disks. It lingers a little
    // after the engine stops, on the haze that is still settling.
    const sheetStrength = Math.max(strength, 0.45 * haze);
    // Scroll rate in sheet radii per second, so the filaments leave the center
    // at the speed the grains do however wide the sheet has grown.
    const scrollRate = DUST_SHEET_RADIAL_FREQ * DUST_SHEET_SCROLL_FRACTION * speedKmS / Math.max(1e-6, sheetRadiusKm);
    const innerHole = Math.min(0.35, Math.max(0.04, footprintKm / Math.max(1e-6, sheetRadiusKm)));
    for (const layer of item.sheet) {
      layer.visible = sheetStrength > 0.01;
      if (!layer.visible) continue;
      layer.position.set(0, 0, (DUST_SHEET_LIFT_M + layer.userData.heightM) * DUST_M_TO_KM);
      layer.scale.set(sheetRadiusKm, sheetRadiusKm, 1);
      const lu = layer.material.uniforms;
      lu.uTime.value = t;
      lu.uTau.value = DUST_SHEET_TAU * layer.userData.weight * Math.sqrt(sheetStrength);
      lu.uScroll.value = scrollRate * layer.userData.scroll;
      lu.uInnerHole.value = innerHole;
      dustApplySun(layer.material, hasSun, displayScale);
    }

    // Haze: the fines still in the air, over a bank that keeps growing. Its
    // top is what casts the sheet's shadow on the ground.
    const hazeRadiusKm = Math.min((options.hazeRadiusM ?? DUST_HAZE_RADIUS_M) * DUST_M_TO_KM,
      sheetRadiusKm * (1.5 + 1.5 * haze));
    const hazeTopKm = DUST_HAZE_HEIGHT_M * DUST_M_TO_KM * Math.min(1, 2 * haze);
    item.haze.visible = haze > 0.005;
    if (item.haze.visible) {
      const hu = item.haze.material.uniforms;
      hu.uTime.value = t;
      hu.uScale.value.set(hazeRadiusKm, hazeRadiusKm, DUST_HAZE_HEIGHT_M * DUST_M_TO_KM);
      hu.uTau.value = DUST_HAZE_TAU * hazeOpacity * haze;
      dustApplySun(item.haze.material, hasSun, displayScale);
    }

    // Scour crater: grows with the mass already moved and stays once it is
    // there. The rim darkens the ground under it; the floor adds a little
    // light, so it has to carry the display gain like the dust does.
    const scourM = DUST_SCOUR_MAX_RADIUS_M * Math.min(1, Math.sqrt(Math.max(0, eroded) / erodedMax));
    for (const part of [item.scourRim, item.scourFloor]) {
      part.visible = scourM > 0.05;
      if (!part.visible) continue;
      part.position.set(0, 0, (part === item.scourRim ? 0.5 : 0.6) * DUST_DECAL_LIFT_M * DUST_M_TO_KM);
      part.scale.setScalar(scourM * DUST_M_TO_KM);
    }
    item.scourRim.material.uniforms.uOpacity.value = item.scourRim.userData.maxOpacity;
    item.scourFloor.material.uniforms.uOpacity.value = item.scourFloor.userData.maxOpacity;

    dustUpdateShadow(item, hasSun, sheetRadiusKm, hazeTopKm, Math.max(sheetStrength, 0.4 * haze));
  }

  // Haze level at `t` for spacecraft `s`, interpolated on the timeline.
  function dustHazeAt(t, s) {
    const { i, f } = frames.locate(t);
    const S = frames.sats;
    if (frames.count < 2) return hazeTrack[s];
    const a = hazeTrack[i * S + s], b = hazeTrack[(i + 1) * S + s];
    return a + f * (b - a);
  }

  return {
    group,
    update(t) {
      if (!visible) return;
      const lighting = dustLighting(options);
      const hasSun = dustSunDirection(lighting);
      const displayScale = dustDisplayScale(lighting);
      for (const item of emitters) dustUpdateEmitter(item, t, hasSun, displayScale);
    },
    setVisible(v) {
      visible = !!v;
      group.visible = visible;
    },
    get visible() { return visible; },
  };
}
