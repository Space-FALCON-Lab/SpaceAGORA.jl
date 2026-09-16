// Thruster plumes: one exhaust plume per thruster glyph, driven by the firing
// levels the run recorded (`frames.thruster_level`, 0 to 1 per thruster in the
// sidecar's `spacecraft[i].thrusters` order).
//
// A plume hangs on the same link group as its cone glyph, so it follows the
// link's pose, and points along the thruster's exhaust direction. Its rated
// thrust sets how far it reaches and how it looks: a descent engine burns long
// and orange-white, an attitude jet puffs short and blue-white. Both are two
// additive cones - a bright core inside a translucent shroud - with a shader
// that fades them toward the tip and flickers them in time, so nothing writes
// depth and the plumes never occlude the vehicle.
//
// How much of that look survives depends on the air around the nozzle. A plume
// is luminous because the exhaust has ambient gas to shock, mix with and burn
// against; in vacuum a hypergolic engine shows almost nothing -- the Apollo
// films of the descent engine record a bare bell against the ground, not a
// flame -- so the page fades the shroud out and leaves a faint blue-white core
// as the ambient density falls (`frames.density`, or the scene's lack of an
// atmosphere model when the run saved none).
//
// Every top-level name here is prefixed `plume`/`PLUME_`: the CDN page builder
// concatenates all viewer modules into one script scope.
import * as THREE from 'three';
import { decodeFloat32 } from 'viewer/data.js';

const PLUME_MIN_LEVEL = 0.0005;        // below this a thruster counts as idle
// A jet's level is the fraction of its rating the controller is asking for. An
// attitude jet holding a vehicle steady asks for well under a percent of it,
// which a linear mapping would draw as nothing at all, so the plume's reach and
// brightness follow the level through this exponent: monotone, full at level 1,
// and still visible at the fractions a real duty-cycled jet fires at.
const PLUME_LEVEL_GAMMA = 0.4;
// Reach: 1.9 m for a 1 kN thruster, times (rating / 1 kN)^0.4. The exponent is
// flatter than the thrust ratio so a 45 kN descent engine reaches about six times
// as far as a 445 N attitude jet rather than ten, which keeps the jets' puffs
// clear of the vehicle while the engine's plume still dominates.
const PLUME_REFERENCE_THRUST_N = 1000.0;
const PLUME_REFERENCE_LENGTH_M = 1.9;
const PLUME_LENGTH_EXPONENT = 0.4;
const PLUME_MIN_LENGTH_M = 0.25;
const PLUME_WIDTH_FRACTION = 0.16;     // tip radius as a fraction of the reach
const PLUME_CORE_LENGTH = 0.55;        // the bright core, as a fraction of the shroud
const PLUME_CORE_WIDTH = 0.45;
// Color by rated thrust: blue-white at and below a few hundred newtons (an
// attitude jet), orange-white from about 10 kN up (a main engine).
const PLUME_COLOR_PIVOT_N = 300.0;
const PLUME_COLOR_DECADES = 1.5;
const PLUME_COOL_COLOR = new THREE.Color(0.62, 0.78, 1.0);
const PLUME_HOT_COLOR = new THREE.Color(1.0, 0.72, 0.38);
const PLUME_CORE_COLOR = new THREE.Color(1.0, 0.96, 0.9);
// Ambient density and what a plume looks like in it. A descent engine at sea
// level draws a bright sooty column; by about 70 km (1e-5 kg/m^3) it is a thin
// diffuse glow; below roughly 1e-9 kg/m^3 there is nothing left to excite and
// the exhaust is invisible. The look follows the logarithm of the density
// between those two, which is where the transition actually happens.
const PLUME_DENSITY_FULL_KG_M3 = 1e-5;
const PLUME_DENSITY_VACUUM_KG_M3 = 1e-9;
// What is left in vacuum: a few percent of the core's brightness, almost none
// of the shroud's, and the color of the recombining hypergolic gas itself.
const PLUME_VACUUM_CORE_SCALE = 0.08;
const PLUME_VACUUM_SHROUD_SCALE = 0.015;
const PLUME_VACUUM_COLOR = new THREE.Color(0.74, 0.86, 1.0);
// A plume is emissive, and this shader writes display values rather than the
// radiance the lit materials write: how bright an exhaust looks is a property
// of the exhaust, not of how the camera around it is exposed. So the exposure
// the lighting module sets -- about ten on a page lit by a grazing lunar sun,
// where it is there to lift a dark ground to a mid tone -- is deliberately not
// applied on top; it would lift the faint vacuum core back into a flame. It is
// still read, because an exposure below one means a scene being compressed
// rather than lifted, and the plume is dimmed with it. The sun's irradiance
// dims it the same way, so a plume at Saturn keeps its contrast against the
// darker ground there.
const PLUME_SCALE_MIN = 0.05;
const PLUME_REFERENCE_IRRADIANCE_W_M2 = 1361.0;
const PLUME_IDLE_GLYPH_OPACITY = 0.35;
const PLUME_IDLE_GLYPH_EMISSIVE = 0.0;
const PLUME_LIVE_GLYPH_EMISSIVE = 0.9;

// Unit cone: apex at the origin, base of radius 1 at y = 1, so the mesh scale
// alone sets the reach and the width and `vT` runs 0 (throat) to 1 (tip).
function plumeUnitCone() {
  const g = new THREE.ConeGeometry(1, 1, 16, 4, true);
  g.rotateX(Math.PI);
  g.translate(0, 0.5, 0);
  return g;
}

const PLUME_VERTEX = `
  #include <common>
  varying float vT;
  varying float vAng;
  #include <logdepthbuf_pars_vertex>
  void main() {
    vT = clamp(position.y, 0.0, 1.0);
    vAng = atan(position.z, position.x);
    vec4 mvPosition = modelViewMatrix * vec4(position, 1.0);
    gl_Position = projectionMatrix * mvPosition;
    #include <logdepthbuf_vertex>
  }
`;

// Brightness falls off along the plume and flickers in time; two incommensurate
// sine terms plus a cheap hash keep the noise from reading as a pulse.
const PLUME_FRAGMENT = `
  #include <common>
  uniform vec3 uColor;
  uniform float uLevel;
  uniform float uTime;
  uniform float uSeed;
  uniform float uIntensity;
  uniform float uVisibility;
  varying float vT;
  varying float vAng;
  #include <logdepthbuf_pars_fragment>
  float plumeHash(vec2 p) {
    return fract(sin(dot(p, vec2(12.9898, 78.233))) * 43758.5453);
  }
  void main() {
    #include <logdepthbuf_fragment>
    float t = uTime * 9.0 + uSeed;
    float flicker = 0.80
      + 0.12 * sin(t * 3.7 + vT * 7.0)
      + 0.08 * sin(t * 11.3 + vAng * 2.0)
      + 0.10 * plumeHash(vec2(floor(t * 4.0), floor(vT * 12.0) + uSeed));
    // A gentle taper: the plume has to stay bright well downstream, because a
    // thruster that sits inside a 3D model only shows from the model's skin out.
    float taper = pow(1.0 - vT, 0.9);
    float throat = smoothstep(0.0, 0.08, vT);
    float a = uIntensity * uVisibility * uLevel * taper * throat * flicker;
    if (a <= 0.002) discard;
    gl_FragColor = vec4(uColor * (0.7 + 0.6 * uLevel), a);
  }
`;

function plumeMaterial(color, intensity, seed) {
  return new THREE.ShaderMaterial({
    uniforms: {
      uColor: { value: color.clone() },
      uLevel: { value: 0 },
      uTime: { value: 0 },
      uSeed: { value: seed },
      uIntensity: { value: intensity },
      uVisibility: { value: 1 },
    },
    vertexShader: PLUME_VERTEX,
    fragmentShader: PLUME_FRAGMENT,
    transparent: true,
    depthWrite: false,
    blending: THREE.AdditiveBlending,
    side: THREE.DoubleSide,
  });
}

// Blue-white for a small jet, orange-white for a main engine, on a log scale
// in the rated thrust so the two ends of a vehicle's set read apart.
function plumeColorFor(maxThrustN) {
  const f = Math.max(0, Math.min(1, Math.log10(Math.max(maxThrustN, 1) / PLUME_COLOR_PIVOT_N) / PLUME_COLOR_DECADES));
  return PLUME_COOL_COLOR.clone().lerp(PLUME_HOT_COLOR, f);
}

// Level as the plume draws it (see PLUME_LEVEL_GAMMA); the panel and the plots
// keep reporting the raw level.
function plumeVisualLevel(level) {
  return Math.pow(Math.max(0, Math.min(1, level)), PLUME_LEVEL_GAMMA);
}

function plumeLengthFor(maxThrustN) {
  const ratio = Math.max(maxThrustN, 0) / PLUME_REFERENCE_THRUST_N;
  return Math.max(PLUME_MIN_LENGTH_M, PLUME_REFERENCE_LENGTH_M * Math.pow(ratio, PLUME_LENGTH_EXPONENT));
}

/**
 * How atmospheric a plume should look at an ambient density of `density`
 * kg/m^3: 1 at and above PLUME_DENSITY_FULL_KG_M3, 0 at and below
 * PLUME_DENSITY_VACUUM_KG_M3, logarithmic between them. A non-finite density
 * means the run saved none and the caller decides from the scene instead.
 */
export function plumeAtmosphereFactor(density) {
  if (!Number.isFinite(density) || density <= PLUME_DENSITY_VACUUM_KG_M3) return 0;
  if (density >= PLUME_DENSITY_FULL_KG_M3) return 1;
  const span = Math.log10(PLUME_DENSITY_FULL_KG_M3 / PLUME_DENSITY_VACUUM_KG_M3);
  return Math.min(1, Math.max(0, Math.log10(density / PLUME_DENSITY_VACUUM_KG_M3) / span));
}

// The lighting handle, which main.js passes as a function because the module
// is built after this one (it has to light the bodies that already exist).
// Missing or half-built handles are fine: the plumes then keep the brightness
// they had before there was a lighting module at all.
function plumeLighting(options) {
  const l = options.lighting;
  return typeof l === 'function' ? l() : (l || null);
}

// Brightness the page draws the plume at, relative to the look it was tuned
// for: never above it, and dimmed on a page with a weaker sun or a compressed
// exposure (see the note on PLUME_SCALE_MIN). 1 without a lighting handle.
function plumeExposureScale(lighting) {
  if (!lighting) return 1;
  const exposure = Number.isFinite(lighting.exposure) && lighting.exposure > 0 ? lighting.exposure : 1;
  const irradiance = Number.isFinite(lighting.sunIrradiance) && lighting.sunIrradiance > 0
    ? lighting.sunIrradiance : PLUME_REFERENCE_IRRADIANCE_W_M2;
  const scale = Math.min(1, irradiance / PLUME_REFERENCE_IRRADIANCE_W_M2) * Math.min(1, exposure);
  return Math.min(1, Math.max(PLUME_SCALE_MIN, scale));
}

// One plume: a translucent shroud with a brighter core, apex on the thruster
// and axis along its exhaust direction, in the link group's own (meter) frame.
function plumeBuild(geometry, spec, index) {
  const group = new THREE.Group();
  group.name = `plume:${index}`;
  const length = plumeLengthFor(spec.max_thrust_n);
  const radius = PLUME_WIDTH_FRACTION * length;
  const color = plumeColorFor(spec.max_thrust_n);
  const seed = 7.3 * index + 1.1;
  const shroudColor = color.clone();
  const coreColor = PLUME_CORE_COLOR.clone().lerp(color, 0.35);
  const shroud = new THREE.Mesh(geometry, plumeMaterial(shroudColor, 0.75, seed));
  shroud.scale.set(radius, length, radius);
  const core = new THREE.Mesh(geometry, plumeMaterial(coreColor, 1.0, seed + 3.7));
  core.scale.set(PLUME_CORE_WIDTH * radius, PLUME_CORE_LENGTH * length, PLUME_CORE_WIDTH * radius);
  for (const m of [shroud, core]) {
    m.frustumCulled = false;
    m.renderOrder = 6;
    group.add(m);
  }
  group.position.set(spec.location_m[0], spec.location_m[1], spec.location_m[2]);
  const dir = new THREE.Vector3(spec.direction[0], spec.direction[1], spec.direction[2]);
  if (dir.lengthSq() > 0) group.quaternion.setFromUnitVectors(new THREE.Vector3(0, 1, 0), dir.normalize());
  group.visible = false;
  return { group, shroud, core, length, radius, shroudColor, coreColor };
}

/**
 * createPlumes(sidecar, frames, lod, options)
 *   options.raw      the undecoded frames payload (carries thruster_level / thruster_counts)
 *   options.enabled  start visible (default true)
 *   options.lighting the lighting handle (or a function returning it), read
 *                    each update for `sunIrradiance` and `exposure`; absent
 *                    keeps the brightness the page had before it existed
 *   options.vacuum   force the vacuum look (true) or the atmospheric one
 *                    (false); by default the run's density decides
 * Returns { update(t), setVisible(v), visible, levelsAt(t, s, out), counts(s), available }.
 */
export function createPlumes(sidecar, frames, lod, options = {}) {
  const raw = options.raw || {};
  const counts = Array.isArray(raw.thruster_counts) ? raw.thruster_counts.map((n) => Number(n) || 0) : null;
  const data = counts && raw.thruster_level ? decodeFloat32(raw.thruster_level) : null;
  const total = counts ? counts.reduce((a, b) => a + b, 0) : 0;
  const offsets = [];
  {
    let acc = 0;
    for (const n of counts || []) { offsets.push(acc); acc += n; }
  }
  const available = !!data && total > 0 && data.length >= frames.count * total;

  // Levels of spacecraft `s` at time t, linear between the bracketing frames.
  function levelsAt(t, s, out) {
    if (!available || !counts[s]) return null;
    const n = counts[s], o = offsets[s];
    const { i, f } = frames.locate(t);
    if (frames.count < 2) {
      for (let k = 0; k < n; k++) out[k] = data[o + k];
      return out;
    }
    const a = i * total + o, b = (i + 1) * total + o;
    for (let k = 0; k < n; k++) out[k] = data[a + k] + f * (data[b + k] - data[a + k]);
    return out;
  }

  // How much air the plume has around it. A run that saved density answers per
  // frame and per spacecraft; one that did not falls back to the scene, where
  // a missing atmosphere block means the run flew with `NoAtmosphereModel` --
  // vacuum -- and any other model keeps the full atmospheric look.
  const densityAvailable = typeof frames.hasScalar === 'function' && frames.hasScalar('density');
  const sceneAtmosphere = options.vacuum === true ? 0
    : options.vacuum === false ? 1
    : (sidecar && sidecar.atmosphere ? 1 : 0);
  function plumeAtmosphereAt(t, s) {
    if (options.vacuum !== undefined || !densityAvailable) return sceneAtmosphere;
    return plumeAtmosphereFactor(frames.scalarAt('density', t, s));
  }

  const geometry = plumeUnitCone();
  const items = [];
  let count = 0;
  for (let s = 0; s < frames.sats; s++) {
    const item = lod.items[s];
    const specs = item && item.spec ? item.spec.thrusters || [] : [];
    const plumes = [];
    // lod.js skips a glyph whose link group is missing, so walk its cone list
    // with its own cursor to keep plume and glyph paired.
    const glyphs = item ? item.glyphs.thrusters : [];
    let gi = 0;
    if (available && counts[s] === specs.length) {
      specs.forEach((spec) => {
        const group = item.linkGroups[spec.link - 1];
        if (!group) { plumes.push(null); return; }
        const cone = glyphs[gi++] || null;
        const built = plumeBuild(geometry, spec, count++);
        group.add(built.group);
        // The static cone glyph stays, dimmed while its thruster is idle.
        if (cone && cone.material) {
          cone.material.transparent = true;
          cone.material.needsUpdate = true;
          cone.userData.plumeBaseEmissive = cone.material.emissiveIntensity ?? 0.35;
        }
        plumes.push({ ...built, cone });
      });
    }
    items.push({ plumes, buffer: new Float32Array(Math.max(1, counts ? counts[s] : 1)) });
  }

  let visible = (options.enabled ?? true) && available;
  let clock = 0;
  let last = null;

  function applyIdleGlyphs(item) {
    for (const p of item.plumes) {
      if (!p || !p.cone || !p.cone.material) continue;
      p.cone.material.opacity = 1;
      p.cone.material.emissiveIntensity = p.cone.userData.plumeBaseEmissive ?? 0.35;
    }
  }

  return {
    available,
    group: null,   // the plumes live on the assemblies' link groups, not a group of their own
    counts(s) { return counts && counts[s] ? counts[s] : 0; },
    levelsAt,
    get visible() { return visible; },
    setVisible(v) {
      visible = !!v && available;
      if (!visible) {
        for (const item of items) {
          for (const p of item.plumes) if (p) p.group.visible = false;
          applyIdleGlyphs(item);
        }
      }
    },
    // One frame: advance the flicker clock, then set every plume's level.
    update(t) {
      const now = (typeof performance !== 'undefined' ? performance.now() : Date.now()) / 1000;
      if (last !== null) clock += Math.min(0.25, Math.max(0, now - last));
      last = now;
      if (!visible) return;
      const lighting = plumeLighting(options);
      const exposureScale = plumeExposureScale(lighting);
      for (let s = 0; s < items.length; s++) {
        const item = items[s];
        if (!item.plumes.length) continue;
        // Ambient density: full plume in air, a faint core in vacuum.
        const atmos = plumeAtmosphereAt(t, s);
        const shroudScale = exposureScale * (PLUME_VACUUM_SHROUD_SCALE + (1 - PLUME_VACUUM_SHROUD_SCALE) * atmos);
        const coreScale = exposureScale * (PLUME_VACUUM_CORE_SCALE + (1 - PLUME_VACUUM_CORE_SCALE) * atmos);
        const levels = levelsAt(t, s, item.buffer);
        for (let k = 0; k < item.plumes.length; k++) {
          const p = item.plumes[k];
          if (!p) continue;
          const level = levels && Number.isFinite(levels[k]) ? Math.max(0, Math.min(1, levels[k])) : 0;
          const firing = level > PLUME_MIN_LEVEL;
          p.group.visible = firing;
          const cone = p.cone;
          if (cone && cone.material) {
            cone.material.opacity = firing ? 1 : PLUME_IDLE_GLYPH_OPACITY;
            cone.material.emissiveIntensity = firing
              ? PLUME_LIVE_GLYPH_EMISSIVE * plumeVisualLevel(level)
              : PLUME_IDLE_GLYPH_EMISSIVE;
          }
          if (!firing) continue;
          // Reach and width grow with the level; the core trails the shroud.
          const shown = plumeVisualLevel(level);
          const grow = 0.3 + 0.7 * shown;
          const wide = p.radius * (0.65 + 0.35 * shown);
          p.shroud.scale.set(wide, p.length * grow, wide);
          p.core.scale.set(
            PLUME_CORE_WIDTH * wide,
            PLUME_CORE_LENGTH * p.length * grow,
            PLUME_CORE_WIDTH * wide,
          );
          // Size carries the magnitude; brightness keeps a floor so a jet on a
          // fraction of a percent of its rating still reads as firing.
          const glow = 0.3 + 0.7 * shown;
          p.shroud.material.uniforms.uLevel.value = glow;
          p.core.material.uniforms.uLevel.value = glow;
          // In vacuum only the recombining core is left, and it is blue-white.
          p.shroud.material.uniforms.uVisibility.value = shroudScale;
          p.core.material.uniforms.uVisibility.value = coreScale;
          p.shroud.material.uniforms.uColor.value.copy(p.shroudColor).lerp(PLUME_VACUUM_COLOR, 1 - atmos);
          p.core.material.uniforms.uColor.value.copy(p.coreColor).lerp(PLUME_VACUUM_COLOR, 1 - atmos);
          p.shroud.material.uniforms.uTime.value = clock;
          p.core.material.uniforms.uTime.value = clock;
        }
      }
    },
  };
}
