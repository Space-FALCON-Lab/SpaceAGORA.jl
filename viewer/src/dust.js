// Blowing regolith under a landing vehicle: the dust sheet the descent engine's
// plume drives off the surface, the ground it obscures, and the scour mark it
// leaves. Driven entirely by the `frames.plume` block the plume-surface
// effector writes (see src/dynamics/coupled/force_torque_models/
// plume_surface_interaction.jl), so the sheet only appears where the physics
// says soil is actually moving.
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
// The sheet is a GPU particle system: the geometry carries only per-particle
// seeds and every position is derived in the vertex shader from those seeds and
// the current time, so a frame update writes uniforms and nothing else.
import * as THREE from 'three';
import { rotateByConjugate } from 'viewer/data.js';

const DUST_M_TO_KM = 1e-3;
const DUST_PARTICLE_COUNT = 4096;
const DUST_MAX_EMITTERS = 8;
const DUST_LIFETIME_S = 2.4;
const DUST_ELEVATION_MIN_RAD = 1.0 * Math.PI / 180;   // the sheet hugs the ground: a
const DUST_ELEVATION_MAX_RAD = 3.0 * Math.PI / 180;   // few degrees above the local slope
const DUST_GRAVITY_M_S2 = 1.62;                        // lunar, unless the payload says otherwise
const DUST_SIZE_KM = 4.5e-4;                           // a puff about half a meter across
const DUST_PROJECTION_SCALE = 820;                     // pixels per unit angle, ~800 px tall view at 45 deg
const DUST_POINT_MIN_PX = 1.2;
const DUST_POINT_MAX_PX = 34.0;
const DUST_COLOR = 0xc9bda6;
const DUST_SCOUR_COLOR = 0x6b6153;
const DUST_MAX_OPACITY = 0.5;
const DUST_OBSCURATION_OPACITY = 0.42;
const DUST_SCOUR_OPACITY = 0.6;
const DUST_SCOUR_MAX_RADIUS_M = 14.0;
const DUST_SHEET_LIFT_M = 0.7;                         // the saltating layer has thickness, and this
                                                       // keeps the sheet out of the terrain patches' facets
const DUST_DECAL_LIFT_M = 1.6;                         // clears the terrain patches' triangle facets
const DUST_PLUME_TAN_HALF_ANGLE = 0.4663;              // tan(25 deg), the plume core the effector models
const DUST_OBSCURATION_FOOTPRINTS = 2.5;               // haze disk radius, in plume footprint radii
const DUST_UP_AXIS = new THREE.Vector3(0, 0, 1);

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
    vec4 mv = modelViewMatrix * vec4(p, 1.0);
    gl_Position = projectionMatrix * mv;
    #include <logdepthbuf_vertex>
    gl_PointSize = clamp(uSizeKm * uProjScale / max(1e-7, -mv.z), ${DUST_POINT_MIN_PX.toFixed(1)}, ${DUST_POINT_MAX_PX.toFixed(1)});
  }
`;

const DUST_FRAGMENT_SHADER = `
  #include <common>
  #include <logdepthbuf_pars_fragment>
  uniform vec3 uColor;
  uniform float uOpacity;
  varying float vAlpha;
  void main() {
    #include <logdepthbuf_fragment>
    vec2 d = gl_PointCoord - vec2(0.5);
    float r2 = dot(d, d);
    if (r2 > 0.25) discard;
    float soft = 1.0 - smoothstep(0.0, 0.25, r2);
    gl_FragColor = vec4(uColor, vAlpha * uOpacity * soft);
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
    },
    vertexShader: DUST_VERTEX_SHADER,
    fragmentShader: DUST_FRAGMENT_SHADER,
    transparent: true,
    depthWrite: false,
    blending: THREE.AdditiveBlending,
  });
}

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

// A flat disk lying on the ground with a soft rim: the obscuration haze (which
// adds light, like dust catching the Sun) and the scour mark (which darkens).
function dustDecal(color, opacity, renderOrder, core, additive) {
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
    blending: additive ? THREE.AdditiveBlending : THREE.NormalBlending,
  });
  const mesh = new THREE.Mesh(new THREE.CircleGeometry(1, 64), material);
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
 * `gravity_m_s2`, `projectionScale`, `maxEmitters`.
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

  const emitters = [];
  for (let s = 0; s < count; s++) {
    const pivot = new THREE.Group();
    const material = dustMaterial(gravityKmS2, projScale);
    const points = new THREE.Points(dustSeedGeometry(DUST_PARTICLE_COUNT), material);
    points.frustumCulled = false;
    points.renderOrder = 30;
    const obscuration = dustDecal(DUST_COLOR, DUST_OBSCURATION_OPACITY, 22, 0.0, true);
    const scour = dustDecal(DUST_SCOUR_COLOR, DUST_SCOUR_OPACITY, 21, 0.55, false);
    pivot.add(scour, obscuration, points);
    pivot.visible = false;
    group.add(pivot);
    emitters.push({ s, pivot, points, material, obscuration, scour });
  }
  if (parent && typeof parent.add === 'function') parent.add(group);

  const posKm = new Float64Array(3), qAtt = new Float32Array(4), qPlanet = new Float32Array(4);
  const bodyKm = new Float64Array(3);
  const axis = new THREE.Vector3(), impact = new THREE.Vector3(), up = new THREE.Vector3();
  const quat = new THREE.Quaternion();
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

  function dustUpdateEmitter(item, t) {
    const s = item.s;
    const erosion = frames.plumeAt('erosion_kg_s', t, s);
    if (!(erosion > 0)) { item.pivot.visible = false; return; }
    const height = frames.plumeAt('height_m', t, s);
    const ejecta = frames.plumeAt('ejecta_mps', t, s);
    const eroded = frames.plumeAt('eroded_kg', t, s);
    if (!Number.isFinite(height) || !dustImpactPoint(t, s, height)) { item.pivot.visible = false; return; }

    const strength = Math.min(1, erosion / erosionMax);
    const speedKmS = Math.max(1e-4, (Number.isFinite(ejecta) ? ejecta : 0) * DUST_M_TO_KM);
    const reachKm = speedKmS * DUST_LIFETIME_S;

    item.pivot.visible = true;
    item.pivot.position.copy(impact);
    item.pivot.quaternion.setFromUnitVectors(DUST_UP_AXIS, up);

    const u = item.material.uniforms;
    u.uTime.value = t;
    u.uSpeed.value = speedKmS;
    u.uFadeKm.value = 0.22 * reachKm;
    u.uOpacity.value = DUST_MAX_OPACITY * Math.sqrt(strength);
    item.points.geometry.setDrawRange(0, Math.max(16, Math.round(DUST_PARTICLE_COUNT * strength)));

    // Ground obscuration: a haze disk over the sheet, opacity following the rate.
    const lift = DUST_DECAL_LIFT_M * DUST_M_TO_KM;
    const footprintKm = Math.max(1e-5, height * DUST_M_TO_KM * DUST_PLUME_TAN_HALF_ANGLE);
    item.obscuration.visible = true;
    item.obscuration.position.set(0, 0, lift);
    item.obscuration.scale.setScalar(DUST_OBSCURATION_FOOTPRINTS * footprintKm);
    item.obscuration.material.uniforms.uOpacity.value = item.obscuration.userData.maxOpacity * strength;

    // Scour mark: grows with the mass already moved and stays once it is there.
    const scourM = DUST_SCOUR_MAX_RADIUS_M * Math.min(1, Math.sqrt(Math.max(0, eroded) / erodedMax));
    if (scourM > 0.05) {
      item.scour.visible = true;
      item.scour.position.set(0, 0, 0.6 * lift);
      item.scour.scale.setScalar(scourM * DUST_M_TO_KM);
      item.scour.material.uniforms.uOpacity.value = item.scour.userData.maxOpacity;
    } else {
      item.scour.visible = false;
    }
  }

  return {
    group,
    update(t) {
      if (!visible) return;
      for (const item of emitters) dustUpdateEmitter(item, t);
    },
    setVisible(v) {
      visible = !!v;
      group.visible = visible;
    },
    get visible() { return visible; },
  };
}
