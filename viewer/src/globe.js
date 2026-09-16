// Central body: an oblate textured sphere driven by the sampled rotation table.
//
// Shading. The sphere carries a custom material rather than a Lambert one so
// that `setSun` can shade it with the same lunar reflectance the site terrain
// uses (viewer/terrain.js): a Lommel-Seeliger disk function with an opposition
// surge, which is what gives the Moon its flat, barely limb-darkened face and
// its brightening toward full. Until `setSun` is called -- a run without
// `frames.sun_dir` never calls it -- the material reproduces the Lambert
// shading the page had, so nothing else changes.
//
// Frames. three.js SphereGeometry puts its pole on +y and its texture seam on
// -x; SpaceAGORA's body-fixed frame has the pole on +z and longitude 0 on +x.
// `GEOMETRY_TO_BODY` is the fixed rotation between the two, applied on an
// inner mesh; the outer group carries q_pi(t), the body-to-inertial rotation
// from the sidecar (scalar-last, the same convention as the state
// quaternions: `rot(q)` in Julia is the passive inertial-to-body matrix, so
// the active body-to-inertial rotation is q itself, which is what
// Object3D.quaternion applies).
import * as THREE from 'three';
import { RotationTable } from 'viewer/data.js';
import { TERRAIN_REFLECTANCE_GLSL } from 'viewer/terrain.js';

// geometry (x, y, z) -> body (x, -z, y): pole +y -> +z, and u = 0.75 (east) -> +y.
export const GEOMETRY_TO_BODY = new THREE.Quaternion().setFromRotationMatrix(
  new THREE.Matrix4().set(
    1, 0, 0, 0,
    0, 0, -1, 0,
    0, 1, 0, 0,
    0, 0, 0, 1,
  ),
);

// Decode the image ourselves so an 8k or 16k texture can be downscaled to
// whatever this GPU allows (many mobile GPUs stop at 4096) instead of failing
// to upload, and so the page can report the size actually used.
function loadTextureWithinLimit(url, maxSize) {
  const texture = new THREE.Texture();
  const img = new Image();
  img.onload = () => {
    let source = img;
    if (img.width > maxSize) {
      const canvas = document.createElement('canvas');
      canvas.width = maxSize;
      canvas.height = Math.round(img.height * maxSize / img.width);
      canvas.getContext('2d').drawImage(img, 0, 0, canvas.width, canvas.height);
      source = canvas;
      console.info(`globe texture ${img.width}px downscaled to ${maxSize}px (GPU limit)`);
    }
    texture.image = source;
    texture.userData.width = source.width;
    texture.needsUpdate = true;
  };
  img.onerror = (e) => console.warn('globe texture failed to load', e);
  img.src = url;
  return texture;
}

const FALLBACK_COLORS = {
  earth: 0x3b6fb6, mars: 0xb5633a, venus: 0xd9b26a, titan: 0xc9a24a, moon: 0x9a9a9a,
};

const GLOBE_VERTEX = `
#include <common>
varying vec2 vGlobeUv;
varying vec3 vBodyDir;
varying vec3 vGlobeNormal;
varying vec3 vGlobeView;
#include <shadowmap_pars_vertex>
#include <logdepthbuf_pars_vertex>
void main() {
  vGlobeUv = uv;
  vBodyDir = normalize(vec3(position.x, -position.z, position.y));
  vec3 objectNormal = normal;
  vec3 transformedNormal = normalMatrix * objectNormal;
  vGlobeNormal = inverseTransformDirection(transformedNormal, viewMatrix);
  vec4 worldPosition = modelMatrix * vec4(position, 1.0);
  vGlobeView = cameraPosition - worldPosition.xyz;
  vec4 mvPosition = modelViewMatrix * vec4(position, 1.0);
  gl_Position = projectionMatrix * mvPosition;
  #include <shadowmap_vertex>
  #include <logdepthbuf_vertex>
}`;

// Built on the first call, not at load: the reflectance it quotes lives in
// terrain.js, and the CDN page concatenates the modules with globe.js first.
function globeFragment() {
  return `
#include <common>
#include <packing>
#include <lights_pars_begin>
#include <shadowmap_pars_fragment>
#include <shadowmask_pars_fragment>
#include <logdepthbuf_pars_fragment>
uniform sampler2D uMap;
uniform float uHasMap;
uniform float uMapOffset;
uniform vec3 uColor;
uniform vec4 uHole;
uniform float uHasHole;
uniform vec3 uSunWorld;
uniform float uLunar;
varying vec2 vGlobeUv;
varying vec3 vBodyDir;
varying vec3 vGlobeNormal;
varying vec3 vGlobeView;
${TERRAIN_REFLECTANCE_GLSL}
void main() {
  #include <logdepthbuf_fragment>
  if (uHasHole > 0.5) {
    float latDeg = degrees(asin(clamp(vBodyDir.z, -1.0, 1.0)));
    float lonDeg = degrees(atan(vBodyDir.y, vBodyDir.x));
    float lonMin = uHole.z, lonMax = uHole.w;
    if (lonDeg < lonMin - 180.0) lonDeg += 360.0;
    if (lonDeg > lonMax + 180.0) lonDeg -= 360.0;
    if (latDeg >= uHole.x && latDeg <= uHole.y && lonDeg >= lonMin && lonDeg <= lonMax) discard;
  }
  vec3 N = normalize(vGlobeNormal);
  vec3 V = normalize(vGlobeView);
  float mu = max(dot(N, V), 1.0e-3);
  vec3 albedo = uHasMap > 0.5 ? texture2D(uMap, vec2(vGlobeUv.x + uMapOffset, vGlobeUv.y)).rgb : uColor;
  float vehicleShadow = getShadowMask();
  vec3 direct = vec3(0.0);
  #if NUM_DIR_LIGHTS > 0
  for (int i = 0; i < NUM_DIR_LIGHTS; i++) {
    vec3 L = inverseTransformDirection(directionalLights[i].direction, viewMatrix);
    float mu0 = dot(N, L);
    if (mu0 <= 0.0) continue;
    float shade = dot(L, uSunWorld) > 0.99 ? vehicleShadow : 1.0;
    float response = uLunar > 0.5 ? saLunarReflectance(mu0, mu, dot(L, V)) : mu0;
    direct += directionalLights[i].color * response * shade;
  }
  #endif
  vec3 indirect = getAmbientLightIrradiance(ambientLightColor);
  #if NUM_HEMI_LIGHTS > 0
  vec3 nView = normalize(mat3(viewMatrix) * N);
  for (int i = 0; i < NUM_HEMI_LIGHTS; i++) indirect += getHemisphereLightIrradiance(hemisphereLights[i], nView);
  #endif
  gl_FragColor = vec4(albedo * RECIPROCAL_PI * (direct + indirect), 1.0);
  #include <tonemapping_fragment>
  #include <colorspace_fragment>
}`;
}

// Bodies whose surface scatters like regolith: an airless, dark, porous
// surface. Everything else keeps the Lambert shading, which is no worse for a
// cloud deck or an ocean than it ever was.
const GLOBE_REGOLITH_BODIES = ['moon', 'luna', 'mercury'];

// `textureEntry` is `{ url, lon_left_deg }` from the bundler (or null). The
// sphere's default mapping puts u = 0.5 on longitude 0, i.e. an image whose
// left edge is 180 W; an image with its left edge at L degrees east is
// shifted by -(L + 180) / 360 with wrap-around.
export function createGlobe(planet, textureEntry, options = {}) {
  const Re = planet.equatorial_radius_m / 1000, Rp = planet.polar_radius_m / 1000;
  const group = new THREE.Group();
  group.name = 'globe';

  const geometry = new THREE.SphereGeometry(1, 128, 64);
  const key = (planet.texture || planet.name || '').toLowerCase();
  const textureDataUrl = textureEntry ? (typeof textureEntry === 'string' ? textureEntry : textureEntry.url) : null;
  const uniforms = THREE.UniformsUtils.merge([THREE.UniformsLib.lights]);
  uniforms.uMap = { value: null };
  uniforms.uHasMap = { value: 0 };
  uniforms.uMapOffset = { value: 0 };
  uniforms.uColor = { value: new THREE.Color(FALLBACK_COLORS[key] ?? 0x888888) };
  uniforms.uSunWorld = { value: new THREE.Vector3(1, 0, 0) };
  uniforms.uLunar = { value: 0 };
  // A terrain patch replaces the sphere inside `options.hole` (a lat/lon box,
  // degrees): the fragment shader discards the sphere there. Geometry (x, y, z)
  // maps to body (x, -z, y).
  uniforms.uHole = { value: new THREE.Vector4(0, 0, 0, 0) };
  uniforms.uHasHole = { value: options.hole ? 1 : 0 };
  if (options.hole) {
    const h = options.hole;
    uniforms.uHole.value.set(h.lat_min, h.lat_max, h.lon_min, h.lon_max);
  }
  const material = new THREE.ShaderMaterial({
    uniforms, vertexShader: GLOBE_VERTEX, fragmentShader: globeFragment(), lights: true,
  });
  if (textureDataUrl) {
    const texture = loadTextureWithinLimit(textureDataUrl, options.maxTextureSize || 4096);
    texture.colorSpace = THREE.SRGBColorSpace;
    texture.anisotropy = options.anisotropy || 4;
    texture.generateMipmaps = true;
    texture.minFilter = THREE.LinearMipmapLinearFilter;
    const lonLeft = (textureEntry && typeof textureEntry === 'object' && Number.isFinite(textureEntry.lon_left_deg)) ? textureEntry.lon_left_deg : -180;
    texture.wrapS = THREE.RepeatWrapping;
    // The sphere's own mapping puts u = 0.5 on longitude 0; the shift is applied
    // in the shader, since a custom material has no texture matrix behind it.
    texture.offset.x = -(lonLeft + 180) / 360;
    uniforms.uMap.value = texture;
    uniforms.uHasMap.value = 1;
    uniforms.uMapOffset.value = texture.offset.x;
    material.map = texture;   // the exposure metering in lighting.js reads the globe texture off the material
  }
  const mesh = new THREE.Mesh(geometry, material);
  // The path tracer uploads materials by their `color`, which a ShaderMaterial
  // has none of: hand it the plain material this one replaces.
  mesh.userData.baseMaterial = new THREE.MeshLambertMaterial({ map: uniforms.uMap.value, color: uniforms.uHasMap.value ? 0xffffff : uniforms.uColor.value.getHex() });
  // The terrain hole is a fragment `discard`, which the shadow pass does not
  // run: a casting globe would drop its own sphere over the terrain patches.
  mesh.receiveShadow = true;
  mesh.castShadow = false;
  mesh.scale.set(Re, Rp, Re); // geometry y is the pole
  mesh.quaternion.copy(GEOMETRY_TO_BODY);
  group.add(mesh);

  // Graticule so a bare-color fallback still reads as a rotating body.
  const grid = new THREE.Group();
  const gridMaterial = new THREE.LineBasicMaterial({ color: 0xffffff, transparent: true, opacity: textureDataUrl ? 0.08 : 0.25 });
  for (let lat = -60; lat <= 60; lat += 30) {
    const r = Math.cos(THREE.MathUtils.degToRad(lat)), z = Math.sin(THREE.MathUtils.degToRad(lat));
    const pts = [];
    for (let k = 0; k <= 128; k++) {
      const a = (2 * Math.PI * k) / 128;
      pts.push(new THREE.Vector3(Re * r * Math.cos(a), Re * r * Math.sin(a), Rp * z));
    }
    grid.add(new THREE.Line(new THREE.BufferGeometry().setFromPoints(pts), gridMaterial));
  }
  for (let lon = 0; lon < 360; lon += 30) {
    const pts = [];
    for (let k = 0; k <= 64; k++) {
      const a = -Math.PI / 2 + (Math.PI * k) / 64;
      const r = Math.cos(a), z = Math.sin(a), l = THREE.MathUtils.degToRad(lon);
      pts.push(new THREE.Vector3(Re * r * Math.cos(l), Re * r * Math.sin(l), Rp * z));
    }
    grid.add(new THREE.Line(new THREE.BufferGeometry().setFromPoints(pts), gridMaterial));
  }
  grid.visible = options.graticule ?? true;
  group.add(grid);

  // Prime meridian and pole markers make the rotation direction obvious.
  const axis = new THREE.Line(
    new THREE.BufferGeometry().setFromPoints([new THREE.Vector3(0, 0, -1.15 * Rp), new THREE.Vector3(0, 0, 1.15 * Rp)]),
    new THREE.LineBasicMaterial({ color: 0x88ccff, transparent: true, opacity: 0.6 }),
  );
  axis.visible = options.axis ?? true;
  group.add(axis);

  const table = new RotationTable(planet.rotation, planet.spin_rad_s);
  const q = new Float32Array(4);
  const lunar = options.reflectance ? (options.reflectance === 'lunar' ? 1 : 0) : (GLOBE_REGOLITH_BODIES.includes(key) ? 1 : 0);

  return {
    group,
    mesh,
    grid,
    axis,
    radiusKm: Re,
    // Scene-space unit vector toward the Sun (lighting.js `direction`), once a
    // frame. Never called on a run without `frames.sun_dir`, which leaves the
    // material on its Lambert branch.
    setSun(dirScene) {
      if (!dirScene) return;
      uniforms.uSunWorld.value.copy(dirScene).normalize();
      uniforms.uLunar.value = lunar;
    },
    // Body-to-inertial quaternion at elapsed time t (seconds since epoch),
    // into a THREE.Quaternion or a 4-element array. (A typed array's `set`
    // takes an array, not four numbers, so the two cases differ.)
    rotationAt(t, out) {
      table.at(t, q);
      if (out.isQuaternion) return out.set(q[0], q[1], q[2], q[3]);
      out[0] = q[0]; out[1] = q[1]; out[2] = q[2]; out[3] = q[3];
      return out;
    },
    update(t) {
      table.at(t, q);
      group.quaternion.set(q[0], q[1], q[2], q[3]);
    },
  };
}
