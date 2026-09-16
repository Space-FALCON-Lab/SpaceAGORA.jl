// Site terrain: nested square patches of the surface around a landing site,
// displaced by the digital elevation grids the page carries and draped with
// imagery that sharpens toward the center (from the global mosaic down to
// sub-meter mosaics). The patches sit in the globe group, so they rotate
// with the body; the globe itself is cut open under the outermost patch by
// the hole the globe material discards. All levels are drawn all the time:
// the finest level covers a few hundred meters, so the resolution rises by
// itself as the camera closes in.
//
// The ground is lit by a lunar reflectance rather than a Lambert law, off a
// normal map the page derives from the height grids (finer than the 96 x 96
// mesh), and it shadows itself: the fragment shader marches toward the Sun
// through the same height grids, uploaded as textures, and multiplies the
// horizon test into the vehicle's shadow map. `setSun` drives all three; a run
// without `frames.sun_dir` never calls it and the ground keeps the Lambert
// look it had, lit by whatever lights the scene holds.
//
// Payload contract (viewer_bundle.jl `terrain_payload`):
//   terrain.site               { lat_deg, lon_deg, height_m, name }
//   terrain.reference_radius_m
//   terrain.grids[]            { name, rows, cols, lat_min, lat_max, lon_min, lon_max, heights (base64 Float32, row-major north to south) }, finest first
//   terrain.imagery[]          { lat_min, lat_max, lon_min, lon_max, width, height, m_per_px, url (data: JPEG), albedo_url (data: JPEG, optional) }, coarse to fine
// `albedo_url` is `url` with the mosaic's own illumination divided out
// (scripts/dev/terrain/fetch_moon_site.py `derive_albedo`); it is what the page
// drapes when it is there, since the page lights the surface itself.
import * as THREE from 'three';
import { decodeFloat32 } from 'viewer/data.js';

const PATCH_SEGMENTS = 96;
const TERRAIN_M_TO_KM = 1e-3;

// ---------------------------------------------------------------------------
// Lunar reflectance
// ---------------------------------------------------------------------------
// Regolith is not a Lambert surface. It is very nearly backscattering: the Moon
// shows almost no limb darkening at full phase and brightens sharply as the
// phase angle closes, which is why a Lambert ball never looks like the Moon.
// What is used here is the lunar-Lambert law -- Lommel-Seeliger scattering
// blended with Lambert by the phase-dependent coefficient L(g) fitted by
// McEwen (1991), "Photometric functions for photoclinometry and other
// applications", Icarus 92, 298, the same form ISIS applies when it
// photometrically normalizes LROC images (Sato et al. 2014, JGR Planets 119,
// 1775) -- times a single Henyey-Greenstein particle phase function and Hapke's
// shadow-hiding opposition surge B(g) = B0 / (1 + tan(g/2)/h) (Hapke 1986,
// Icarus 67, 264; Helfenstein & Veverka 1987 fit B0 near 1 and h near 0.06 for
// the Moon). The whole phase term is normalized to 1 at g = 45 degrees so a
// scene under an ordinary phase keeps the brightness the Lambert page had and
// the opposition surge reads as the brightening it is.
const TERRAIN_HG_ASYMMETRY = -0.29;     // Henyey-Greenstein xi, backscattering
const TERRAIN_OPPOSITION_B0 = 1.0;      // surge amplitude
const TERRAIN_OPPOSITION_H = 0.06;      // surge angular width
const TERRAIN_PHASE_REF_DEG = 45.0;     // phase the reflectance is normalized at

function terrainPhaseTerm(gRad) {
  const xi = TERRAIN_HG_ASYMMETRY;
  const hg = (1 - xi * xi) / Math.pow(Math.max(1 + 2 * xi * Math.cos(gRad) + xi * xi, 1e-4), 1.5);
  const surge = TERRAIN_OPPOSITION_B0 / (1 + Math.tan(0.5 * gRad) / TERRAIN_OPPOSITION_H);
  return hg * (1 + surge);
}
const TERRAIN_PHASE_NORM = terrainPhaseTerm(TERRAIN_PHASE_REF_DEG * Math.PI / 180);

// Shared with globe.js, which shades the rest of the body the same way.
// `mu0` = cos(incidence), `mu` = cos(emission), `cosg` = cos(phase).
export const TERRAIN_REFLECTANCE_GLSL = `
float saLunarReflectance(float mu0, float mu, float cosg) {
  if (mu0 <= 0.0) return 0.0;
  float g = acos(clamp(cosg, -1.0, 1.0));
  float gDeg = degrees(g);
  // McEwen (1991) limb-darkening coefficient
  float L = clamp(1.0 - 1.9e-2 * gDeg + 2.42e-4 * gDeg * gDeg - 1.46e-6 * gDeg * gDeg * gDeg, 0.0, 1.0);
  float disk = L * (2.0 * mu0 / max(mu0 + max(mu, 1e-3), 1e-3)) + (1.0 - L) * mu0;
  float xi = ${TERRAIN_HG_ASYMMETRY.toFixed(4)};
  float hg = (1.0 - xi * xi) / pow(max(1.0 + 2.0 * xi * cosg + xi * xi, 1e-4), 1.5);
  float surge = ${TERRAIN_OPPOSITION_B0.toFixed(4)} / (1.0 + tan(min(0.5 * g, 1.5533)) / ${TERRAIN_OPPOSITION_H.toFixed(4)});
  return disk * hg * (1.0 + surge) / ${TERRAIN_PHASE_NORM.toFixed(6)};
}
`;

// Ray-marched terrain shadows. The march is a horizon test in the site's local
// east/north/up frame: the height of the ray above the reference sphere is
// h0 + s cos z + s^2 sin^2 z / 2R (the sphere curving away under a straight
// ray), and the ground under it comes from the same grids the mesh was
// displaced with. Steps grow geometrically from a few meters to a few
// kilometers, so a boulder-scale ledge and a crater rim a kilometer up-sun both
// land inside 64 samples; the sample stops at the first solid hit. The soft
// edge is physical: at distance s the Sun's half-degree disk is blurred over
// s * 4.65e-3, which is what the smoothstep width is set to.
const TERRAIN_SHADOW_STEPS = 64;
const TERRAIN_SHADOW_FIRST_KM = 0.008;   // first sample, up-sun of the fragment
const TERRAIN_SHADOW_REACH_KM = 6.0;     // last sample
const TERRAIN_SHADOW_BIAS_M = 2.0;       // clears the grids' own sampling noise
const TERRAIN_SUN_ANGULAR_RADIUS = 0.00465;   // radians, the Sun from 1 au
// Normal map: one texture per imagery level, sampled at the resolution of the
// finest height grid that covers the level, so the relief the DEM resolves
// shades correctly however coarse the drawn mesh is.
const TERRAIN_NORMAL_MIN = 64;
const TERRAIN_NORMAL_MAX = 512;

// Height grids as textures. A grid is uploaded as 16-bit fixed point across the
// red and green bytes of an RGBA8 texture -- 0.05 m over the LOLA window, 2 mm
// over the NAC one -- read with nearest filtering, which needs no float-texture
// extension anywhere and keeps the march at one tap a step.
function terrainHeightTexture(grid) {
  const n = grid.rows * grid.cols;
  let lo = Infinity, hi = -Infinity;
  for (let k = 0; k < n; k++) { const v = grid.data[k]; if (v < lo) lo = v; if (v > hi) hi = v; }
  if (!(lo < hi)) { lo = Number.isFinite(lo) ? lo : 0; hi = lo + 1; }
  const range = hi - lo;
  const bytes = new Uint8Array(4 * n);
  for (let k = 0; k < n; k++) {
    const q = Math.max(0, Math.min(65535, Math.round((grid.data[k] - lo) / range * 65535)));
    bytes[4 * k] = q >> 8; bytes[4 * k + 1] = q & 255; bytes[4 * k + 3] = 255;
  }
  const texture = new THREE.DataTexture(bytes, grid.cols, grid.rows, THREE.RGBAFormat);
  texture.magFilter = THREE.NearestFilter;
  texture.minFilter = THREE.NearestFilter;
  texture.generateMipmaps = false;
  texture.needsUpdate = true;
  return { texture, base: lo, range };
}

// Per-level normal map. The height field is sampled once per texel (plus a ring
// outside the level, for the differences) and turned into a unit normal in the
// local east/north/up frame, the ordinary tangent-space encoding: 8 bits a
// component is about half a degree of slope, which is finer than the grids
// themselves resolve.
function terrainNormalTexture(level, heightAt, radiusM, fallbackM) {
  const n = level.size;
  const dlat = (level.lat_max - level.lat_min) / n, dlon = (level.lon_max - level.lon_min) / n;
  const latOf = (i) => level.lat_max - (i + 0.5) * dlat;         // row 0 is north, i.e. v = 0
  const lonOf = (j) => level.lon_min + (j + 0.5) * dlon;
  const h = new Float32Array((n + 2) * (n + 2));
  for (let i = -1; i <= n; i++) {
    for (let j = -1; j <= n; j++) {
      const v = heightAt(latOf(i), lonOf(j));
      h[(i + 1) * (n + 2) + (j + 1)] = Number.isNaN(v) ? fallbackM : v;
    }
  }
  const bytes = new Uint8Array(4 * n * n);
  const dN = THREE.MathUtils.degToRad(dlat) * radiusM;
  for (let i = 0; i < n; i++) {
    const dE = THREE.MathUtils.degToRad(dlon) * radiusM * Math.cos(THREE.MathUtils.degToRad(latOf(i)));
    for (let j = 0; j < n; j++) {
      const c = (i + 1) * (n + 2) + (j + 1);
      const dhdE = (h[c + 1] - h[c - 1]) / (2 * dE);
      const dhdN = (h[c - (n + 2)] - h[c + (n + 2)]) / (2 * dN);   // row i-1 is north
      const inv = 1 / Math.sqrt(dhdE * dhdE + dhdN * dhdN + 1);
      const k = 4 * (i * n + j);
      bytes[k] = Math.round((-dhdE * inv * 0.5 + 0.5) * 255);
      bytes[k + 1] = Math.round((-dhdN * inv * 0.5 + 0.5) * 255);
      bytes[k + 2] = Math.round((inv * 0.5 + 0.5) * 255);
      bytes[k + 3] = 255;
    }
  }
  const texture = new THREE.DataTexture(bytes, n, n, THREE.RGBAFormat);
  texture.magFilter = THREE.LinearFilter;
  texture.minFilter = THREE.LinearMipmapLinearFilter;
  texture.generateMipmaps = true;
  texture.needsUpdate = true;
  return texture;
}

// The grid lookup is generated rather than looped: a sampler array indexed by a
// running variable is not portable, and there are only ever a couple of grids.
function terrainGridGlsl(count) {
  let declarations = '', lookup = '';
  for (let k = 0; k < count; k++) {
    declarations += `uniform sampler2D uGridTex${k};\nuniform vec4 uGridBox${k};\nuniform vec2 uGridRange${k};\n`;
    lookup += `
  {
    vec2 uv = vec2((lon - uGridBox${k}.z) / (uGridBox${k}.w - uGridBox${k}.z), (uGridBox${k}.y - lat) / (uGridBox${k}.y - uGridBox${k}.x));
    if (uv.x > 0.0 && uv.x < 1.0 && uv.y > 0.0 && uv.y < 1.0) {
      vec4 t = texture2D(uGridTex${k}, uv);
      return uGridRange${k}.x + (t.r * 65280.0 + t.g * 255.0) / 65535.0 * uGridRange${k}.y;
    }
  }`;
  }
  return `${declarations}
float terrainGridHeight(float lat, float lon) {${lookup}
  return -1.0e9;
}
`;
}

const TERRAIN_VERTEX = `
#include <common>
uniform mat4 uModel;
uniform vec3 uSunLocal;
attribute vec2 aLatLon;
attribute float aHeightM;
varying vec2 vUv;
varying vec2 vLatLon;
varying float vHeightM;
varying vec3 vSunENU;
varying vec3 vViewDir;
varying vec3 vUpWorld;
varying vec3 vEastWorld;
varying vec3 vNorthWorld;
#include <shadowmap_pars_vertex>
#include <logdepthbuf_pars_vertex>
void main() {
  vUv = uv;
  vLatLon = aLatLon;
  vHeightM = aHeightM;
  float latR = radians(aLatLon.x), lonR = radians(aLatLon.y);
  float cl = cos(latR), sl = sin(latR), co = cos(lonR), so = sin(lonR);
  vec3 up = vec3(cl * co, cl * so, sl);            // object space is the body-fixed frame
  vec3 east = vec3(-so, co, 0.0);
  vec3 north = vec3(-sl * co, -sl * so, cl);
  vSunENU = vec3(dot(uSunLocal, east), dot(uSunLocal, north), dot(uSunLocal, up));
  mat3 m = mat3(uModel);
  vUpWorld = m * up; vEastWorld = m * east; vNorthWorld = m * north;
  vec3 objectNormal = normal;
  vec3 transformedNormal = normalMatrix * objectNormal;
  vec4 worldPosition = modelMatrix * vec4(position, 1.0);
  vViewDir = cameraPosition - worldPosition.xyz;
  vec4 mvPosition = modelViewMatrix * vec4(position, 1.0);
  gl_Position = projectionMatrix * mvPosition;
  #include <shadowmap_vertex>
  #include <logdepthbuf_vertex>
}`;

function terrainFragment(gridCount) {
  return `
#include <common>
#include <packing>
#include <lights_pars_begin>
#include <shadowmap_pars_fragment>
#include <shadowmask_pars_fragment>
#include <logdepthbuf_pars_fragment>
uniform sampler2D uMap;
uniform sampler2D uNormalMap;
uniform vec3 uSunWorld;
uniform float uLunar;
uniform float uRefRadiusKm;
uniform vec4 uMarch;          // first step (km), step ratio, bias (m), march on/off
varying vec2 vUv;
varying vec2 vLatLon;
varying float vHeightM;
varying vec3 vSunENU;
varying vec3 vViewDir;
varying vec3 vUpWorld;
varying vec3 vEastWorld;
varying vec3 vNorthWorld;
${TERRAIN_REFLECTANCE_GLSL}
${terrainGridGlsl(gridCount)}
// 1 where the Sun is visible from this fragment, 0 where terrain up-sun hides it.
float terrainSunVisibility(float cosLat) {
  float cosZ = vSunENU.z;
  if (cosZ <= 0.0) return 0.0;
  float horiz2 = max(1.0 - cosZ * cosZ, 0.0);
  float h0 = max(vHeightM, terrainGridHeight(vLatLon.x, vLatLon.y));
  float occlusion = 0.0;
  float s = uMarch.x;
  for (int i = 0; i < ${TERRAIN_SHADOW_STEPS}; i++) {
    float hRay = h0 + 1000.0 * (s * cosZ + s * s * horiz2 / (2.0 * uRefRadiusKm)) + uMarch.z;
    float lat = vLatLon.x + degrees(s * vSunENU.y / uRefRadiusKm);
    float lon = vLatLon.y + degrees(s * vSunENU.x / (uRefRadiusKm * max(cosLat, 1.0e-4)));
    float ground = terrainGridHeight(lat, lon);
    if (ground > -1.0e8) {
      float w = max(500.0 * s * ${TERRAIN_SUN_ANGULAR_RADIUS}, 0.25);
      occlusion = max(occlusion, smoothstep(-w, w, ground - hRay));
      if (occlusion > 0.995) break;
    }
    s *= uMarch.y;
  }
  return 1.0 - occlusion;
}
void main() {
  #include <logdepthbuf_fragment>
  vec3 up = normalize(vUpWorld), east = normalize(vEastWorld), north = normalize(vNorthWorld);
  vec3 slope = texture2D(uNormalMap, vUv).xyz * 2.0 - 1.0;
  vec3 N = normalize(slope.x * east + slope.y * north + slope.z * up);
  vec3 V = normalize(vViewDir);
  float mu = max(dot(N, V), 1.0e-3);
  vec3 albedo = texture2D(uMap, vUv).rgb;
  float vehicleShadow = getShadowMask();
  float sunVisible = uMarch.w > 0.5 ? terrainSunVisibility(cos(radians(vLatLon.x))) : 1.0;
  vec3 direct = vec3(0.0);
  #if NUM_DIR_LIGHTS > 0
  for (int i = 0; i < NUM_DIR_LIGHTS; i++) {
    vec3 L = inverseTransformDirection(directionalLights[i].direction, viewMatrix);
    float mu0 = dot(N, L);
    if (mu0 <= 0.0) continue;
    // The Sun is the light the terrain horizon and the shadow map belong to; a
    // fill light from somewhere else (earthshine) is not shadowed by either.
    float shade = dot(L, uSunWorld) > 0.99 ? vehicleShadow * sunVisible : 1.0;
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

function gridHeight(g, latDeg, lonDeg) {
  let lon = lonDeg;
  while (lon < g.lon_min - 1e-12 && lon + 360 <= g.lon_max + 1e-9) lon += 360;
  while (lon > g.lon_max + 1e-12 && lon - 360 >= g.lon_min - 1e-9) lon -= 360;
  if (latDeg < g.lat_min || latDeg > g.lat_max || lon < g.lon_min || lon > g.lon_max) return NaN;
  const rows = g.rows, cols = g.cols;
  const dlat = (g.lat_max - g.lat_min) / rows, dlon = (g.lon_max - g.lon_min) / cols;
  let fr = (g.lat_max - latDeg) / dlat - 0.5, fc = (lon - g.lon_min) / dlon - 0.5;
  fr = Math.min(rows - 1, Math.max(0, fr)); fc = Math.min(cols - 1, Math.max(0, fc));
  const r0 = Math.min(rows - 2, Math.floor(fr)), c0 = Math.min(cols - 2, Math.floor(fc));
  const tr = fr - r0, tc = fc - c0;
  const h = g.data;
  const h00 = h[r0 * cols + c0], h01 = h[r0 * cols + c0 + 1], h10 = h[(r0 + 1) * cols + c0], h11 = h[(r0 + 1) * cols + c0 + 1];
  return (1 - tr) * ((1 - tc) * h00 + tc * h01) + tr * ((1 - tc) * h10 + tc * h11);
}

// Texels across a level: the resolution of the finest grid that covers its
// center, rounded to a power of two and clamped, so the normal map is as sharp
// as the DEM under it and no sharper.
function terrainNormalSize(level, grids) {
  const lat = 0.5 * (level.lat_min + level.lat_max), lon = 0.5 * (level.lon_min + level.lon_max);
  let cell = Infinity;
  for (const g of grids) {
    if (lat < g.lat_min || lat > g.lat_max || lon < g.lon_min || lon > g.lon_max) continue;
    cell = Math.min(cell, (g.lat_max - g.lat_min) / g.rows);
    break;    // grids are finest first
  }
  if (!Number.isFinite(cell) || !(cell > 0)) return TERRAIN_NORMAL_MIN;
  const want = Math.pow(2, Math.round(Math.log2((level.lat_max - level.lat_min) / cell)));
  return Math.max(TERRAIN_NORMAL_MIN, Math.min(TERRAIN_NORMAL_MAX, want));
}

export function createTerrain(spec, planet, options = {}) {
  const group = new THREE.Group();
  group.name = 'terrain';
  if (!spec || !spec.grids || spec.grids.length === 0) {
    return { group, heightAt: () => NaN, levels: [], hole: null, site: null, setSun() {}, setVisible() {} };
  }
  const R = (spec.reference_radius_m || planet.equatorial_radius_m) * TERRAIN_M_TO_KM;
  const grids = spec.grids.map((g) => ({ ...g, data: decodeFloat32(g.heights) }));
  const heightAt = (latDeg, lonDeg) => {
    for (const g of grids) { const h = gridHeight(g, latDeg, lonDeg); if (!Number.isNaN(h)) return h; }
    return NaN;
  };
  const bodyPoint = (latDeg, lonDeg, hM) => {
    const φ = THREE.MathUtils.degToRad(latDeg), λ = THREE.MathUtils.degToRad(lonDeg);
    const r = R + hM * TERRAIN_M_TO_KM;
    return new THREE.Vector3(r * Math.cos(φ) * Math.cos(λ), r * Math.cos(φ) * Math.sin(λ), r * Math.sin(φ));
  };
  const heightTextures = grids.map(terrainHeightTexture);
  const marchRatio = Math.pow(TERRAIN_SHADOW_REACH_KM / TERRAIN_SHADOW_FIRST_KM, 1 / (TERRAIN_SHADOW_STEPS - 1));
  const gridFragment = terrainFragment(grids.length);
  const levels = [];
  const imagery = (spec.imagery || []).slice().sort((a, b) => (b.lat_max - b.lat_min) - (a.lat_max - a.lat_min)); // coarse first
  // Each level is drawn only where no finer level exists (the finer box is cut
  // out of it), so a coarse chord never pokes through the fine surface; a
  // skirt hangs from every level's edges to hide the seams between them.
  const SKIRT_M = 80;
  imagery.forEach((lvl, k) => {
    const n = PATCH_SEGMENTS;
    const inner = imagery[k + 1] || null;
    const geometry = new THREE.BufferGeometry();
    const positions = [], uvs = [], index = [], latlons = [], heights = [];
    const heightOf = (lat, lon) => { const h = heightAt(lat, lon); return Number.isNaN(h) ? (spec.site ? spec.site.height_m : 0) : h; };
    const push = (lat, lon, h, u, vv) => {
      const p = bodyPoint(lat, lon, h);
      positions.push(p.x, p.y, p.z); uvs.push(u, vv); latlons.push(lat, lon); heights.push(h);
      return positions.length / 3 - 1;
    };
    const latOf = (i) => lvl.lat_max - (lvl.lat_max - lvl.lat_min) * i / n;
    const lonOf = (j) => lvl.lon_min + (lvl.lon_max - lvl.lon_min) * j / n;
    for (let i = 0; i <= n; i++) for (let j = 0; j <= n; j++) push(latOf(i), lonOf(j), heightOf(latOf(i), lonOf(j)), j / n, 1 - i / n);
    const insideInner = (lat, lon) => inner && lat > inner.lat_min && lat < inner.lat_max && lon > inner.lon_min && lon < inner.lon_max;
    for (let i = 0; i < n; i++) for (let j = 0; j < n; j++) {
      if (insideInner(0.5 * (latOf(i) + latOf(i + 1)), 0.5 * (lonOf(j) + lonOf(j + 1)))) continue;
      const a = i * (n + 1) + j, b = a + 1, c = a + (n + 1), d = c + 1;
      index.push(a, c, b, b, c, d);
    }
    // skirts along the four outer edges and, for a level with a cutout, the four inner edges
    const skirt = (pts) => {
      for (let m = 0; m + 1 < pts.length; m++) {
        const [i0, j0] = pts[m], [i1, j1] = pts[m + 1];
        const top0 = i0 * (n + 1) + j0, top1 = i1 * (n + 1) + j1;
        const b0 = push(latOf(i0), lonOf(j0), heightOf(latOf(i0), lonOf(j0)) - SKIRT_M, j0 / n, 1 - i0 / n);
        const b1 = push(latOf(i1), lonOf(j1), heightOf(latOf(i1), lonOf(j1)) - SKIRT_M, j1 / n, 1 - i1 / n);
        index.push(top0, top1, b0, top1, b1, b0, top0, b0, top1, top1, b0, b1);   // both windings: skirts are seen from either side
      }
    };
    const edge = (i0, j0, i1, j1) => { const pts = []; const steps = Math.max(Math.abs(i1 - i0), Math.abs(j1 - j0)); for (let m = 0; m <= steps; m++) pts.push([i0 + Math.sign(i1 - i0) * m, j0 + Math.sign(j1 - j0) * m]); return pts; };
    skirt(edge(0, 0, 0, n)); skirt(edge(n, 0, n, n)); skirt(edge(0, 0, n, 0)); skirt(edge(0, n, n, n));
    if (inner) {
      const iA = Math.round((lvl.lat_max - inner.lat_max) / (lvl.lat_max - lvl.lat_min) * n), iB = Math.round((lvl.lat_max - inner.lat_min) / (lvl.lat_max - lvl.lat_min) * n);
      const jA = Math.round((inner.lon_min - lvl.lon_min) / (lvl.lon_max - lvl.lon_min) * n), jB = Math.round((inner.lon_max - lvl.lon_min) / (lvl.lon_max - lvl.lon_min) * n);
      skirt(edge(iA, jA, iA, jB)); skirt(edge(iB, jA, iB, jB)); skirt(edge(iA, jA, iB, jA)); skirt(edge(iA, jB, iB, jB));
    }
    geometry.setAttribute('position', new THREE.Float32BufferAttribute(positions, 3));
    geometry.setAttribute('uv', new THREE.Float32BufferAttribute(uvs, 2));
    geometry.setAttribute('aLatLon', new THREE.Float32BufferAttribute(latlons, 2));
    geometry.setAttribute('aHeightM', new THREE.Float32BufferAttribute(heights, 1));
    geometry.setIndex(index);
    geometry.computeVertexNormals();
    geometry.computeBoundingSphere();
    // The draped image: the albedo-normalized copy when the site carries one,
    // since the page lights the ground itself, and the mosaic otherwise.
    const texture = new THREE.TextureLoader().load(lvl.albedo_url || lvl.url);
    texture.colorSpace = THREE.SRGBColorSpace;
    texture.anisotropy = options.anisotropy || 4;
    texture.generateMipmaps = true;
    texture.minFilter = THREE.LinearMipmapLinearFilter;
    const size = terrainNormalSize(lvl, grids);
    const normalMap = terrainNormalTexture({ ...lvl, size }, heightAt, (spec.reference_radius_m || planet.equatorial_radius_m), spec.site ? spec.site.height_m : 0);
    const uniforms = THREE.UniformsUtils.merge([THREE.UniformsLib.lights]);
    uniforms.uMap = { value: texture };
    uniforms.uNormalMap = { value: normalMap };
    uniforms.uSunWorld = { value: new THREE.Vector3(1, 0, 0) };
    uniforms.uSunLocal = { value: new THREE.Vector3(1, 0, 0) };
    uniforms.uLunar = { value: 0 };
    uniforms.uRefRadiusKm = { value: R };
    uniforms.uModel = { value: new THREE.Matrix4() };
    uniforms.uMarch = { value: new THREE.Vector4(TERRAIN_SHADOW_FIRST_KM, marchRatio, TERRAIN_SHADOW_BIAS_M, 0) };
    heightTextures.forEach((g, gi) => {
      uniforms[`uGridTex${gi}`] = { value: g.texture };
      uniforms[`uGridBox${gi}`] = { value: new THREE.Vector4(grids[gi].lat_min, grids[gi].lat_max, grids[gi].lon_min, grids[gi].lon_max) };
      uniforms[`uGridRange${gi}`] = { value: new THREE.Vector2(g.base, g.range) };
    });
    const material = new THREE.ShaderMaterial({
      uniforms, vertexShader: TERRAIN_VERTEX, fragmentShader: gridFragment, lights: true, side: THREE.DoubleSide,
    });
    material.map = texture;   // the exposure metering in lighting.js reads the ground texture off the material
    const mesh = new THREE.Mesh(geometry, material);
    uniforms.uModel.value = mesh.matrixWorld;
    // The path tracer uploads materials by their `color`, which a ShaderMaterial
    // has none of: hand it the plain material this one replaces.
    mesh.userData.baseMaterial = new THREE.MeshStandardMaterial({ map: texture, roughness: 1, side: THREE.DoubleSide });
    mesh.name = `terrain-level-${k}`;
    mesh.castShadow = true;
    mesh.receiveShadow = true;
    mesh.renderOrder = 1 + k;
    group.add(mesh);
    levels.push({ mesh, spec: lvl, m_per_px: lvl.m_per_px, normalSize: size, albedo: !!lvl.albedo_url });
  });
  // the globe is cut under the outermost level
  const outer = imagery[0];
  const hole = outer ? { lat_min: outer.lat_min, lat_max: outer.lat_max, lon_min: outer.lon_min, lon_max: outer.lon_max } : null;
  // site marker: a thin ring on the ground
  let siteMarker = null;
  if (spec.site) {
    const h = Number.isFinite(spec.site.height_m) ? spec.site.height_m : heightAt(spec.site.lat_deg, spec.site.lon_deg) || 0;
    const center = bodyPoint(spec.site.lat_deg, spec.site.lon_deg, h + 0.5);
    const ring = new THREE.Mesh(new THREE.RingGeometry(0.004, 0.005, 48), new THREE.MeshBasicMaterial({ color: 0xf2b950, side: THREE.DoubleSide, transparent: true, opacity: 0.8, depthWrite: false }));
    ring.position.copy(center);
    ring.lookAt(center.clone().multiplyScalar(2));
    ring.renderOrder = 20;
    ring.name = 'site-marker';
    ring.userData.uiOnly = true;
    group.add(ring);
    siteMarker = ring;
  }
  const lunar = (options.reflectance || 'lunar') === 'lunar' ? 1 : 0;
  const terrainWorldQuaternion = new THREE.Quaternion();
  const terrainSunLocal = new THREE.Vector3();
  let sunPlaced = false;
  return {
    group,
    levels,
    hole,
    site: spec.site || null,
    referenceRadiusKm: R,
    heightAt,
    setVisible(v) { group.visible = v; },
    // Scene-space unit vector toward the Sun (lighting.js `direction`), once a
    // frame. The reflectance and the horizon test both want it in the
    // body-fixed frame the grids are in, which is this group's own.
    setSun(dirScene) {
      if (!dirScene) return;
      group.getWorldQuaternion(terrainWorldQuaternion);
      terrainSunLocal.copy(dirScene).applyQuaternion(terrainWorldQuaternion.invert()).normalize();
      sunPlaced = true;
      for (const level of levels) {
        const u = level.mesh.material.uniforms;
        u.uSunWorld.value.copy(dirScene).normalize();
        u.uSunLocal.value.copy(terrainSunLocal);
        u.uLunar.value = lunar;
        u.uMarch.value.w = options.terrainShadows === false ? 0 : 1;
      }
    },
    get sunPlaced() { return sunPlaced; },
    get modelStatus() {
      // `options.sun` is the run's own sun direction: without it setSun is never
      // called and both the reflectance and the horizon test stay switched off.
      const lit = sunPlaced || options.sun;
      const shading = lit ? `${lunar ? 'lunar' : 'Lambert'} reflectance, ray-marched shadows` : 'Lambert shading';
      return `${levels.length} levels, ${grids.length} grids, finest ${levels.length ? levels[levels.length - 1].m_per_px.toFixed(2) : '–'} m/px, ${levels.some((l) => l.albedo) ? 'albedo imagery' : 'mosaic imagery'}, ${shading}`;
    },
  };
}
