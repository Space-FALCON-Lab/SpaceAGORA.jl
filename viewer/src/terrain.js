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
//   terrain.grids[]            { name, rows, cols, lat_min, lat_max, lon_min, lon_max, and the heights, row-major north to south: `heights_i16` (base64 Int16 steps of `height_scale_m` about `height_base_m`) or `heights` (base64 Float32) }, finest first
//   terrain.imagery[]          { lat_min, lat_max, lon_min, lon_max, width, height, m_per_px, url (data: JPEG), albedo_url (data: JPEG, optional) }, coarse to fine
// `albedo_url` is `url` with the mosaic's own illumination divided out
// (scripts/dev/terrain/fetch_moon_site.py `derive_albedo`); it is what the page
// drapes when it is there, since the page lights the surface itself.
import * as THREE from 'three';
import { decodeBytes, decodeFloat32 } from 'viewer/data.js';

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
// Normal map: one texture per imagery level, sampled at half the cell of the
// finest height grid that covers it, so the relief the DEM resolves shades
// correctly however coarse the drawn mesh is. Over the innermost level that is
// about 1.9 m a texel against the NAC grid's 4 m and the drawn mesh's 5 m.
const TERRAIN_NORMAL_MIN = 64;
const TERRAIN_NORMAL_MAX = 512;
const TERRAIN_NORMAL_OVERSAMPLE = 2;        // texels per DEM cell
const TERRAIN_NORMAL_TEXEL_MIN_M = 1.0;     // never finer than this from the DEM alone

// Micro-relief. A NAC digital terrain model stops at a few meters -- its own
// grid is 4 m here -- but the surface does not: regolith keeps roughening down
// through decimeter craters to the grain, and at the 10 degree sun of a landing
// that roughness is most of what the eye reads as ground. Below the DEM's cell
// the page therefore adds a procedural detail normal: a tiled field of
// multi-octave value noise and small craters, built once, sampled by the finest
// levels only and faded out by its own mipmaps as the camera pulls away. It is
// a texture, not a claim about this particular site; `options.microRelief =
// false` switches it off. The slope it adds (0.12 RMS, about 7 degrees at a
// meter) continues the measured slope-versus-baseline trend of the site's own
// grid (10.9 degrees RMS at 4 m, 8.0 at 8 m, 6.0 at 16 m).
//
// A detail normal map cannot simply be minified. Mip filtering averages the
// normals, but the shading term they feed is not linear in the normal and it
// clamps at the terminator, so the average of the shading is not the shading of
// the average: at a 10 degree sun a slope of 0.12 RMS against a terminator at
// tan(10.6 deg) = 0.187 tips a large tail of texels past it, and once a pixel
// covers many texels -- which it does a few tens of meters out at the grazing
// view of a landing, where anisotropic filtering under-samples the stretched
// axis -- the far field reads as speckle rather than as grain. The shader therefore keeps only
// the part of the detail slope the pixel footprint resolves and carries the
// rest into the shading response, the treatment of Toksvig, "Mipmapping Normal
// Maps" (NVIDIA, 2004) and of the LEAN/CLEAN mapping family after it. See
// terrainFragment below for the two lines that do it.
const TERRAIN_DETAIL_SIZE = 256;            // texels across one tile
const TERRAIN_DETAIL_TILE_M = 8.0;          // meters across one tile
const TERRAIN_DETAIL_SLOPE = 0.12;          // RMS slope of the detail field
const TERRAIN_DETAIL_CRATERS = 140;         // craters per tile
const TERRAIN_DETAIL_MAX_MPP = 2.0;         // levels at least this sharp get it
const TERRAIN_DETAIL_FOOTPRINT = 2.0;       // pixel footprint, in detail texels, that costs half the slope variance

// Heights as the bundler wrote them: Int16 steps about a base (half the bytes,
// still a centimeter) or plain Float32.
function terrainHeights(grid) {
  if (grid.heights_i16 != null) {
    const bytes = decodeBytes(grid.heights_i16);
    const q = new Int16Array(bytes.buffer, 0, bytes.length >> 1);
    const out = new Float32Array(q.length);
    const base = grid.height_base_m || 0, scale = grid.height_scale_m || 1;
    for (let k = 0; k < q.length; k++) out[k] = base + q[k] * scale;
    return out;
  }
  return decodeFloat32(grid.heights);
}

// A seamless tile of micro-relief, as a tangent-space normal map: a few octaves
// of value noise (wavelengths from half the tile down to four texels)
// plus small craters -- a parabolic bowl inside a raised rim, the shape a
// simple impact leaves -- and the whole field scaled to TERRAIN_DETAIL_SLOPE
// RMS slope. Everything wraps, so the tile repeats without a seam.
function terrainDetailTexture(anisotropy) {
  const n = TERRAIN_DETAIL_SIZE, texelM = TERRAIN_DETAIL_TILE_M / n;
  const h = new Float32Array(n * n);
  let seed = 20250916;
  const rand = () => { seed = (seed * 1664525 + 1013904223) >>> 0; return seed / 4294967296; };
  // value noise: a lattice of random values per octave, bilinear, wrapped
  for (let octave = 2; octave <= 64; octave *= 2) {
    const m = octave, amp = TERRAIN_DETAIL_TILE_M / octave;   // 1/f: amplitude follows wavelength
    const lattice = new Float32Array(m * m);
    for (let k = 0; k < m * m; k++) lattice[k] = rand() - 0.5;
    for (let i = 0; i < n; i++) {
      const fy = i / n * m, y0 = Math.floor(fy), ty = fy - y0;
      for (let j = 0; j < n; j++) {
        const fx = j / n * m, x0 = Math.floor(fx), tx = fx - x0;
        const x1 = (x0 + 1) % m, y1 = (y0 + 1) % m;
        const sx = tx * tx * (3 - 2 * tx), sy = ty * ty * (3 - 2 * ty);
        const a = lattice[y0 % m * m + x0 % m], b = lattice[y0 % m * m + x1];
        const c = lattice[y1 * m + x0 % m], d = lattice[y1 * m + x1];
        h[i * n + j] += amp * ((1 - sy) * ((1 - sx) * a + sx * b) + sy * ((1 - sx) * c + sx * d));
      }
    }
  }
  // craters: radii power-law distributed between two texels and a tenth of the tile
  for (let c = 0; c < TERRAIN_DETAIL_CRATERS; c++) {
    const rM = 2 * texelM * Math.pow(0.1 * TERRAIN_DETAIL_TILE_M / (2 * texelM), Math.pow(rand(), 2.2));
    const r = rM / texelM, depth = 0.18 * rM, cx = rand() * n, cy = rand() * n;
    const reach = Math.ceil(1.5 * r);
    for (let di = -reach; di <= reach; di++) {
      for (let dj = -reach; dj <= reach; dj++) {
        const d = Math.hypot(di + (cy - Math.floor(cy)), dj + (cx - Math.floor(cx))) / r;
        if (d > 1.5) continue;
        const i = (Math.floor(cy) + di + n) % n, j = (Math.floor(cx) + dj + n) % n;
        // bowl inside the rim, rim crest at d = 1, skirt out to 1.5
        h[i * n + j] += d <= 1 ? depth * (d * d - 0.75) : 0.25 * depth * (1.5 - d) / 0.5;
      }
    }
  }
  // slopes, scaled to the wanted RMS, encoded the usual way
  const sx = new Float32Array(n * n), sy = new Float32Array(n * n);
  let sum = 0;
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) {
      const e = h[i * n + (j + 1) % n], w = h[i * n + (j - 1 + n) % n];
      const s = h[((i + 1) % n) * n + j], no = h[((i - 1 + n) % n) * n + j];
      const dE = (e - w) / (2 * texelM), dN = (no - s) / (2 * texelM);   // row i-1 is north
      sx[i * n + j] = dE; sy[i * n + j] = dN;
      sum += dE * dE + dN * dN;
    }
  }
  const rms = Math.sqrt(sum / (n * n));
  const gain = rms > 0 ? TERRAIN_DETAIL_SLOPE * Math.SQRT2 / rms : 0;
  const bytes = new Uint8Array(4 * n * n);
  for (let k = 0; k < n * n; k++) {
    const dE = sx[k] * gain, dN = sy[k] * gain;
    const inv = 1 / Math.sqrt(dE * dE + dN * dN + 1);
    bytes[4 * k] = Math.round((-dE * inv * 0.5 + 0.5) * 255);
    bytes[4 * k + 1] = Math.round((-dN * inv * 0.5 + 0.5) * 255);
    bytes[4 * k + 2] = Math.round((inv * 0.5 + 0.5) * 255);
    bytes[4 * k + 3] = 255;
  }
  const texture = new THREE.DataTexture(bytes, n, n, THREE.RGBAFormat);
  texture.wrapS = THREE.RepeatWrapping;
  texture.wrapT = THREE.RepeatWrapping;
  texture.magFilter = THREE.LinearFilter;
  texture.minFilter = THREE.LinearMipmapLinearFilter;
  texture.anisotropy = anisotropy || 4;
  texture.generateMipmaps = true;
  texture.needsUpdate = true;
  return texture;
}

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
function terrainNormalTexture(level, heightAt, radiusM, fallbackM, anisotropy) {
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
  texture.anisotropy = anisotropy || 4;
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
uniform sampler2D uDetailMap;
uniform vec3 uDetail;         // tiles across this level, the strength of their slope, its RMS
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
  // Micro-relief rides on top as a slope, added to the slope the DEM gives
  // before the two are turned back into a normal -- but only the part of it
  // this pixel resolves. Toksvig, "Mipmapping Normal Maps" (NVIDIA, 2004):
  // averaging a normal map shortens the mean normal, and what a filtered
  // normal loses in length has to come back as a wider shading response, or
  // the non-linear, terminator-clamped term it feeds turns to noise under
  // minification. fp is the footprint in detail texels along its major axis
  // -- the axis anisotropic filtering under-samples, and the one the grazing
  // view of a landing stretches -- and Toksvig's mean-normal factor
  // 1/sqrt(1 + (fp/f0)^2) is how much slope survives it. The remainder leaves
  // as detailRms, the RMS slope this fragment could not resolve, which the
  // direct term below spends as terminator width instead of as geometry.
  float detailRms = 0.0;
  if (uDetail.y > 0.0) {
    vec2 duv = vUv * uDetail.x;
    float fp = max(length(dFdx(duv)), length(dFdy(duv))) * ${TERRAIN_DETAIL_SIZE}.0 / ${TERRAIN_DETAIL_FOOTPRINT.toFixed(1)};
    float keep = inversesqrt(1.0 + fp * fp);
    vec3 grain = texture2D(uDetailMap, duv).xyz * 2.0 - 1.0;
    vec2 ds = -grain.xy / max(grain.z, 0.2) * (uDetail.y * keep);
    slope += slope.z * vec3(-ds.x, -ds.y, 0.0);
    detailRms = uDetail.z * sqrt(max(1.0 - keep * keep, 0.0));
  }
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
    // The slope the footprint could not resolve widens the terminator rather
    // than vanishing. A Gaussian slope of RMS detailRms moves cos(i) by a
    // Gaussian of width w = detailRms * |L tangential|, and the mean of the
    // clamped cosine over it is c*Phi(c/w) + w*phi(c/w), which the square root
    // (c + sqrt(c^2 + 2 w^2 / pi)) / 2 follows to within about five percent --
    // exact at the terminator, asymptotic to max(c, 0) well away from it, and
    // an error function cheaper. Shadowing between the facets is not modeled.
    float w = detailRms * length(L - dot(L, N) * N);
    mu0 = w > 0.0 ? 0.5 * (mu0 + sqrt(mu0 * mu0 + 0.63662 * w * w)) : max(mu0, 0.0);
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

// Texels across a level: half the cell of the finest grid that covers its
// center (the extra factor is what lets the interpolated field, and the
// micro-relief on top of it, carry detail between DEM samples), rounded to a
// power of two and clamped.
function terrainNormalSize(level, grids, radiusM) {
  const lat = 0.5 * (level.lat_min + level.lat_max), lon = 0.5 * (level.lon_min + level.lon_max);
  let cell = Infinity;
  for (const g of grids) {
    if (lat < g.lat_min || lat > g.lat_max || lon < g.lon_min || lon > g.lon_max) continue;
    cell = Math.min(cell, (g.lat_max - g.lat_min) / g.rows);
    break;    // grids are finest first
  }
  if (!Number.isFinite(cell) || !(cell > 0)) return TERRAIN_NORMAL_MIN;
  const spanM = THREE.MathUtils.degToRad(level.lat_max - level.lat_min) * radiusM;
  const cellM = THREE.MathUtils.degToRad(cell) * radiusM;
  const texelM = Math.max(cellM / TERRAIN_NORMAL_OVERSAMPLE, TERRAIN_NORMAL_TEXEL_MIN_M);
  const want = Math.pow(2, Math.round(Math.log2(spanM / texelM)));
  return Math.max(TERRAIN_NORMAL_MIN, Math.min(TERRAIN_NORMAL_MAX, want));
}

export function createTerrain(spec, planet, options = {}) {
  const group = new THREE.Group();
  group.name = 'terrain';
  if (!spec || !spec.grids || spec.grids.length === 0) {
    return { group, heightAt: () => NaN, levels: [], hole: null, site: null, setSun() {}, setVisible() {} };
  }
  const R = (spec.reference_radius_m || planet.equatorial_radius_m) * TERRAIN_M_TO_KM;
  const grids = spec.grids.map((g) => ({ ...g, data: terrainHeights(g) }));
  const heightAt = (latDeg, lonDeg) => {
    for (const g of grids) { const h = gridHeight(g, latDeg, lonDeg); if (!Number.isNaN(h)) return h; }
    return NaN;
  };
  const bodyPoint = (latDeg, lonDeg, hM) => {
    const φ = THREE.MathUtils.degToRad(latDeg), λ = THREE.MathUtils.degToRad(lonDeg);
    const r = R + hM * TERRAIN_M_TO_KM;
    return new THREE.Vector3(r * Math.cos(φ) * Math.cos(λ), r * Math.cos(φ) * Math.sin(λ), r * Math.sin(φ));
  };
  const radiusM = spec.reference_radius_m || planet.equatorial_radius_m;
  const heightTextures = grids.map(terrainHeightTexture);
  const detailMap = options.microRelief === false ? null : terrainDetailTexture(options.anisotropy);
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
    const size = terrainNormalSize(lvl, grids, radiusM);
    const normalMap = terrainNormalTexture({ ...lvl, size }, heightAt, radiusM, spec.site ? spec.site.height_m : 0, options.anisotropy);
    // Micro-relief on the levels sharp enough to show it, one tile every
    // TERRAIN_DETAIL_TILE_M of ground.
    const spanM = THREE.MathUtils.degToRad(lvl.lat_max - lvl.lat_min) * radiusM;
    const detail = detailMap && lvl.m_per_px <= TERRAIN_DETAIL_MAX_MPP
      ? new THREE.Vector3(spanM / TERRAIN_DETAIL_TILE_M, 1, TERRAIN_DETAIL_SLOPE)
      : new THREE.Vector3(1, 0, 0);
    const uniforms = THREE.UniformsUtils.merge([THREE.UniformsLib.lights]);
    uniforms.uMap = { value: texture };
    uniforms.uNormalMap = { value: normalMap };
    uniforms.uDetailMap = { value: detailMap || normalMap };
    uniforms.uDetail = { value: detail };
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
    material.extensions = { derivatives: true };   // the micro-relief footprint; core on WebGL2, an extension below it
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
    levels.push({ mesh, spec: lvl, m_per_px: lvl.m_per_px, normalSize: size, normalTexelM: spanM / size, microRelief: detail.y > 0, albedo: !!lvl.albedo_url });
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
      const finest = levels.length ? levels[levels.length - 1] : null;
      const relief = finest ? `normals ${finest.normalTexelM.toFixed(1)} m/texel${finest.microRelief ? ' + micro-relief' : ''}` : 'no normals';
      return `${levels.length} levels, ${grids.length} grids, finest ${finest ? finest.m_per_px.toFixed(2) : '–'} m/px, ${levels.some((l) => l.albedo) ? 'albedo imagery' : 'mosaic imagery'}, ${relief}, ${shading}`;
    },
  };
}
