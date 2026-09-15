// Atmosphere layer (phase 6): a limb glow at the entry interface, translucent
// density shells from the run's sampled profile, and, for models that vary
// horizontally, a density map draped at one altitude. Everything is built in
// the body-fixed frame and parented to the globe group so it rotates with the
// body. Units are kilometers like the rest of the scene.
import * as THREE from 'three';
import { GEOMETRY_TO_BODY } from 'viewer/globe.js';
import { inferno } from 'viewer/colormaps.js';

const LIMB_VERTEX = `
  #include <common>
  #include <logdepthbuf_pars_vertex>
  varying vec3 vNormal;
  varying vec3 vViewPos;
  void main() {
    vNormal = normalize(normalMatrix * normal);
    vec4 mvPosition = modelViewMatrix * vec4(position, 1.0);
    vViewPos = mvPosition.xyz;
    gl_Position = projectionMatrix * mvPosition;
    #include <logdepthbuf_vertex>
  }`;

const LIMB_FRAGMENT = `
  #include <common>
  #include <logdepthbuf_pars_fragment>
  uniform vec3 glowColor;
  uniform float intensity;
  uniform float power;
  varying vec3 vNormal;
  varying vec3 vViewPos;
  void main() {
    #include <logdepthbuf_fragment>
    float facing = abs(dot(normalize(vNormal), normalize(-vViewPos)));
    float rim = pow(1.0 - facing, power);
    gl_FragColor = vec4(glowColor * rim * intensity, rim * intensity);
  }`;

const GLOW_COLORS = { earth: 0x7fb8ff, mars: 0xf0b48c, venus: 0xf6dfa8, titan: 0xe6c27a, moon: 0xffffff };

function oblateSphere(Re, Rp, altitudeKm, segments = 96) {
  const geometry = new THREE.SphereGeometry(1, segments, segments / 2);
  return { geometry, scale: new THREE.Vector3(Re + altitudeKm, Rp + altitudeKm, Re + altitudeKm) };
}

// Interpolate log-density from the sampled profile at an altitude in meters.
function profileDensity(profile, hM) {
  const h = profile.altitude_m, d = profile.density_kg_m3;
  if (!h || h.length < 2) return NaN;
  if (hM <= h[0]) return d[0];
  if (hM >= h[h.length - 1]) return d[d.length - 1];
  let lo = 0, hi = h.length - 1;
  while (hi - lo > 1) { const mid = (lo + hi) >> 1; if (h[mid] <= hM) lo = mid; else hi = mid; }
  const f = (hM - h[lo]) / (h[hi] - h[lo]);
  const a = Math.log(Math.max(d[lo], 1e-30)), b = Math.log(Math.max(d[hi], 1e-30));
  return Math.exp(a + f * (b - a));
}

export function createAtmosphere(spec, planet, options = {}) {
  const Re = planet.equatorial_radius_m / 1000, Rp = planet.polar_radius_m / 1000;
  const eiKm = spec.ei_altitude_m / 1000;
  const key = (planet.texture || planet.name || '').toLowerCase();
  const color = new THREE.Color(GLOW_COLORS[key] ?? 0x9fc5ff);
  const group = new THREE.Group();
  group.name = 'atmosphere';
  group.quaternion.copy(GEOMETRY_TO_BODY);

  // Limb glow at the entry interface (both sides, so it reads from inside too).
  const limb = oblateSphere(Re, Rp, eiKm);
  const limbMesh = new THREE.Mesh(limb.geometry, new THREE.ShaderMaterial({
    uniforms: { glowColor: { value: color }, intensity: { value: options.intensity ?? 0.9 }, power: { value: 2.2 } },
    vertexShader: LIMB_VERTEX,
    fragmentShader: LIMB_FRAGMENT,
    transparent: true,
    depthWrite: false,
    blending: THREE.AdditiveBlending,
    side: THREE.DoubleSide,
  }));
  limbMesh.scale.copy(limb.scale);
  limbMesh.renderOrder = 2;
  group.add(limbMesh);

  // Density shells between 0.3 and 1.0 of the EI, opacity from the profile.
  const layers = new THREE.Group();
  const layerInfo = [];
  const profile = spec.profile || {};
  const hasProfile = profile.altitude_m && profile.altitude_m.length > 1;
  if (hasProfile) {
    const rhoEi = profileDensity(profile, spec.ei_altitude_m);
    const rhoLow = profileDensity(profile, 0.3 * spec.ei_altitude_m);
    const span = Math.log(Math.max(rhoLow, 1e-30)) - Math.log(Math.max(rhoEi, 1e-30));
    const K = options.layers ?? 6;
    for (let k = 0; k < K; k++) {
      const hM = spec.ei_altitude_m * (0.3 + 0.7 * k / (K - 1));
      const rho = profileDensity(profile, hM);
      const a = span > 0 ? (Math.log(Math.max(rho, 1e-30)) - Math.log(Math.max(rhoEi, 1e-30))) / span : 0;
      const opacity = 0.02 + 0.16 * Math.min(1, Math.max(0, a));
      const shell = oblateSphere(Re, Rp, hM / 1000, 64);
      const mesh = new THREE.Mesh(shell.geometry, new THREE.MeshBasicMaterial({
        color, transparent: true, opacity, depthWrite: false, blending: THREE.AdditiveBlending, side: THREE.FrontSide,
      }));
      mesh.scale.copy(shell.scale);
      mesh.renderOrder = 1;
      layers.add(mesh);
      layerInfo.push({ altitude_km: hM / 1000, density: rho, opacity });
    }
  }
  layers.visible = options.layersVisible ?? true;
  group.add(layers);

  // Density map at one altitude (lat/lon grid -> inferno texture on a shell).
  let mapMesh = null;
  let mapInfo = null;
  const map = spec.map;
  if (map && map.density_kg_m3 && map.density_kg_m3.length === map.lat_deg.length * map.lon_deg.length) {
    const W = map.lon_deg.length, H = map.lat_deg.length;
    const vals = map.density_kg_m3.map((v) => Math.log10(Math.max(v, 1e-30)));
    const finite = vals.filter(Number.isFinite);
    const lo = Math.min(...finite), hi = Math.max(...finite);
    const canvas = document.createElement('canvas');
    canvas.width = W; canvas.height = H;
    const ctx = canvas.getContext('2d');
    const img = ctx.createImageData(W, H);
    for (let i = 0; i < H; i++) {
      const row = H - 1 - i; // canvas row 0 is north; lat_deg is ascending
      for (let j = 0; j < W; j++) {
        const t = hi > lo ? (vals[i * W + j] - lo) / (hi - lo) : 0.5;
        const c = inferno(t);
        const o = 4 * (row * W + j);
        img.data[o] = Math.round(255 * c[0]); img.data[o + 1] = Math.round(255 * c[1]); img.data[o + 2] = Math.round(255 * c[2]); img.data[o + 3] = 255;
      }
    }
    ctx.putImageData(img, 0, 0);
    const texture = new THREE.CanvasTexture(canvas);
    texture.colorSpace = THREE.SRGBColorSpace;
    texture.magFilter = THREE.LinearFilter;
    texture.wrapS = THREE.RepeatWrapping;
    texture.offset.x = 0; // grid runs 180 W .. 180 E like the default sphere mapping
    const shell = oblateSphere(Re, Rp, map.altitude_m / 1000, 96);
    mapMesh = new THREE.Mesh(shell.geometry, new THREE.MeshBasicMaterial({ map: texture, transparent: true, opacity: 0.6, depthWrite: false }));
    mapMesh.scale.copy(shell.scale);
    mapMesh.renderOrder = 1;
    mapMesh.visible = options.mapVisible ?? true;
    group.add(mapMesh);
    mapInfo = { altitude_km: map.altitude_m / 1000, min: Math.pow(10, lo), max: Math.pow(10, hi), varies: hi - lo > 0.02 };
  }

  return {
    group,
    limb: limbMesh,
    layers,
    map: mapMesh,
    info: {
      model: spec.model,
      ei_km: eiKm,
      profile: hasProfile ? { points: profile.altitude_m.length, surface: profile.density_kg_m3[0], atEi: profileDensity(profile, spec.ei_altitude_m) } : null,
      layers: layerInfo,
      map: mapInfo,
    },
    densityAtAltitude(hM) { return hasProfile ? profileDensity(profile, hM) : NaN; },
    setVisible(v) { group.visible = v; },
    setLayersVisible(v) { layers.visible = v; },
    setMapVisible(v) { if (mapMesh) mapMesh.visible = v; },
    setLimbVisible(v) { limbMesh.visible = v; },
  };
}
