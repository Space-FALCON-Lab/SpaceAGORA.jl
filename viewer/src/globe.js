// Central body: an oblate textured sphere driven by the sampled rotation table.
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

// `textureEntry` is `{ url, lon_left_deg }` from the bundler (or null). The
// sphere's default mapping puts u = 0.5 on longitude 0, i.e. an image whose
// left edge is 180 W; an image with its left edge at L degrees east is
// shifted by -(L + 180) / 360 with wrap-around.
export function createGlobe(planet, textureEntry, options = {}) {
  const Re = planet.equatorial_radius_m / 1000, Rp = planet.polar_radius_m / 1000;
  const group = new THREE.Group();
  group.name = 'globe';

  const geometry = new THREE.SphereGeometry(1, 128, 64);
  const material = new THREE.MeshLambertMaterial({ color: 0xffffff });
  const key = (planet.texture || planet.name || '').toLowerCase();
  const textureDataUrl = textureEntry ? (typeof textureEntry === 'string' ? textureEntry : textureEntry.url) : null;
  if (textureDataUrl) {
    const texture = loadTextureWithinLimit(textureDataUrl, options.maxTextureSize || 4096);
    texture.colorSpace = THREE.SRGBColorSpace;
    texture.anisotropy = options.anisotropy || 4;
    texture.generateMipmaps = true;
    texture.minFilter = THREE.LinearMipmapLinearFilter;
    const lonLeft = (textureEntry && typeof textureEntry === 'object' && Number.isFinite(textureEntry.lon_left_deg)) ? textureEntry.lon_left_deg : -180;
    texture.wrapS = THREE.RepeatWrapping;
    texture.offset.x = -(lonLeft + 180) / 360;
    material.map = texture;
  } else {
    material.color.setHex(FALLBACK_COLORS[key] ?? 0x888888);
  }
  // A terrain patch replaces the sphere inside `options.hole` (a lat/lon box,
  // degrees): the fragment shader discards the sphere there. Geometry (x, y, z)
  // maps to body (x, -z, y).
  if (options.hole) {
    const h = options.hole;
    material.onBeforeCompile = (shader) => {
      shader.uniforms.uHole = { value: new THREE.Vector4(h.lat_min, h.lat_max, h.lon_min, h.lon_max) };
      shader.vertexShader = shader.vertexShader
        .replace('#include <common>', '#include <common>\nvarying vec3 vBodyDir;')
        .replace('#include <begin_vertex>', '#include <begin_vertex>\nvBodyDir = normalize(vec3(position.x, -position.z, position.y));');
      shader.fragmentShader = shader.fragmentShader
        .replace('#include <common>', '#include <common>\nuniform vec4 uHole;\nvarying vec3 vBodyDir;')
        .replace('#include <clipping_planes_fragment>', `#include <clipping_planes_fragment>
  {
    float latDeg = degrees(asin(clamp(vBodyDir.z, -1.0, 1.0)));
    float lonDeg = degrees(atan(vBodyDir.y, vBodyDir.x));
    float lonMin = uHole.z, lonMax = uHole.w;
    if (lonDeg < lonMin - 180.0) lonDeg += 360.0;
    if (lonDeg > lonMax + 180.0) lonDeg -= 360.0;
    if (latDeg >= uHole.x && latDeg <= uHole.y && lonDeg >= lonMin && lonDeg <= lonMax) discard;
  }`);
    };
  }
  const mesh = new THREE.Mesh(geometry, material);
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

  return {
    group,
    mesh,
    grid,
    axis,
    radiusKm: Re,
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
