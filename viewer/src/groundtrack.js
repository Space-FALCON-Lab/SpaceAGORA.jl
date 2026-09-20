// Ground tracks: the sub-satellite point of every spacecraft drawn on the
// body's own surface.
//
// Frames. The track is a body-fixed curve, so it is built once in the
// planet's body-fixed frame and added to the globe group, which already
// carries q_pi(t). It is then correct in both viewer frames without any
// per-frame rotation: in the inertial view the track turns with the globe
// under the spacecraft, and in the planet-fixed view it stands still while
// the orbit precesses over it.
//
// Each sample is the geodetic sub-satellite point (viewer/data.js `geodetic`
// inverted onto the ellipsoid), lifted by a small fraction of the equatorial
// radius so the line does not z-fight with the globe mesh or disappear into
// a terrain patch. The history up to the current time is drawn solid and a
// marker rides the head; the rest of the pass is left undrawn so the track
// reads as something being laid down rather than a static map.
import * as THREE from 'three';

import { geodetic, rotateByConjugate } from 'viewer/data.js';
import { spacecraftColor } from 'viewer/spacecraft.js';

// The track sits this fraction of the equatorial radius above the ellipsoid.
const GROUNDTRACK_LIFT_FRACTION = 0.0015;
// Samples per spacecraft. A 96 h, eight-spacecraft run decimates to a few
// thousand frames; more points than this buy nothing at globe scale.
const GROUNDTRACK_MAX_POINTS = 4000;
const GROUNDTRACK_MARKER_PIXELS = 5;

// Body-fixed point on the ellipsoid surface below `rBody` (kilometers),
// lifted by `lift`. Geodetic latitude, so the track sits under the vehicle
// rather than under its geocentric direction (the two differ by up to about
// 0.19 degrees of latitude on Earth).
function groundtrackSurfacePoint(rBody, Re, Rp, lift, out) {
  const g = geodetic(rBody, Re, Rp);
  const lat = g.latDeg * Math.PI / 180, lon = g.lonDeg * Math.PI / 180;
  const e2 = 1 - (Rp * Rp) / (Re * Re);
  const s = Math.sin(lat), c = Math.cos(lat);
  const N = Re / Math.sqrt(1 - e2 * s * s);
  out.set((N + lift) * c * Math.cos(lon), (N + lift) * c * Math.sin(lon), (N * (1 - e2) + lift) * s);
  return out;
}

/**
 * createGroundTracks(globe, frames, sidecar, options) -> { group, update(t), setVisible(v), visible }
 *
 * `globe` supplies the rotation table and the ellipsoid radii; `frames` the
 * inertial positions; `sidecar` the spacecraft list (for the count only —
 * the colors match the markers in viewer/spacecraft.js).
 */
export function createGroundTracks(globe, frames, sidecar, options = {}) {
  const group = new THREE.Group();
  group.name = 'ground-tracks';
  const S = frames.sats;
  const Re = options.equatorialRadiusKm ?? globe.radiusKm;
  const Rp = options.polarRadiusKm ?? globe.radiusKm;
  const lift = GROUNDTRACK_LIFT_FRACTION * Re;

  const stride = Math.max(1, Math.ceil(frames.count / GROUNDTRACK_MAX_POINTS));
  const times = [];
  for (let k = 0; k < frames.count; k += stride) times.push(frames.t[k]);
  const n = times.length;

  // Body-fixed surface points, one pass over the frames.
  const tracks = [];
  const rIn = new Float64Array(3), rBody = new Float64Array(3);
  const q = new Float32Array(4);
  const p = new THREE.Vector3();
  for (let s = 0; s < S; s++) {
    const pts = new Float32Array(n * 3);
    const valid = new Uint8Array(n);
    for (let i = 0; i < n; i++) {
      frames.positionAt(times[i], s, rIn);
      if (!Number.isFinite(rIn[0])) continue;
      globe.rotationAt(times[i], q);
      rotateByConjugate(q, rIn, rBody);
      groundtrackSurfacePoint(rBody, Re, Rp, lift, p);
      pts[3 * i] = p.x; pts[3 * i + 1] = p.y; pts[3 * i + 2] = p.z;
      valid[i] = 1;
    }
    const color = options.colorFor ? options.colorFor(s, S) : spacecraftColor(s, S);
    const geometry = new THREE.BufferGeometry();
    geometry.setAttribute('position', new THREE.BufferAttribute(pts, 3));
    geometry.setDrawRange(0, 0);
    const line = new THREE.Line(geometry, new THREE.LineBasicMaterial({ color, transparent: true, opacity: 0.85 }));
    line.frustumCulled = false;
    group.add(line);

    const headGeometry = new THREE.BufferGeometry();
    headGeometry.setAttribute('position', new THREE.BufferAttribute(new Float32Array(3), 3));
    const head = new THREE.Points(headGeometry, new THREE.PointsMaterial({
      color, size: GROUNDTRACK_MARKER_PIXELS * (options.pixelRatio ?? 1), sizeAttenuation: false, transparent: true, opacity: 0.95,
    }));
    head.frustumCulled = false;
    group.add(head);

    tracks.push({ pts, valid, line, head });
  }

  let visible = options.enabled ?? true;
  group.visible = visible;

  return {
    group,
    count: S,
    // Draw the track up to elapsed time `t`.
    update(t) {
      if (!visible) return;
      let hi = 0;
      while (hi < n && times[hi] <= t) hi++;
      for (const tr of tracks) {
        // A leading run of frames where the spacecraft has not started yet
        // (NaN position) must not be drawn as a line to the origin.
        let lo = 0;
        while (lo < hi && !tr.valid[lo]) lo++;
        let end = hi;
        while (end > lo && !tr.valid[end - 1]) end--;
        const drawn = Math.max(0, end - lo);
        tr.line.geometry.setDrawRange(lo, drawn);
        tr.line.visible = drawn >= 2;
        if (drawn >= 1) {
          const a = tr.head.geometry.getAttribute('position').array;
          a[0] = tr.pts[3 * (end - 1)]; a[1] = tr.pts[3 * (end - 1) + 1]; a[2] = tr.pts[3 * (end - 1) + 2];
          tr.head.geometry.getAttribute('position').needsUpdate = true;
        }
        tr.head.visible = drawn >= 1;
      }
    },
    setVisible(v) { visible = !!v; group.visible = visible; },
    get visible() { return visible; },
  };
}
