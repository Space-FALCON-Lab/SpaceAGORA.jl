// Reference polylines beside the flown trajectories: a planned RPO path in a
// target's RTN frame, a body-frame path, or an inertial one. RTN and body
// paths are re-expressed every frame from the target's state.
import * as THREE from 'three';
import { decodeFloat32 } from 'viewer/data.js';

export function createPaths(specs, frames, options = {}) {
  const group = new THREE.Group();
  group.name = 'reference-paths';
  const items = [];
  for (const spec of specs || []) {
    const pts = decodeFloat32(spec.points_km);
    const n = spec.count || pts.length / 3;
    if (n < 2) continue;
    const geometry = new THREE.BufferGeometry();
    geometry.setAttribute('position', new THREE.BufferAttribute(new Float32Array(n * 3), 3));
    const material = spec.dashed
      ? new THREE.LineDashedMaterial({ color: new THREE.Color(spec.color || '#7fe0ff'), dashSize: options.dashKm ?? 0.002, gapSize: options.dashKm ?? 0.002, transparent: true, opacity: 0.9 })
      : new THREE.LineBasicMaterial({ color: new THREE.Color(spec.color || '#7fe0ff'), transparent: true, opacity: 0.9 });
    const line = new THREE.Line(geometry, material);
    line.frustumCulled = false;
    group.add(line);
    items.push({ spec, pts, n, line, target: Math.max(0, (spec.target || 1) - 1) });
  }
  const r = new Float64Array(3), v = new Float64Array(3), q = new Float32Array(4);
  const R = new THREE.Vector3(), T = new THREE.Vector3(), N = new THREE.Vector3(), p = new THREE.Vector3(), quat = new THREE.Quaternion();
  return {
    group,
    items,
    // anchorKm: floating origin the group sits at (same as the marker group)
    update(t, anchorKm) {
      group.position.set(anchorKm[0], anchorKm[1], anchorKm[2]);
      for (const it of items) {
        const arr = it.line.geometry.getAttribute('position').array;
        frames.positionAt(t, it.target, r);
        const present = Number.isFinite(r[0]);
        it.line.visible = present;
        if (!present) continue;
        if (it.spec.frame === 'rtn') {
          frames.velocityAt(t, it.target, v);
          R.set(r[0], r[1], r[2]).normalize();
          N.set(r[0], r[1], r[2]).cross(new THREE.Vector3(v[0], v[1], v[2])).normalize();
          T.copy(N).cross(R);
        } else if (it.spec.frame === 'body') {
          if (frames.quaternionAt(t, it.target, q)) quat.set(q[0], q[1], q[2], q[3]); else quat.identity();
        }
        for (let k = 0; k < it.n; k++) {
          p.set(it.pts[3 * k], it.pts[3 * k + 1], it.pts[3 * k + 2]);
          if (it.spec.frame === 'rtn') {
            const x = p.x, y = p.y, z = p.z;
            p.set(R.x * x + T.x * y + N.x * z, R.y * x + T.y * y + N.y * z, R.z * x + T.z * y + N.z * z);
          } else if (it.spec.frame === 'body') {
            p.applyQuaternion(quat);
          }
          arr[3 * k] = r[0] + p.x - anchorKm[0];
          arr[3 * k + 1] = r[1] + p.y - anchorKm[1];
          arr[3 * k + 2] = r[2] + p.z - anchorKm[2];
        }
        it.line.geometry.getAttribute('position').needsUpdate = true;
        if (it.spec.dashed) it.line.computeLineDistances();
      }
    },
    setVisible(v) { group.visible = v; },
  };
}
