// Reference trajectories: a "ghost" of a spacecraft driven by a state table
// that did not come from the integrator (a SPICE reconstruction, a
// telemetry record, a planned trajectory), drawn beside the flown one so the
// two can be compared on the same timeline. Each ghost copies the geometry
// of the spacecraft it refers to (the 3D model override or the link boxes),
// rendered translucent, with its own full-history line, a ring marker and a
// label. Positions and times are Float64 with their own time grid; the
// group sits at the floating origin like the assemblies.
import * as THREE from 'three';
import { decodeFloat64, decodeFloat32, slerp, velocityAlignedQuaternion } from 'viewer/data.js';
import { loadModelObject } from 'viewer/lod.js';
import { makeLabelSprite } from 'viewer/spacecraft.js';

const GHOST_M_TO_KM = 1e-3;

function ringSprite(color) {
  const canvas = document.createElement('canvas');
  canvas.width = 64; canvas.height = 64;
  const ctx = canvas.getContext('2d');
  ctx.strokeStyle = color;
  ctx.lineWidth = 6;
  ctx.beginPath();
  ctx.arc(32, 32, 22, 0, 2 * Math.PI);
  ctx.stroke();
  const texture = new THREE.CanvasTexture(canvas);
  texture.colorSpace = THREE.SRGBColorSpace;
  const sprite = new THREE.Sprite(new THREE.SpriteMaterial({ map: texture, transparent: true, depthTest: false, opacity: 0.95 }));
  return sprite;
}

// Translucent copy of every material so the ghost reads as a ghost.
function ghostMaterials(object, color, opacity) {
  const tint = new THREE.Color(color);
  object.traverse((child) => {
    if (!child.isMesh) return;
    const mats = Array.isArray(child.material) ? child.material : [child.material];
    const out = mats.map((m) => {
      const g = m && m.clone ? m.clone() : new THREE.MeshStandardMaterial();
      g.transparent = true;
      g.opacity = opacity;
      g.depthWrite = false;
      if (g.color) g.color.lerp(tint, 0.35);
      if ('emissive' in g && g.emissive) { g.emissive.copy(tint); g.emissiveIntensity = 0.12; }
      return g;
    });
    child.material = Array.isArray(child.material) ? out : out[0];
    child.renderOrder = 2;
  });
}

// Wireframe boxes from the target spacecraft's links, for runs without a model.
function ghostBoxes(spec, color, opacity) {
  const group = new THREE.Group();
  for (const link of spec.links) {
    const geometry = new THREE.BoxGeometry(Math.max(link.dims_m[0], 1e-3), Math.max(link.dims_m[1], 1e-3), Math.max(link.dims_m[2], 1e-3));
    const mesh = new THREE.Mesh(geometry, new THREE.MeshStandardMaterial({ color, transparent: true, opacity, depthWrite: false, roughness: 0.6 }));
    mesh.add(new THREE.LineSegments(new THREE.EdgesGeometry(geometry), new THREE.LineBasicMaterial({ color, transparent: true, opacity: Math.min(1, opacity + 0.3) })));
    mesh.position.set(link.r_m[0], link.r_m[1], link.r_m[2]);
    mesh.quaternion.set(link.q[0], link.q[1], link.q[2], link.q[3]);
    group.add(mesh);
  }
  return group;
}

class StateTable {
  constructor(spec) {
    this.t = decodeFloat64(spec.t_s);
    this.pos = decodeFloat64(spec.pos_km);
    this.vel = spec.vel_kms ? decodeFloat32(spec.vel_kms) : null;
    this.q = spec.q ? decodeFloat32(spec.q) : null;
    this.count = Math.min(spec.count || this.t.length, this.t.length, this.pos.length / 3 | 0);
  }

  locate(time) {
    const n = this.count;
    if (n < 2) return { i: 0, f: 0 };
    if (time <= this.t[0]) return { i: 0, f: 0 };
    if (time >= this.t[n - 1]) return { i: n - 2, f: 1 };
    let lo = 0, hi = n - 1;
    while (hi - lo > 1) { const mid = (lo + hi) >> 1; if (this.t[mid] <= time) lo = mid; else hi = mid; }
    const dt = this.t[hi] - this.t[lo];
    return { i: lo, f: dt > 0 ? (time - this.t[lo]) / dt : 0 };
  }

  // false when `time` lies outside the table (the ghost is then not drawn).
  covers(time) { return this.count > 0 && time >= this.t[0] - 1e-9 && time <= this.t[this.count - 1] + 1e-9; }

  positionAt(time, out) {
    const { i, f } = this.locate(time);
    const a = 3 * i, b = 3 * (i + 1);
    if (this.count < 2) { out[0] = this.pos[a]; out[1] = this.pos[a + 1]; out[2] = this.pos[a + 2]; return out; }
    for (let k = 0; k < 3; k++) out[k] = this.pos[a + k] + f * (this.pos[b + k] - this.pos[a + k]);
    return out;
  }

  velocityAt(time, out) {
    const { i, f } = this.locate(time);
    const a = 3 * i, b = 3 * (i + 1);
    if (this.vel) {
      if (this.count < 2) { out[0] = this.vel[a]; out[1] = this.vel[a + 1]; out[2] = this.vel[a + 2]; return out; }
      for (let k = 0; k < 3; k++) out[k] = this.vel[a + k] + f * (this.vel[b + k] - this.vel[a + k]);
      return out;
    }
    if (this.count < 2) { out[0] = 0; out[1] = 0; out[2] = 0; return out; }
    const dt = this.t[i + 1] - this.t[i] || 1;
    for (let k = 0; k < 3; k++) out[k] = (this.pos[b + k] - this.pos[a + k]) / dt;
    return out;
  }

  quaternionAt(time, out) {
    if (!this.q) return null;
    const { i, f } = this.locate(time);
    if (this.count < 2) { for (let k = 0; k < 4; k++) out[k] = this.q[k]; return out; }
    return slerp(this.q, 4 * i, this.q, 4 * (i + 1), f, out);
  }
}

// specs: payload.references; sidecar/frames: the run; models: payload.models
// (keyed by spacecraft id), so a ghost can copy the target's model override.
export function createReferences(specs, sidecar, frames, models, options = {}) {
  const root = new THREE.Group();
  root.name = 'references';
  const items = [];
  const showPx = options.showPx ?? 3;
  const hidePx = options.hidePx ?? 0.6 * showPx;
  const opacity = options.opacity ?? 0.45;

  for (const spec of specs || []) {
    const table = new StateTable(spec);
    if (table.count < 1) continue;
    const target = Math.min(Math.max(0, (spec.target || 1) - 1), frames.sats - 1);
    const targetSpec = sidecar.spacecraft[target] || sidecar.spacecraft[0];
    const color = spec.color || '#ff8c69';
    const alpha = spec.opacity ?? opacity;

    // Full-history line: the whole reference trajectory, dimmed.
    const n = table.count;
    const lineGeometry = new THREE.BufferGeometry();
    lineGeometry.setAttribute('position', new THREE.BufferAttribute(new Float32Array(n * 3), 3));
    const line = new THREE.Line(lineGeometry, new THREE.LineBasicMaterial({ color: new THREE.Color(color), transparent: true, opacity: 0.55 }));
    line.frustumCulled = false;
    line.visible = spec.trail !== false;
    root.add(line);

    const marker = ringSprite(color);
    root.add(marker);
    const label = makeLabelSprite(spec.name || 'reference', color);
    root.add(label);

    // The ghost body: assembly-scaled group carrying the model copy or the boxes.
    const body = new THREE.Group();
    body.name = `reference:${spec.name}`;
    body.scale.setScalar(GHOST_M_TO_KM);
    body.visible = false;
    root.add(body);
    const model = models && models[String(targetSpec.id)];
    let status = null;
    if (model && model.url) {
      status = 'loading';
      loadModelObject(model, spec.name || 'reference', (object, info) => {
        ghostMaterials(object, color, alpha);
        body.add(object);
        status = info;
      }, (message) => {
        status = `failed: ${message}`;
        body.add(ghostBoxes(targetSpec, color, alpha));
      });
    } else {
      body.add(ghostBoxes(targetSpec, color, alpha));
      status = 'boxes';
    }

    items.push({
      spec, table, target, targetSpec, color, line, marker, label, body,
      radiusKm: Math.max(targetSpec.bounding_radius_m, 0.1) * GHOST_M_TO_KM,
      visible: false,
      get modelStatus() { return status; },
    });
  }

  const pos = new Float64Array(3), vel = new Float64Array(3), q = new Float32Array(4), tp = new Float64Array(3);
  const worldPos = new THREE.Vector3();
  const anchor = new Float64Array(3);
  let enabled = true;
  let labelsWanted = true;
  let lastAnchor = [NaN, NaN, NaN];

  function projectedPixels(item, camera, viewportHeight, groupMatrix) {
    worldPos.set(pos[0] - anchor[0], pos[1] - anchor[1], pos[2] - anchor[2]).applyMatrix4(groupMatrix);
    const dist = Math.max(1e-9, camera.position.distanceTo(worldPos));
    const halfHeight = Math.tan(THREE.MathUtils.degToRad(camera.fov) / 2);
    return (item.radiusKm / dist) * (viewportHeight / 2) / halfHeight;
  }

  function uploadLine(item) {
    const arr = item.line.geometry.getAttribute('position').array;
    const src = item.table.pos, n = item.table.count;
    for (let k = 0; k < n; k++) {
      arr[3 * k] = src[3 * k] - anchor[0];
      arr[3 * k + 1] = src[3 * k + 1] - anchor[1];
      arr[3 * k + 2] = src[3 * k + 2] - anchor[2];
    }
    item.line.geometry.getAttribute('position').needsUpdate = true;
    item.line.geometry.setDrawRange(0, n);
  }

  return {
    group: root,
    items,
    // groupMatrix: this group's matrixWorld; anchorKm: floating origin the group sits at
    update(t, camera, viewportHeight, groupMatrix, anchorKm) {
      anchor[0] = anchorKm[0]; anchor[1] = anchorKm[1]; anchor[2] = anchorKm[2];
      root.position.set(anchor[0], anchor[1], anchor[2]);
      const anchorMoved = anchor[0] !== lastAnchor[0] || anchor[1] !== lastAnchor[1] || anchor[2] !== lastAnchor[2];
      if (anchorMoved) lastAnchor = [anchor[0], anchor[1], anchor[2]];
      for (const item of items) {
        if (anchorMoved) uploadLine(item);
        const present = enabled && item.table.covers(t);
        item.marker.visible = present;
        item.label.visible = present && labelsWanted;
        if (!present) { item.body.visible = false; item.visible = false; continue; }
        item.table.positionAt(t, pos);
        const rel = [pos[0] - anchor[0], pos[1] - anchor[1], pos[2] - anchor[2]];
        item.marker.position.set(rel[0], rel[1], rel[2]);
        item.label.position.set(rel[0], rel[1], rel[2]);
        worldPos.set(rel[0], rel[1], rel[2]).applyMatrix4(groupMatrix);
        const dist = camera.position.distanceTo(worldPos);
        const h = 0.028 * dist;
        item.label.scale.set(h * item.label.userData.aspect, h, 1);
        const px = projectedPixels(item, camera, viewportHeight, groupMatrix);
        const show = item.visible ? px >= hidePx : px >= showPx;
        item.visible = show;
        item.body.visible = show;
        // The ring stays as a locator until the body is large enough to read.
        const ring = 0.018 * dist;
        item.marker.scale.set(ring, ring, 1);
        item.marker.material.opacity = show ? 0.35 : 0.95;
        if (!show) continue;
        item.body.position.set(rel[0], rel[1], rel[2]);
        if (item.table.quaternionAt(t, q) === null) {
          item.table.velocityAt(t, vel);
          velocityAlignedQuaternion(pos, vel, q);
        }
        item.body.quaternion.set(q[0], q[1], q[2], q[3]);
      }
    },
    // Distance (km) between reference k and the spacecraft it refers to, or NaN.
    separationKm(t, k) {
      const item = items[k];
      if (!item || !item.table.covers(t)) return NaN;
      frames.positionAt(t, item.target, tp);
      if (!Number.isFinite(tp[0])) return NaN;
      item.table.positionAt(t, pos);
      return Math.hypot(pos[0] - tp[0], pos[1] - tp[1], pos[2] - tp[2]);
    },
    setVisible(v) { enabled = v; root.visible = v; },
    setLabelsVisible(v) { labelsWanted = v; },
    modelStatus(k) { const it = items[k]; return it ? it.modelStatus : null; },
  };
}
