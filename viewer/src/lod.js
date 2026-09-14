// Close-up level of detail: one box assembly per spacecraft, built from the
// sidecar's link boxes, placed by the saved attitude (or a velocity-aligned
// one) and, when the run recorded them, by the per-step link poses.
//
// Units: the assembly group is scaled by 1e-3 so its children are authored
// in metres while the scene is in kilometres. Visibility is decided by how
// many pixels the spacecraft's bounding radius covers on screen, with
// hysteresis so the switch does not flicker.
import * as THREE from 'three';
import { STLLoader } from 'three/addons/loaders/STLLoader.js';
import { OBJLoader } from 'three/addons/loaders/OBJLoader.js';
import { GLTFLoader } from 'three/addons/loaders/GLTFLoader.js';
import { decodeBytes, velocityAlignedQuaternion } from 'viewer/data.js';

const M_TO_KM = 1e-3;
const ROOT_COLOR = 0xb9c4d2;
const LINK_COLOR = 0x6f8fb8;
const EDGE_COLOR = 0x1a2230;
const THRUSTER_COLOR = 0xff8a3d;
const FACET_COLOR = 0x5fd4ff;
const STL_COLOR = 0xc8cfd8;
const ARM_COLOR = 0xe8b04a;
const ARM_JOINT_COLOR = 0x3a4658;

// Robot arm: one cylinder per link along its own vector, placed each frame
// from the integrated arm state (inertial, relative to the spacecraft), so it
// lives in an unrotated sibling group rather than inside the body-frame assembly.
function buildArm(spec) {
  const group = new THREE.Group();
  group.name = 'arm';
  group.scale.setScalar(M_TO_KM);
  const links = [];
  for (const link of spec.links) {
    const v = new THREE.Vector3(link.vector_m[0], link.vector_m[1], link.vector_m[2]);
    const length = Math.max(v.length(), 1e-3);
    const holder = new THREE.Group(); // placed at the joint origin with the link quaternion
    const cyl = new THREE.Mesh(
      new THREE.CylinderGeometry(link.radius_m, link.radius_m, length, 16),
      new THREE.MeshStandardMaterial({ color: ARM_COLOR, roughness: 0.5, metalness: 0.3 }),
    );
    // CylinderGeometry runs along +y centred on the origin; move it to span 0..length along the link vector.
    cyl.quaternion.setFromUnitVectors(new THREE.Vector3(0, 1, 0), v.clone().normalize());
    cyl.position.copy(v).multiplyScalar(0.5);
    holder.add(cyl);
    const joint = new THREE.Mesh(new THREE.SphereGeometry(1.25 * link.radius_m, 12, 8), new THREE.MeshStandardMaterial({ color: ARM_JOINT_COLOR, roughness: 0.6 }));
    holder.add(joint);
    holder.userData.comOffset = new THREE.Vector3(link.com_offset_m[0], link.com_offset_m[1], link.com_offset_m[2]);
    group.add(holder);
    links.push(holder);
  }
  const tip = new THREE.Mesh(new THREE.SphereGeometry(0.8 * (spec.links[spec.links.length - 1]?.radius_m || 0.02), 12, 8), new THREE.MeshStandardMaterial({ color: 0xff6b6b, emissive: 0xff6b6b, emissiveIntensity: 0.4 }));
  group.add(tip);
  return { group, links, tip };
}

function boxMesh(dims, color) {
  const geometry = new THREE.BoxGeometry(Math.max(dims[0], 1e-3), Math.max(dims[1], 1e-3), Math.max(dims[2], 1e-3));
  const mesh = new THREE.Mesh(geometry, new THREE.MeshStandardMaterial({ color, roughness: 0.6, metalness: 0.15 }));
  const edges = new THREE.LineSegments(new THREE.EdgesGeometry(geometry), new THREE.LineBasicMaterial({ color: EDGE_COLOR }));
  mesh.add(edges);
  return mesh;
}

function orientToDirection(object, dir) {
  // Cones and planes are authored along +y / +z; rotate that axis onto `dir`.
  const d = new THREE.Vector3(dir[0], dir[1], dir[2]);
  if (d.lengthSq() === 0) return;
  d.normalize();
  object.quaternion.setFromUnitVectors(object.userData.axis, d);
}

function buildAssembly(spec, models, scLength) {
  const group = new THREE.Group();
  group.name = `assembly:${spec.name}`;
  group.scale.setScalar(M_TO_KM);
  const linkGroups = [];
  const glyphs = { thrusters: [], facets: [] };
  const r = Math.max(spec.bounding_radius_m, 0.1);
  spec.links.forEach((link, k) => {
    const g = new THREE.Group();
    g.position.set(link.r_m[0], link.r_m[1], link.r_m[2]);
    g.quaternion.set(link.q[0], link.q[1], link.q[2], link.q[3]);
    g.userData.boxes = boxMesh(link.dims_m, k === 0 ? ROOT_COLOR : LINK_COLOR);
    g.add(g.userData.boxes);
    group.add(g);
    linkGroups.push(g);
  });

  // 3D model override (STL, OBJ, glTF/GLB): the mesh replaces the boxes, glyphs stay on their links.
  const model = models && models[String(spec.id)];
  group.userData.modelStatus = model && model.url ? 'loading' : null;
  if (model && model.url) {
    const install = (object) => {
      let meshes = 0;
      object.traverse((child) => { if (child.isMesh) meshes++; });
      group.userData.modelStatus = `${model.format} (${meshes} mesh${meshes === 1 ? '' : 'es'}, ×${model.scale || 1})`;
      const s = model.scale || 1;
      object.scale.setScalar(s);
      const r = model.rotation_deg || [0, 0, 0];
      object.rotation.set(THREE.MathUtils.degToRad(r[0]), THREE.MathUtils.degToRad(r[1]), THREE.MathUtils.degToRad(r[2]), 'XYZ');
      // body = R * S * (v - c): translate by -R*S*c so the bounding-box centre sits on the spacecraft.
      const c = model.center || [0, 0, 0];
      object.position.set(-c[0] * s, -c[1] * s, -c[2] * s).applyEuler(object.rotation);
      object.traverse((child) => {
        if (child.isMesh) {
          if (!child.material || model.format === 'obj' || model.format === 'stl') {
            child.material = new THREE.MeshStandardMaterial({ color: STL_COLOR, roughness: 0.55, metalness: 0.2 });
          }
          child.castShadow = false;
        }
      });
      linkGroups[0].add(object);
      linkGroups.forEach((g) => { g.userData.boxes.visible = false; });
      group.userData.model = object;
    };
    try {
      const bytes = decodeBytes(model.url.split(',')[1] || '');
      const format = model.format || 'stl';
      if (format === 'stl') {
        const geometry = new STLLoader().parse(bytes.buffer);
        geometry.computeVertexNormals();
        install(new THREE.Mesh(geometry));
      } else if (format === 'obj') {
        install(new OBJLoader().parse(new TextDecoder().decode(bytes)));
      } else if (format === 'glb' || format === 'gltf') {
        const loader = new GLTFLoader();
        const payload = format === 'glb' ? bytes.buffer : new TextDecoder().decode(bytes);
        loader.parse(payload, '', (gltf) => install(gltf.scene), (err) => {
          group.userData.modelStatus = `failed: ${err && err.message ? err.message : err}`;
          console.warn(`glTF for ${spec.name} failed to parse; showing boxes instead.`, err);
        });
      } else {
        group.userData.modelStatus = `failed: unknown format ${format}`;
      }
    } catch (err) {
      group.userData.modelStatus = `failed: ${err && err.message ? err.message : err}`;
      console.warn(`Model for ${spec.name} could not be parsed; showing boxes instead.`, err);
    }
  }

  const coneH = Math.max(0.05, 0.12 * r), coneR = 0.35 * coneH;
  for (const th of spec.thrusters) {
    const g = linkGroups[th.link - 1];
    if (!g) continue;
    const cone = new THREE.Mesh(
      new THREE.ConeGeometry(coneR, coneH, 12),
      new THREE.MeshStandardMaterial({ color: THRUSTER_COLOR, emissive: THRUSTER_COLOR, emissiveIntensity: 0.35 }),
    );
    cone.userData.axis = new THREE.Vector3(0, 1, 0);
    // ConeGeometry points along +y with its base at -h/2; put the apex on the location.
    cone.geometry.translate(0, -coneH / 2, 0);
    cone.position.set(th.location_m[0], th.location_m[1], th.location_m[2]);
    orientToDirection(cone, th.direction);
    g.add(cone);
    glyphs.thrusters.push(cone);
  }
  for (const f of spec.facets) {
    const g = linkGroups[f.link - 1];
    if (!g) continue;
    const side = Math.sqrt(Math.max(f.area_m2, 1e-6));
    const plane = new THREE.Mesh(
      new THREE.PlaneGeometry(side, side),
      new THREE.MeshBasicMaterial({ color: FACET_COLOR, transparent: true, opacity: 0.28, side: THREE.DoubleSide, depthWrite: false }),
    );
    plane.userData.axis = new THREE.Vector3(0, 0, 1);
    plane.position.set(f.cp_m[0], f.cp_m[1], f.cp_m[2]);
    orientToDirection(plane, f.normal);
    const normalLine = new THREE.ArrowHelper(new THREE.Vector3(0, 0, 1), new THREE.Vector3(0, 0, 0), 0.6 * side, FACET_COLOR, 0.15 * side, 0.08 * side);
    plane.add(normalLine);
    g.add(plane);
    glyphs.facets.push(plane);
  }
  // Body axes so attitude reads at a glance: x red, y green, z blue.
  const axes = new THREE.AxesHelper(1.4 * r);
  axes.material.transparent = true;
  axes.material.opacity = 0.7;
  linkGroups[0].add(axes);
  group.userData.axes = axes;
  return { group, linkGroups, glyphs };
}

export function createAssemblies(sidecar, frames, options = {}) {
  const S = frames.sats;
  const specs = sidecar.spacecraft;
  const root = new THREE.Group();
  root.name = 'assemblies';
  const items = [];
  const showPx = options.showPx ?? 3;
  const hidePx = options.hidePx ?? 0.6 * showPx;
  const hideMarkerPx = options.hideMarkerPx ?? 12;

  // Assemblies are built for at most `assemblyLimit` spacecraft (a large
  // ensemble or constellation keeps markers only beyond that).
  const assemblyLimit = options.assemblyLimit ?? 256;
  for (let s = 0; s < S; s++) {
    const spec = specs[s] || specs[0];
    const radiusKm = Math.max(spec.bounding_radius_m, 0.1) * M_TO_KM;
    if (s >= assemblyLimit) { items.push({ spec, group: null, linkGroups: [], glyphs: { thrusters: [], facets: [] }, radiusKm }); continue; }
    const built = buildAssembly(spec, options.models, S);
    built.group.visible = false;
    root.add(built.group);
    let arm = null;
    if (spec.arm && spec.arm.links && spec.arm.links.length > 0 && frames.armLinkCount(s) === spec.arm.links.length) {
      arm = buildArm(spec.arm);
      arm.group.visible = false;
      root.add(arm.group);
    }
    items.push({ spec, ...built, arm, radiusKm });
  }

  const visible = new Uint8Array(S);      // assembly drawn
  const markerHidden = new Uint8Array(S); // marker suppressed because the assembly is large on screen
  const pxSize = new Float32Array(S);
  const pos = new Float64Array(3), vel = new Float64Array(3), q = new Float32Array(4);
  const lpBuf = new Float32Array(7 * 64);
  const armBuf = new Float32Array(7 * 64);
  const armQ = new THREE.Quaternion();
  const armCom = new THREE.Vector3();
  const worldPos = new THREE.Vector3();
  let enabled = options.enabled ?? true;
  let thrustersVisible = true, facetsVisible = true, axesVisible = true;

  const anchor = new Float64Array(3);

  function projectedPixels(item, camera, viewportHeight, groupMatrix) {
    worldPos.set(pos[0] - anchor[0], pos[1] - anchor[1], pos[2] - anchor[2]).applyMatrix4(groupMatrix);
    const dist = Math.max(1e-9, camera.position.distanceTo(worldPos));
    const halfHeight = Math.tan(THREE.MathUtils.degToRad(camera.fov) / 2);
    return (item.radiusKm / dist) * (viewportHeight / 2) / halfHeight;
  }

  return {
    group: root,
    items,
    visible,
    markerHidden,
    pxSize,
    // groupMatrix: this group's matrixWorld; anchorKm: Float64 [x, y, z] the group sits at
    update(t, camera, viewportHeight, groupMatrix, anchorKm) {
      anchor[0] = anchorKm[0]; anchor[1] = anchorKm[1]; anchor[2] = anchorKm[2];
      root.position.set(anchor[0], anchor[1], anchor[2]);
      for (let s = 0; s < S; s++) {
        const item = items[s];
        frames.positionAt(t, s, pos);
        const present = Number.isFinite(pos[0]) && Number.isFinite(pos[1]) && Number.isFinite(pos[2]);
        const px = enabled && present && item.group ? projectedPixels(item, camera, viewportHeight, groupMatrix) : 0;
        pxSize[s] = px;
        const wasVisible = visible[s] === 1;
        const show = wasVisible ? px >= hidePx : px >= showPx;
        visible[s] = show ? 1 : 0;
        markerHidden[s] = px >= hideMarkerPx ? 1 : 0;
        if (item.group) item.group.visible = show;
        if (item.arm) item.arm.group.visible = show;
        if (!show) continue;
        item.group.position.set(pos[0] - anchor[0], pos[1] - anchor[1], pos[2] - anchor[2]);
        if (frames.quaternionAt(t, s, q) === null) {
          frames.velocityAt(t, s, vel);
          velocityAlignedQuaternion(pos, vel, q);
        }
        item.group.quaternion.set(q[0], q[1], q[2], q[3]);
        if (item.arm) {
          // Same origin as the assembly, no body rotation: arm poses are inertial.
          item.arm.group.position.copy(item.group.position);
          const ap = frames.armPosesAt(t, s, armBuf);
          if (ap) {
            let tipPos = null;
            for (let k = 0; k < item.arm.links.length; k++) {
              const h = item.arm.links[k], o = 7 * k;
              armQ.set(ap[o + 3], ap[o + 4], ap[o + 5], ap[o + 6]);
              // COM = joint origin + R * com_offset  ->  joint origin = COM - R * com_offset
              armCom.set(ap[o], ap[o + 1], ap[o + 2]);
              h.position.copy(armCom).sub(h.userData.comOffset.clone().applyQuaternion(armQ));
              h.quaternion.copy(armQ);
              const spec = item.spec.arm.links[k];
              tipPos = h.position.clone().add(new THREE.Vector3(spec.vector_m[0], spec.vector_m[1], spec.vector_m[2]).applyQuaternion(armQ));
            }
            if (tipPos) item.arm.tip.position.copy(tipPos);
          }
        }
        const n = frames.linkCount(s);
        if (n > 0 && n + 1 === item.linkGroups.length) {
          const lp = frames.linkPosesAt(t, s, lpBuf);
          if (lp) {
            for (let k = 0; k < n; k++) {
              const g = item.linkGroups[k + 1], o = 7 * k;
              g.position.set(lp[o], lp[o + 1], lp[o + 2]);
              g.quaternion.set(lp[o + 3], lp[o + 4], lp[o + 5], lp[o + 6]);
            }
          }
        }
      }
    },
    setEnabled(v) { enabled = v; if (!v) { visible.fill(0); markerHidden.fill(0); for (const it of items) { if (it.group) it.group.visible = false; if (it.arm) it.arm.group.visible = false; } } },
    setThrustersVisible(v) { thrustersVisible = v; for (const it of items) for (const c of it.glyphs.thrusters) c.visible = v; },
    setFacetsVisible(v) { facetsVisible = v; for (const it of items) for (const f of it.glyphs.facets) f.visible = v; },
    setAxesVisible(v) { axesVisible = v; for (const it of items) if (it.group) it.group.userData.axes.visible = v; },
    get enabled() { return enabled; },
    modelStatus(s) { const it = items[s]; return it && it.group ? it.group.userData.modelStatus : null; },
  };
}
