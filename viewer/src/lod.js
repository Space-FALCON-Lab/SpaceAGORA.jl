// Close-up level of detail: one box assembly per spacecraft, built from the
// sidecar's link boxes, placed by the saved attitude (or a velocity-aligned
// one) and, when the run recorded them, by the per-step link poses.
//
// Units: the assembly group is scaled by 1e-3 so its children are authored
// in meters while the scene is in kilometers. Visibility is decided by how
// many pixels the spacecraft's bounding radius covers on screen, with
// hysteresis so the switch does not flicker.
import * as THREE from 'three';
import { STLLoader } from 'three/addons/loaders/STLLoader.js';
import { OBJLoader } from 'three/addons/loaders/OBJLoader.js';
import { GLTFLoader } from 'three/addons/loaders/GLTFLoader.js';
import { decodeBytes, velocityAlignedQuaternion, RotationTable } from 'viewer/data.js';
import { INFERNO_GLSL } from 'viewer/colormaps.js';

const M_TO_KM = 1e-3;
const ROOT_COLOR = 0xb9c4d2;
const LINK_COLOR = 0x6f8fb8;
const EDGE_COLOR = 0x1a2230;
const THRUSTER_COLOR = 0xff8a3d;
const FACET_COLOR = 0x5fd4ff;
const STL_COLOR = 0xc8cfd8;
const ARM_COLOR = 0xe8b04a;
const ARM_JOINT_COLOR = 0x3a4658;

// Heating overlay: every heatable face is shaded by the free-molecular
// incident energy flux 0.5 rho V^3 cos(theta) (W/m^2, full accommodation),
// theta the angle between the outward normal and the airspeed, through
// inferno on a log scale. Faces turned away from the flow stay cold.
const HEAT_VERTEX = `
  #include <common>
  varying vec3 vWorldNormal;
  #include <logdepthbuf_pars_vertex>
  void main() {
    vWorldNormal = normalize(mat3(modelMatrix) * normal);
    vec4 mvPosition = modelViewMatrix * vec4(position, 1.0);
    gl_Position = projectionMatrix * mvPosition;
    #include <logdepthbuf_vertex>
  }`;
const HEAT_FRAGMENT = `
  #include <common>
  uniform vec3 flowDir;
  uniform float q0;
  uniform float logLo;
  uniform float logHi;
  varying vec3 vWorldNormal;
  #include <logdepthbuf_pars_fragment>
  ${INFERNO_GLSL}
  void main() {
    #include <logdepthbuf_fragment>
    float c = max(dot(normalize(vWorldNormal), flowDir), 0.0);
    float q = q0 * c;
    float t = q > 0.0 ? clamp((log(q) / 2.302585 - logLo) / max(logHi - logLo, 1e-6), 0.0, 1.0) : 0.0;
    vec3 col = q > 0.0 ? inferno(t) : vec3(0.10, 0.11, 0.14);
    gl_FragColor = vec4(col, 1.0);
  }`;

function makeHeatMaterial() {
  return new THREE.ShaderMaterial({
    uniforms: { flowDir: { value: new THREE.Vector3(1, 0, 0) }, q0: { value: 0 }, logLo: { value: 0 }, logHi: { value: 1 } },
    vertexShader: HEAT_VERTEX, fragmentShader: HEAT_FRAGMENT, side: THREE.DoubleSide,
  });
}

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
    // CylinderGeometry runs along +y centerd on the origin; move it to span 0..length along the link vector.
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

// Parse a model override (an entry of payload.models: data URL, format,
// scale, rotation_deg, center) and hand back an object with the scale,
// rotation and centeing applied, authored in meters. `onReady(object,
// status)` runs synchronously for STL/OBJ and after the parse for glTF;
// `onFail(message)` when the bytes cannot be decoded. Shared by the
// assemblies and the reference ghosts so both draw the same geometry.
// Pose parts of a parsed model: every vertex inside an articulation's region
// (an axis-aligned box in model units, null bounds unbounded) is rotated by
// angle_deg about the axis through the pivot, in the model's own frame, and
// its normal with it. The bundler applies the same rotation to the aero mesh.
function articulateObject(object, articulations) {
  if (!articulations || articulations.length === 0) return;
  object.updateMatrixWorld(true);
  const p = new THREE.Vector3(), n = new THREE.Vector3(), pivot = new THREE.Vector3(), axis = new THREE.Vector3();
  const q = new THREE.Quaternion(), toLocal = new THREE.Matrix4(), normalToWorld = new THREE.Matrix3(), normalToLocal = new THREE.Matrix3();
  for (const art of articulations) {
    const lo = art.region.min.map((v) => (v == null ? -Infinity : v)), hi = art.region.max.map((v) => (v == null ? Infinity : v));
    axis.set(art.axis[0], art.axis[1], art.axis[2]).normalize();
    q.setFromAxisAngle(axis, THREE.MathUtils.degToRad(art.angle_deg));
    pivot.set(art.pivot[0], art.pivot[1], art.pivot[2]);
    object.traverse((child) => {
      if (!child.isMesh) return;
      const geometry = child.geometry;
      const pos = geometry.getAttribute('position');
      const nrm = geometry.getAttribute('normal');
      const toWorld = child.matrixWorld;
      toLocal.copy(toWorld).invert();
      normalToWorld.getNormalMatrix(toWorld);
      normalToLocal.copy(normalToWorld).invert();
      let touched = false;
      for (let i = 0; i < pos.count; i++) {
        p.fromBufferAttribute(pos, i).applyMatrix4(toWorld);
        if (p.x < lo[0] || p.x > hi[0] || p.y < lo[1] || p.y > hi[1] || p.z < lo[2] || p.z > hi[2]) continue;
        p.sub(pivot).applyQuaternion(q).add(pivot).applyMatrix4(toLocal);
        pos.setXYZ(i, p.x, p.y, p.z);
        if (nrm) {
          n.fromBufferAttribute(nrm, i).applyMatrix3(normalToWorld).applyQuaternion(q).applyMatrix3(normalToLocal).normalize();
          nrm.setXYZ(i, n.x, n.y, n.z);
        }
        touched = true;
      }
      if (touched) { pos.needsUpdate = true; if (nrm) nrm.needsUpdate = true; geometry.computeBoundingSphere(); geometry.computeBoundingBox(); }
    });
  }
}

export function loadModelObject(model, label, onReady, onFail) {
  const install = (object) => {
    articulateObject(object, model.articulations);
    let meshes = 0;
    object.traverse((child) => { if (child.isMesh) meshes++; });
    const s = model.scale || 1;
    object.scale.setScalar(s);
    const r = model.rotation_deg || [0, 0, 0];
    object.rotation.set(THREE.MathUtils.degToRad(r[0]), THREE.MathUtils.degToRad(r[1]), THREE.MathUtils.degToRad(r[2]), 'XYZ');
    // body = R * S * (v - c): translate by -R*S*c so the bounding-box center sits on the spacecraft.
    const c = model.center || [0, 0, 0];
    object.position.set(-c[0] * s, -c[1] * s, -c[2] * s).applyEuler(object.rotation);
    object.traverse((child) => {
      if (child.isMesh) {
        if (!child.material || model.format === 'obj' || model.format === 'stl') {
          child.material = new THREE.MeshStandardMaterial({ color: STL_COLOR, roughness: 0.55, metalness: 0.2 });
        }
        child.userData.heatable = true;
        child.castShadow = true;
        child.receiveShadow = true;
      }
    });
    onReady(object, `${model.format} (${meshes} mesh${meshes === 1 ? '' : 'es'}, ×${s})`);
  };
  const message = (err) => (err && err.message ? err.message : String(err));
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
      loader.parse(payload, '', (gltf) => {
        try { install(gltf.scene); } catch (err) { onFail(message(err)); }
      }, (err) => onFail(message(err)));
    } else {
      onFail(`unknown format ${format}`);
    }
  } catch (err) {
    onFail(message(err));
  }
}

function boxMesh(dims, color) {
  const geometry = new THREE.BoxGeometry(Math.max(dims[0], 1e-3), Math.max(dims[1], 1e-3), Math.max(dims[2], 1e-3));
  const mesh = new THREE.Mesh(geometry, new THREE.MeshStandardMaterial({ color, roughness: 0.6, metalness: 0.15 }));
  mesh.userData.heatable = true;
  mesh.castShadow = true;
  mesh.receiveShadow = true;
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
    loadModelObject(model, spec.name, (object, status) => {
      linkGroups[0].add(object);
      linkGroups.forEach((g) => { g.userData.boxes.visible = false; });
      group.userData.model = object;
      group.userData.modelStatus = status;
    }, (message) => {
      group.userData.modelStatus = `failed: ${message}`;
      console.warn(`Model for ${spec.name} could not be parsed; showing boxes instead.`, message);
    });
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
    cone.userData.uiOnly = true;
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
    plane.userData.uiOnly = true;
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
    items.push({ spec, ...built, arm, radiusKm, heatMaterial: makeHeatMaterial(), heated: false });
  }

  // Heating overlay: needs density and velocity in the frames and the planet's spin for the airspeed.
  const planet = sidecar.planet;
  const heatingAvailable = frames.hasScalar('density') && !!frames.vel;
  const rotation = planet && planet.rotation ? new RotationTable(planet.rotation, planet.spin_rad_s || [0, 0, 0]) : null;
  const spinRate = planet && planet.spin_rad_s ? planet.spin_rad_s[2] : 0;
  let heating = false;
  // Range of 0.5 rho V^3 over the run (W/m^2), three decades below the peak.
  let heatLogHi = 1, heatLogLo = 0;
  if (heatingAvailable) {
    let peak = 0;
    for (let k = 0; k < frames.count; k++) {
      for (let s = 0; s < S; s++) {
        const rho = frames.scalarAtIndex('density', k, s), v = frames.scalarAtIndex('speed', k, s) * 1000;
        if (Number.isFinite(rho) && rho > 0 && Number.isFinite(v)) peak = Math.max(peak, 0.5 * rho * v * v * v);
      }
    }
    if (peak > 0) { heatLogHi = Math.log10(peak); heatLogLo = heatLogHi - 3; }
  }
  const qPlanet = new Float32Array(4), omega = new THREE.Vector3(), rVec = new THREE.Vector3(), vRel = new THREE.Vector3(), qTmp = new THREE.Quaternion();
  const posA = new Float64Array(3), velA = new Float64Array(3);
  // Airspeed (km/s, inertial axes): the inertial velocity minus the co-rotating atmosphere.
  function airspeedAt(t, s, out) {
    frames.positionAt(t, s, posA);
    frames.velocityAt(t, s, velA);
    vRel.set(velA[0], velA[1], velA[2]);
    if (rotation && spinRate !== 0) {
      rotation.at(t, qPlanet);
      qTmp.set(qPlanet[0], qPlanet[1], qPlanet[2], qPlanet[3]);
      omega.set(0, 0, spinRate).applyQuaternion(qTmp);
      rVec.set(posA[0], posA[1], posA[2]);
      vRel.sub(omega.cross(rVec));
    }
    out[0] = vRel.x; out[1] = vRel.y; out[2] = vRel.z;
    return out;
  }
  // Body attitude [x, y, z, w]: the saved quaternion, else velocity-aligned (what the assembly is drawn with).
  const velB = new Float64Array(3), posB = new Float64Array(3);
  function attitudeAt(t, s, out) {
    if (frames.quaternionAt(t, s, out) !== null) return out;
    frames.positionAt(t, s, posB);
    frames.velocityAt(t, s, velB);
    return velocityAlignedQuaternion(posB, velB, out);
  }
  // Pose of link k (0-based) in the body frame: position (m) and quaternion, from the recorded
  // link poses when the run saved them, else the configured placement.
  const lpPose = new Float32Array(7 * 64);
  function linkPoseAt(t, s, k, outP, outQ) {
    const item = items[s];
    const link = item.spec.links[k];
    if (k > 0 && frames.linkCount(s) === item.spec.links.length - 1) {
      const lp = frames.linkPosesAt(t, s, lpPose);
      if (lp) {
        const o = 7 * (k - 1);
        outP[0] = lp[o]; outP[1] = lp[o + 1]; outP[2] = lp[o + 2];
        outQ[0] = lp[o + 3]; outQ[1] = lp[o + 4]; outQ[2] = lp[o + 5]; outQ[3] = lp[o + 6];
        return true;
      }
    }
    outP[0] = link.r_m[0]; outP[1] = link.r_m[1]; outP[2] = link.r_m[2];
    outQ[0] = link.q[0]; outQ[1] = link.q[1]; outQ[2] = link.q[2]; outQ[3] = link.q[3];
    return false;
  }
  function applyHeat(item, on) {
    if (item.heated === on || !item.group) return;
    item.heated = on;
    item.group.traverse((child) => {
      if (!child.isMesh || !child.userData.heatable) return;
      if (on) { child.userData.baseMaterial = child.material; child.material = item.heatMaterial; }
      else if (child.userData.baseMaterial) { child.material = child.userData.baseMaterial; }
    });
  }
  const airTmp = new Float64Array(3);
  function updateHeat(item, s, t, groupMatrix) {
    airspeedAt(t, s, airTmp);
    vRel.set(airTmp[0], airTmp[1], airTmp[2]);
    const speed = vRel.length() * 1000;
    const rho = frames.scalarAt('density', t, s);
    const u = item.heatMaterial.uniforms;
    u.q0.value = Number.isFinite(rho) && rho > 0 ? 0.5 * rho * speed * speed * speed : 0;
    u.logLo.value = heatLogLo; u.logHi.value = heatLogHi;
    if (speed > 0) u.flowDir.value.copy(vRel).normalize().transformDirection(groupMatrix);
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

  // Face picking: the first heatable mesh (a link box or a model mesh) under
  // the ray among the drawn assemblies, with the face normal and the hit
  // point in that link's frame (meters), so the pick follows the link's pose.
  const relMatrix = new THREE.Matrix4(), normalMatrix = new THREE.Matrix3(), nTmp = new THREE.Vector3(), pTmp = new THREE.Vector3();
  function ancestorsVisible(o) { for (let a = o; a; a = a.parent) if (!a.visible) return false; return true; }
  function pickFace(raycaster) {
    const targets = [];
    for (let s = 0; s < S; s++) if (visible[s] && items[s].group) targets.push(items[s].group);
    if (targets.length === 0) return null;
    const hits = raycaster.intersectObjects(targets, true);
    for (const h of hits) {
      const m = h.object;
      if (!m.isMesh || !m.userData.heatable || !h.face || !ancestorsVisible(m)) continue;
      // Walk up to the assembly group; the child of it on the way is the link group.
      let linkGroup = m, assembly = m.parent;
      while (assembly && assembly.parent !== root) { linkGroup = assembly; assembly = assembly.parent; }
      if (!assembly) continue;
      const s = items.findIndex((it) => it.group === assembly);
      if (s < 0) continue;
      const k = Math.max(0, items[s].linkGroups.indexOf(linkGroup));
      const lg = items[s].linkGroups[k];
      relMatrix.copy(lg.matrixWorld).invert().multiply(m.matrixWorld);
      normalMatrix.getNormalMatrix(relMatrix);
      nTmp.copy(h.face.normal).applyMatrix3(normalMatrix).normalize();
      pTmp.copy(h.point).applyMatrix4(relMatrix.copy(lg.matrixWorld).invert());
      return {
        s, link: k, source: m === lg.userData.boxes ? 'box' : 'model',
        normal: [nTmp.x, nTmp.y, nTmp.z], point: [pTmp.x, pTmp.y, pTmp.z],
      };
    }
    return null;
  }
  let faceMarker = null;
  function clearFaceMarker() {
    if (faceMarker) { faceMarker.parent && faceMarker.parent.remove(faceMarker); faceMarker = null; }
  }
  function setFaceMarker(pick) {
    clearFaceMarker();
    if (!pick) return;
    const item = items[pick.s];
    const lg = item && item.linkGroups[pick.link];
    if (!lg) return;
    const r = Math.max(item.spec.bounding_radius_m, 0.1);
    const n = new THREE.Vector3(pick.normal[0], pick.normal[1], pick.normal[2]);
    const p = new THREE.Vector3(pick.point[0], pick.point[1], pick.point[2]);
    const marker = new THREE.Group();
    marker.name = 'face-marker';
    marker.userData.uiOnly = true;
    const arrow = new THREE.ArrowHelper(n, p, 0.45 * r, 0xf2b950, 0.12 * r, 0.06 * r);
    marker.add(arrow);
    const dot = new THREE.Mesh(new THREE.SphereGeometry(0.025 * r, 12, 8), new THREE.MeshBasicMaterial({ color: 0xf2b950, depthTest: false }));
    dot.position.copy(p);
    dot.renderOrder = 10;
    marker.add(dot);
    lg.add(marker);
    faceMarker = marker;
  }

  return {
    group: root,
    items,
    visible,
    markerHidden,
    pxSize,
    airspeedAt,
    attitudeAt,
    linkPoseAt,
    pickFace,
    setFaceMarker,
    clearFaceMarker,
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
        if (heating && heatingAvailable) { applyHeat(item, true); updateHeat(item, s, t, groupMatrix); } else applyHeat(item, false);
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
    // Heating overlay: on/off, and its legend range in W/m^2 (log scale, three decades below the run's peak).
    heatingAvailable,
    setHeatingVisible(v) { heating = v && heatingAvailable; if (!heating) for (const it of items) applyHeat(it, false); },
    heatRange() { return { lo: Math.pow(10, heatLogLo), hi: Math.pow(10, heatLogHi), log: true }; },
    get enabled() { return enabled; },
    modelStatus(s) { const it = items[s]; return it && it.group ? it.group.userData.modelStatus : null; },
  };
}
