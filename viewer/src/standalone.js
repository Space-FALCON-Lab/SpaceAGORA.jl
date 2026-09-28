// Standalone data path: build the viewer payload in the browser from raw
// trajectory data (a CSV or JSON table, SpaceAGORA's own results CSV
// included), a planet name and an epoch, with no simulation and no Julia.
// The page viewer/standalone.html drives this from a form; the build script
// viewer/build_standalone.py can also embed a data file so the page opens
// straight into the run.
//
// What it reproduces from the Julia bundler: the scene sidecar (planet spec
// with an IAU rotation table, one box link per spacecraft), the frame blocks
// (as typed arrays; data.js accepts those directly), model overrides
// (centerd from the parsed geometry), reference ghosts and paths. Planet
// rotation uses the IAU 2009 pole and prime-meridian polynomials, evaluated
// at the epoch plus elapsed time, so a ground track lands within a fraction
// of a degree of the SPICE-driven one.
import * as THREE from 'three';
import { STLLoader } from 'three/addons/loaders/STLLoader.js';
import { OBJLoader } from 'three/addons/loaders/OBJLoader.js';
import { GLTFLoader } from 'three/addons/loaders/GLTFLoader.js';
import { decodeBytes } from 'viewer/data.js';

// IAU 2009 (Archinal et al. 2011) pole (RA, Dec at J2000, degrees) and prime
// meridian W = w0 + w1 * d (degrees, d = days past J2000 TDB); radii in meters.
// The Moon's periodic pole terms and Earth's nutation are omitted (< 0.1 deg).
export const PLANETS = {
  earth: { name: 'Earth', re: 6378137.0, rp: 6356752.3, ra: 0.0, dec: 90.0, w0: 190.147, w1: 360.9856235 },
  mars: { name: 'Mars', re: 3396190.0, rp: 3376200.0, ra: 317.68143, dec: 52.88650, w0: 176.630, w1: 350.89198226 },
  venus: { name: 'Venus', re: 6051800.0, rp: 6051800.0, ra: 272.76, dec: 67.16, w0: 160.20, w1: -1.4813688 },
  moon: { name: 'Moon', re: 1737400.0, rp: 1737400.0, ra: 269.9949, dec: 66.5392, w0: 38.3213, w1: 13.17635815 },
  titan: { name: 'Titan', re: 2575000.0, rp: 2575000.0, ra: 39.4827, dec: 83.4279, w0: 186.5855, w1: 22.5769768 },
};

const DEG = Math.PI / 180;
const J2000_MS = Date.UTC(2000, 0, 1, 11, 58, 55, 816); // 2000-01-01T12:00:00 TDB is 11:58:55.816 UTC

// Seconds past J2000 for an ISO UTC string (leap seconds ignored: < 1 min, i.e. < 0.3 deg of Earth rotation).
export function secondsPastJ2000(utcIso) {
  const ms = Date.parse(utcIso);
  if (Number.isNaN(ms)) throw new Error(`epoch "${utcIso}" is not an ISO 8601 UTC date`);
  return (ms - J2000_MS) / 1000;
}

// Active planet-fixed -> J2000 quaternion (scalar-last) at seconds past J2000,
// the same object the Julia bundler stores as q_pi: rot(q) is the passive
// J2000 -> body matrix, and the viewer applies q to the globe mesh.
export function iauRotationQuaternion(planet, tPastJ2000) {
  const d = tPastJ2000 / 86400;
  const ra = planet.ra * DEG, dec = planet.dec * DEG, W = (planet.w0 + planet.w1 * d) * DEG;
  const z = [Math.cos(dec) * Math.cos(ra), Math.cos(dec) * Math.sin(ra), Math.sin(dec)];
  const Q = [-Math.sin(ra), Math.cos(ra), 0]; // node of the planet equator on the J2000 equator
  const zxQ = [z[1] * Q[2] - z[2] * Q[1], z[2] * Q[0] - z[0] * Q[2], z[0] * Q[1] - z[1] * Q[0]];
  const x = [Q[0] * Math.cos(W) + zxQ[0] * Math.sin(W), Q[1] * Math.cos(W) + zxQ[1] * Math.sin(W), Q[2] * Math.cos(W) + zxQ[2] * Math.sin(W)];
  const y = [z[1] * x[2] - z[2] * x[1], z[2] * x[0] - z[0] * x[2], z[0] * x[1] - z[1] * x[0]];
  // columns x, y, z of the active body -> inertial matrix
  const m = new THREE.Matrix4().makeBasis(new THREE.Vector3(...x), new THREE.Vector3(...y), new THREE.Vector3(...z));
  const q = new THREE.Quaternion().setFromRotationMatrix(m);
  return [q.x, q.y, q.z, q.w];
}

export function rotationTable(planet, tStartPastJ2000, times, maxSamples = 4096) {
  const n = times.length;
  const stride = Math.max(1, Math.ceil(n / maxSamples));
  const t_s = [], q_pi = [];
  let prev = null;
  for (let k = 0; k < n; k += stride) {
    const q = iauRotationQuaternion(planet, tStartPastJ2000 + times[k]);
    if (prev && prev[0] * q[0] + prev[1] * q[1] + prev[2] * q[2] + prev[3] * q[3] < 0) for (let c = 0; c < 4; c++) q[c] = -q[c];
    t_s.push(times[k]); q_pi.push(q); prev = q;
  }
  if (t_s[t_s.length - 1] !== times[n - 1]) {
    const q = iauRotationQuaternion(planet, tStartPastJ2000 + times[n - 1]);
    if (prev && prev[0] * q[0] + prev[1] * q[1] + prev[2] * q[2] + prev[3] * q[3] < 0) for (let c = 0; c < 4; c++) q[c] = -q[c];
    t_s.push(times[n - 1]); q_pi.push(q);
  }
  return { t_s, q_pi };
}

// ---------------------------------------------------------------------------
// Tables
// ---------------------------------------------------------------------------

// CSV text -> { columns: [names], rows: number of rows, data: {name: Float64Array} }.
export function parseCsv(text) {
  const lines = text.split(/\r?\n/).filter((l) => l.trim().length > 0 && !l.startsWith('#'));
  if (lines.length < 2) throw new Error('the CSV needs a header line and at least one data row');
  const sep = lines[0].includes(',') ? ',' : /\t/.test(lines[0]) ? '\t' : /;/.test(lines[0]) ? ';' : /\s+/;
  const columns = lines[0].split(sep).map((c) => c.trim().replace(/^"|"$/g, ''));
  const rows = lines.length - 1;
  const data = {};
  for (const c of columns) data[c] = new Float64Array(rows);
  for (let r = 0; r < rows; r++) {
    const parts = lines[r + 1].split(sep);
    for (let c = 0; c < columns.length; c++) {
      const v = parts[c] === undefined ? NaN : Number(parts[c]);
      data[columns[c]][r] = v;
    }
  }
  return { columns, rows, data };
}

// A spacecraft's column set from a parsed table: SpaceAGORA's `sc{i}_...`
// layout, or the generic `t,x,y,z[,vx,vy,vz][,qx,qy,qz,qw]` (optionally with
// a `sc{i}_` prefix). Returns null when the prefix has no position columns.
function columnSet(table, prefix) {
  const has = (n) => Object.prototype.hasOwnProperty.call(table.data, n);
  const pick = (...names) => names.find(has) || null;
  const pos = [pick(`${prefix}pos_1`, `${prefix}x`, `${prefix}pos_x`, `${prefix}X`), pick(`${prefix}pos_2`, `${prefix}y`, `${prefix}pos_y`, `${prefix}Y`), pick(`${prefix}pos_3`, `${prefix}z`, `${prefix}pos_z`, `${prefix}Z`)];
  if (pos.some((c) => c === null)) return null;
  const vel = [pick(`${prefix}vel_1`, `${prefix}vx`, `${prefix}vel_x`), pick(`${prefix}vel_2`, `${prefix}vy`, `${prefix}vel_y`), pick(`${prefix}vel_3`, `${prefix}vz`, `${prefix}vel_z`)];
  const q = [pick(`${prefix}q_1`, `${prefix}qx`), pick(`${prefix}q_2`, `${prefix}qy`), pick(`${prefix}q_3`, `${prefix}qz`), pick(`${prefix}q_4`, `${prefix}qw`)];
  const linkPose = [];
  for (let k = 1; has(`${prefix}link_pose_${k}`); k++) linkPose.push(`${prefix}link_pose_${k}`);
  return {
    pos,
    vel: vel.every((c) => c !== null) ? vel : null,
    q: q.every((c) => c !== null) ? q : null,
    mass: pick(`${prefix}mass`, `${prefix}mass_kg`),
    density: pick(`${prefix}density`, `${prefix}density_kg_m3`),
    heat: pick(`${prefix}heat_rate`),
    drag: [pick(`${prefix}drag_1`), pick(`${prefix}drag_2`), pick(`${prefix}drag_3`)].every((c) => c !== null) ? [`${prefix}drag_1`, `${prefix}drag_2`, `${prefix}drag_3`] : null,
    wind: [pick(`${prefix}wind_1`), pick(`${prefix}wind_2`), pick(`${prefix}wind_3`)].every((c) => c !== null) ? [`${prefix}wind_1`, `${prefix}wind_2`, `${prefix}wind_3`] : null,
    linkPose: linkPose.length % 7 === 0 && linkPose.length > 0 ? linkPose : [],
  };
}

// Discover the spacecraft in a table: every `sc{i}_` prefix with positions, else the unprefixed generic set.
export function discoverSpacecraft(table) {
  const found = [];
  const seen = new Set();
  for (const c of table.columns) {
    const m = /^(sc\d+_)/.exec(c);
    if (m && !seen.has(m[1])) { seen.add(m[1]); const set = columnSet(table, m[1]); if (set) found.push({ prefix: m[1], set }); }
  }
  if (found.length === 0) {
    const set = columnSet(table, '');
    if (set) found.push({ prefix: '', set });
  }
  found.sort((a, b) => (parseInt(a.prefix.slice(2), 10) || 0) - (parseInt(b.prefix.slice(2), 10) || 0));
  return found;
}

function timeColumn(table) {
  const names = ['time', 't', 'time_s', 't_s', 'elapsed_s', 'seconds'];
  const c = names.find((n) => Object.prototype.hasOwnProperty.call(table.data, n));
  if (!c) throw new Error(`no time column (one of ${names.join(', ')})`);
  return table.data[c];
}

// Build the frame blocks for one table. `lengthUnit` 'm' or 'km'; `maxFrames` decimates evenly.
export function buildFrames(table, spacecraft, { lengthUnit = 'm', maxFrames = 4000 } = {}) {
  const tAll = timeColumn(table);
  const n = table.rows, S = spacecraft.length;
  const stride = Math.max(1, Math.ceil(n / maxFrames));
  const idx = [];
  for (let r = 0; r < n; r += stride) idx.push(r);
  if (idx[idx.length - 1] !== n - 1) idx.push(n - 1);
  const N = idx.length;
  const toKm = lengthUnit === 'km' ? 1 : 1e-3;
  const toM = lengthUnit === 'km' ? 1e3 : 1;
  const hasVel = spacecraft.every((s) => s.set.vel), hasQ = spacecraft.every((s) => s.set.q), hasMass = spacecraft.every((s) => s.set.mass);
  const hasDensity = spacecraft.every((s) => s.set.density), hasHeat = spacecraft.every((s) => s.set.heat);
  const hasDrag = spacecraft.every((s) => s.set.drag), hasWind = spacecraft.every((s) => s.set.wind);
  const counts = spacecraft.map((s) => s.set.linkPose.length / 7);
  const hasLp = counts.some((c) => c > 0);
  const offsets = []; let total = 0;
  for (const c of counts) { offsets.push(total); total += 7 * c; }
  const t_s = new Float64Array(N), pos = new Float64Array(N * S * 3);
  const vel = hasVel ? new Float32Array(N * S * 3) : null, q = hasQ ? new Float32Array(N * S * 4) : null;
  const mass = hasMass ? new Float32Array(N * S) : null, density = hasDensity ? new Float32Array(N * S) : null, heat = hasHeat ? new Float32Array(N * S) : null;
  const drag = hasDrag ? new Float32Array(N * S) : null, wind = hasWind ? new Float32Array(N * S * 3) : null;
  const lp = hasLp ? new Float32Array(N * total) : null;
  const d = table.data;
  for (let f = 0; f < N; f++) {
    const r = idx[f];
    t_s[f] = tAll[r];
    for (let s = 0; s < S; s++) {
      const set = spacecraft[s].set, b3 = (f * S + s) * 3;
      for (let c = 0; c < 3; c++) pos[b3 + c] = d[set.pos[c]][r] * toKm;
      if (vel) for (let c = 0; c < 3; c++) vel[b3 + c] = d[set.vel[c]][r] * toKm;
      if (q) for (let c = 0; c < 4; c++) q[(f * S + s) * 4 + c] = d[set.q[c]][r];
      if (mass) mass[f * S + s] = d[set.mass][r];
      if (density) density[f * S + s] = d[set.density][r];
      if (heat) heat[f * S + s] = d[set.heat][r];
      if (drag) drag[f * S + s] = Math.hypot(d[set.drag[0]][r], d[set.drag[1]][r], d[set.drag[2]][r]);
      if (wind) for (let c = 0; c < 3; c++) wind[b3 + c] = d[set.wind[c]][r];
      if (lp && counts[s] > 0) for (let k = 0; k < 7 * counts[s]; k++) lp[f * total + offsets[s] + k] = d[set.linkPose[k]][r] * (k % 7 < 3 ? toM : 1);
    }
  }
  return {
    count: N, sats: S, source_rows: n, stride_rows: stride,
    t_dtype: 'f64', t_s, pos_dtype: 'f64', pos_km: pos, vel_kms: vel, q, mass_kg: mass, density_kg_m3: density,
    heat_rate_w_m2: heat, drag_n: drag, wind_ms: wind,
    link_pose: hasLp ? { stride: 7, counts, offsets, total, data: lp } : null, arm_pose: null,
  };
}

// Reference ghost table from a generic/SpaceAGORA CSV (first spacecraft in it).
export function buildReference(table, { name = 'reference', lengthUnit = 'm', target = 1, color = '#ff8c69', opacity = 0.45 } = {}) {
  const sc = discoverSpacecraft(table);
  if (sc.length === 0) throw new Error(`reference "${name}": no position columns`);
  const set = sc[0].set, t = timeColumn(table), n = table.rows, d = table.data;
  const toKm = lengthUnit === 'km' ? 1 : 1e-3;
  const pos = new Float64Array(3 * n), vel = set.vel ? new Float32Array(3 * n) : null, q = set.q ? new Float32Array(4 * n) : null;
  const t_s = new Float64Array(n);
  for (let r = 0; r < n; r++) {
    t_s[r] = t[r];
    for (let c = 0; c < 3; c++) { pos[3 * r + c] = d[set.pos[c]][r] * toKm; if (vel) vel[3 * r + c] = d[set.vel[c]][r] * toKm; }
    if (q) for (let c = 0; c < 4; c++) q[4 * r + c] = d[set.q[c]][r];
  }
  return { name, target, count: n, t_s, pos_km: pos, vel_kms: vel, q, color, opacity, trail: true };
}

// ---------------------------------------------------------------------------
// Scene
// ---------------------------------------------------------------------------

export function buildScene({ planetKey, epochUtc, frames, spacecraftSpecs }) {
  const planet = PLANETS[planetKey];
  if (!planet) throw new Error(`unknown planet "${planetKey}" (one of ${Object.keys(PLANETS).join(', ')})`);
  const t0 = secondsPastJ2000(epochUtc);
  const rotation = rotationTable(planet, t0, Array.from(frames.t_s));
  const spin = planet.w1 * DEG / 86400;
  const spacecraft = spacecraftSpecs.map((spec, i) => {
    const dims = spec.dims_m || [1, 1, 1];
    const links = [{ name: 'root', root: true, dims_m: dims, r_m: [0, 0, 0], q: [0, 0, 0, 1], mass_kg: spec.mass_kg || 100 }];
    for (const l of spec.links || []) links.push({ name: l.name || `link${links.length}`, root: false, dims_m: l.dims_m || [0.5, 0.5, 0.1], r_m: l.r_m || [0, 0, 0], q: l.q || [0, 0, 0, 1], mass_kg: l.mass_kg || 10 });
    const radius = Math.max(...links.map((l) => 0.5 * Math.hypot(...l.dims_m) + Math.hypot(...l.r_m)));
    return { id: spec.id ?? i + 1, name: spec.name || `sc${i + 1}`, links, thrusters: [], facets: [], joints: [], bounding_radius_m: radius, stl_path: null, arm: null };
  });
  return {
    schema: 1,
    epoch: { utc: new Date(Date.parse(epochUtc)).toISOString(), et_start_s: t0 },
    planet: {
      name: planet.name, equatorial_radius_m: planet.re, polar_radius_m: planet.rp, spin_rad_s: [0, 0, spin],
      inertial_frame: 'J2000', texture: planetKey, rotation,
    },
    spacecraft,
    orientation_sim: !!frames.q,
    atmosphere: null,
    results: { feather: '', link_pose: { field: 'link_pose', stride: 7, layout: ['rx', 'ry', 'rz', 'qx', 'qy', 'qz', 'qw'], columns: '' }, arm_pose: null },
  };
}

// ---------------------------------------------------------------------------
// Models
// ---------------------------------------------------------------------------

export function modelFormat(filename) {
  const ext = (filename.split('.').pop() || '').toLowerCase();
  if (ext === 'stl') return { format: 'stl', mime: 'model/stl' };
  if (ext === 'obj') return { format: 'obj', mime: 'model/obj' };
  if (ext === 'glb') return { format: 'glb', mime: 'model/gltf-binary' };
  if (ext === 'gltf') return { format: 'gltf', mime: 'model/gltf+json' };
  throw new Error(`unsupported model format .${ext} (use .stl, .obj, .glb or .gltf)`);
}

// Bounding-box center of a model (model units), parsed with the same loaders the viewer uses.
export function modelCenter(dataUrl, format) {
  return new Promise((resolve, reject) => {
    const bytes = decodeBytes(dataUrl.split(',')[1] || '');
    const finish = (object) => {
      const box = new THREE.Box3().setFromObject(object);
      const c = new THREE.Vector3(); box.getCenter(c);
      resolve([c.x, c.y, c.z]);
    };
    try {
      if (format === 'stl') finish(new THREE.Mesh(new STLLoader().parse(bytes.buffer)));
      else if (format === 'obj') finish(new OBJLoader().parse(new TextDecoder().decode(bytes)));
      else new GLTFLoader().parse(format === 'glb' ? bytes.buffer : new TextDecoder().decode(bytes), '', (g) => finish(g.scene), reject);
    } catch (err) { reject(err); }
  });
}

// { url, format, scale, rotation_deg, center, source } for payload.models[id].
export async function buildModel({ dataUrl, filename, scale = 1, rotation_deg = [0, 0, 0], center = true, articulations = [] }) {
  const { format } = modelFormat(filename);
  const c = center ? await modelCenter(dataUrl, format) : [0, 0, 0];
  return { url: dataUrl, format, scale, rotation_deg, center: c, articulations, source: filename };
}

export function fileToDataUrl(file, mime) {
  return new Promise((resolve, reject) => {
    const reader = new FileReader();
    reader.onload = () => resolve(mime ? reader.result.replace(/^data:[^;]*;/, `data:${mime};`) : reader.result);
    reader.onerror = () => reject(reader.error);
    reader.readAsDataURL(file);
  });
}

// ---------------------------------------------------------------------------
// Assembly
// ---------------------------------------------------------------------------

/**
 * Build the full viewer payload.
 *   spec = { planet: 'venus', epoch: '1993-05-26T00:00:07Z', csv: text | table: {columns, rows, data} | json: {...},
 *            lengthUnit: 'm'|'km', maxFrames, spacecraft: [{ id, name, dims_m, mass_kg, links }],
 *            models: [{ id, dataUrl, filename, scale, rotation_deg }], references: [{ name, csv, lengthUnit, target, color }],
 *            paths: [...], options: { title, trail_orbits, frame, speed } }
 *   textures = window.SPACEAGORA_TEXTURES (planetKey -> { url, lon_left_deg, resolution, ... }) or {}
 */
export async function buildPayload(spec, textures = {}) {
  let table = spec.table || null;
  if (!table && spec.csv) table = parseCsv(spec.csv);
  if (!table && spec.json) table = tableFromJson(spec.json);
  if (!table) throw new Error('no trajectory data (csv, table or json)');
  const found = discoverSpacecraft(table);
  if (found.length === 0) throw new Error('no position columns found (expected sc1_pos_1..3 or x,y,z)');
  const specs = found.map((f, i) => Object.assign({ id: i + 1, name: f.prefix ? f.prefix.replace(/_$/, '') : 'sc1' }, (spec.spacecraft || [])[i] || {}));
  const frames = buildFrames(table, found, { lengthUnit: spec.lengthUnit || 'm', maxFrames: spec.maxFrames || 4000 });
  const scene = buildScene({ planetKey: (spec.planet || 'earth').toLowerCase(), epochUtc: spec.epoch || '2000-01-01T12:00:00Z', frames, spacecraftSpecs: specs });
  const models = {};
  for (const m of spec.models || []) {
    const id = m.id ?? specs[0].id;
    models[String(id)] = await buildModel(m);
  }
  const references = (spec.references || []).map((r) => buildReference(r.table || parseCsv(r.csv), r));
  const key = (spec.planet || 'earth').toLowerCase();
  const tex = textures[key] ? { [key]: textures[key] } : {};
  return { scene, frames, textures: tex, models, paths: spec.paths || [], references, options: Object.assign({ frame: 'inertial' }, spec.options || {}) };
}

// JSON form: { time: [...], spacecraft: [{ name, pos: [[x,y,z],...], vel: [[...]], q: [[x,y,z,w],...] }] }
// or a column dictionary { time: [...], sc1_pos_1: [...], ... }.
export function tableFromJson(doc) {
  if (Array.isArray(doc.spacecraft)) {
    const time = Float64Array.from(doc.time || doc.t);
    const columns = ['time'], data = { time };
    doc.spacecraft.forEach((sc, i) => {
      const p = `sc${i + 1}_`;
      for (let c = 0; c < 3; c++) { columns.push(`${p}pos_${c + 1}`); data[`${p}pos_${c + 1}`] = Float64Array.from(sc.pos.map((v) => v[c])); }
      if (sc.vel) for (let c = 0; c < 3; c++) { columns.push(`${p}vel_${c + 1}`); data[`${p}vel_${c + 1}`] = Float64Array.from(sc.vel.map((v) => v[c])); }
      if (sc.q) for (let c = 0; c < 4; c++) { columns.push(`${p}q_${c + 1}`); data[`${p}q_${c + 1}`] = Float64Array.from(sc.q.map((v) => v[c])); }
    });
    return { columns, rows: time.length, data };
  }
  const columns = Object.keys(doc);
  const data = {};
  for (const c of columns) data[c] = Float64Array.from(doc[c]);
  return { columns, rows: data[columns[0]].length, data };
}

// ---------------------------------------------------------------------------
// Ensembles from several tables
// ---------------------------------------------------------------------------

// Linear interpolation of a column onto a time grid; NaN outside the sample's span.
function resample(tSrc, ySrc, tDst) {
  const out = new Float64Array(tDst.length);
  let j = 0;
  for (let i = 0; i < tDst.length; i++) {
    const t = tDst[i];
    if (t < tSrc[0] || t > tSrc[tSrc.length - 1]) { out[i] = NaN; continue; }
    while (j < tSrc.length - 2 && tSrc[j + 1] < t) j++;
    const dt = tSrc[j + 1] - tSrc[j];
    const f = dt > 0 ? (t - tSrc[j]) / dt : 0;
    out[i] = ySrc[j] + f * (ySrc[j + 1] - ySrc[j]);
  }
  return out;
}

/**
 * Several trajectory files as one ensemble: each file's first spacecraft becomes
 * a pseudo-spacecraft on the first file's time grid; the first file is the
 * nominal unless a sample has `nominal: true`. `spec.samples = [{ label, csv | json | table, scalar }]`.
 */
export async function buildEnsemblePayload(spec, textures = {}) {
  const tables = spec.samples.map((s) => s.table || (s.csv ? parseCsv(s.csv) : tableFromJson(s.json)));
  const sets = tables.map((t, i) => { const f = discoverSpacecraft(t); if (!f.length) throw new Error(`sample ${i + 1}: no position columns`); return f[0].set; });
  const tAxis = timeColumn(tables[0]);
  const n = spec.samples.length;
  const merged = { columns: ['time'], rows: tAxis.length, data: { time: Float64Array.from(tAxis) } };
  for (let i = 0; i < n; i++) {
    const t = timeColumn(tables[i]), d = tables[i].data, p = `sc${i + 1}_`;
    for (let c = 0; c < 3; c++) { merged.columns.push(`${p}pos_${c + 1}`); merged.data[`${p}pos_${c + 1}`] = resample(t, d[sets[i].pos[c]], tAxis); }
    if (sets.every((s) => s.vel)) for (let c = 0; c < 3; c++) { merged.columns.push(`${p}vel_${c + 1}`); merged.data[`${p}vel_${c + 1}`] = resample(t, d[sets[i].vel[c]], tAxis); }
    if (sets.every((s) => s.q)) for (let c = 0; c < 4; c++) { merged.columns.push(`${p}q_${c + 1}`); merged.data[`${p}q_${c + 1}`] = resample(t, d[sets[i].q[c]], tAxis); }
  }
  const nominal = Math.max(0, spec.samples.findIndex((s) => s.nominal));
  const single = Object.assign({}, spec, { table: merged, csv: null, json: null, samples: null,
    spacecraft: spec.samples.map((s, i) => Object.assign({ id: i + 1, name: s.label || `sample ${i + 1}` }, (spec.spacecraft || [])[0] || {})) });
  const payload = await buildPayload(single, textures);
  // one model applies to every sample when given
  if (spec.models && spec.models.length) { const m = payload.models[String(spec.models[0].id ?? 1)]; if (m) for (let i = 0; i < n; i++) payload.models[String(i + 1)] = m; }
  payload.ensemble = {
    count: n, spacecraft_per_sample: 1, scalar_name: spec.scalarName || '',
    samples: spec.samples.map((s, i) => ({ index: i + 1, seed: '', success: true, scalar: Number.isFinite(s.scalar) ? s.scalar : null, label: s.label || `sample ${i + 1}`, directory: '' })),
    nominal,
  };
  return payload;
}
