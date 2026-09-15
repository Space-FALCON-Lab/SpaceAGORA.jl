// Site terrain: nested square patches of the surface around a landing site,
// displaced by the digital elevation grids the page carries and draped with
// imagery that sharpens toward the center (from the global mosaic down to
// sub-meter mosaics). The patches sit in the globe group, so they rotate
// with the body; the globe itself is cut open under the outermost patch by
// the hole the globe material discards. All levels are drawn all the time:
// the finest level covers a few hundred meters, so the resolution rises by
// itself as the camera closes in.
//
// Payload contract (viewer_bundle.jl `terrain_payload`):
//   terrain.site               { lat_deg, lon_deg, height_m, name }
//   terrain.reference_radius_m
//   terrain.grids[]            { name, rows, cols, lat_min, lat_max, lon_min, lon_max, heights (base64 Float32, row-major north to south) }, finest first
//   terrain.imagery[]          { lat_min, lat_max, lon_min, lon_max, width, height, m_per_px, url (data: JPEG) }, coarse to fine
import * as THREE from 'three';
import { decodeFloat32 } from 'viewer/data.js';

const PATCH_SEGMENTS = 96;
const M_TO_KM = 1e-3;

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

export function createTerrain(spec, planet, options = {}) {
  const group = new THREE.Group();
  group.name = 'terrain';
  if (!spec || !spec.grids || spec.grids.length === 0) return { group, heightAt: () => NaN, levels: [], hole: null, site: null };
  const R = (spec.reference_radius_m || planet.equatorial_radius_m) * M_TO_KM;
  const grids = spec.grids.map((g) => ({ ...g, data: decodeFloat32(g.heights) }));
  const heightAt = (latDeg, lonDeg) => {
    for (const g of grids) { const h = gridHeight(g, latDeg, lonDeg); if (!Number.isNaN(h)) return h; }
    return NaN;
  };
  const bodyPoint = (latDeg, lonDeg, hM) => {
    const φ = THREE.MathUtils.degToRad(latDeg), λ = THREE.MathUtils.degToRad(lonDeg);
    const r = R + hM * M_TO_KM;
    return new THREE.Vector3(r * Math.cos(φ) * Math.cos(λ), r * Math.cos(φ) * Math.sin(λ), r * Math.sin(φ));
  };
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
    const positions = [], uvs = [], index = [];
    const heightOf = (lat, lon) => { const h = heightAt(lat, lon); return Number.isNaN(h) ? (spec.site ? spec.site.height_m : 0) : h; };
    const push = (lat, lon, h, u, vv) => { const p = bodyPoint(lat, lon, h); positions.push(p.x, p.y, p.z); uvs.push(u, vv); return positions.length / 3 - 1; };
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
    geometry.setIndex(index);
    geometry.computeVertexNormals();
    geometry.computeBoundingSphere();
    const texture = new THREE.TextureLoader().load(lvl.url);
    texture.colorSpace = THREE.SRGBColorSpace;
    texture.anisotropy = options.anisotropy || 4;
    texture.generateMipmaps = true;
    texture.minFilter = THREE.LinearMipmapLinearFilter;
    const material = new THREE.MeshLambertMaterial({ map: texture, side: THREE.DoubleSide });
    const mesh = new THREE.Mesh(geometry, material);
    mesh.name = `terrain-level-${k}`;
    mesh.renderOrder = 1 + k;
    group.add(mesh);
    levels.push({ mesh, spec: lvl, m_per_px: lvl.m_per_px });
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
    group.add(ring);
    siteMarker = ring;
  }
  return {
    group,
    levels,
    hole,
    site: spec.site || null,
    referenceRadiusKm: R,
    heightAt,
    setVisible(v) { group.visible = v; },
    get modelStatus() { return `${levels.length} levels, ${grids.length} grids, finest ${levels.length ? levels[levels.length - 1].m_per_px.toFixed(2) : '–'} m/px`; },
  };
}
