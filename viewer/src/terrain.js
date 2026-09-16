// Site terrain: a view-dependent quadtree over the region the payload's
// imagery covers, displaced by the digital elevation grids the page carries.
// The tree sits in the globe group, so it rotates with the body; the globe
// itself is cut open under the covered region by the hole the globe material
// discards.
//
// Why a quadtree. A lander sees the ground from 230 km away at powered descent
// initiation and from 4 m away at touchdown. The horizon shrinks faster than
// the resolution demand rises (at 15 km the horizon is 231 km and 11 m/px
// fills the frame; at 300 m it is 32 km and 0.22 m/px), and a tile of side s
// at distance d covers about 1349 s/d pixels of a 1400 px view, so one tree
// walked per frame from the root, splitting a node whose projected size
// exceeds `splitPixels`, spends the triangles and the texels where the camera
// is actually looking. Five fixed squares around the site could not: they left
// the rest of the horizon on the 2.7 km/px global mosaic with a hard square
// seam where they ended.
//
// Texture inheritance. Nodes exist in the payload only where imagery was
// built, so coverage may be a funnel (wide and coarse over the descent
// corridor, narrow and fine at the site) rather than a full pyramid. A node
// with no tile of its own draws the nearest present ancestor's tile through
// the sub-rectangle of its UV range that belongs to it. Every point inside the
// covered region is therefore always textured at the best resolution that
// exists for it, and the transition between levels is a gradual loss of
// sharpness with distance rather than an edge.
//
// Payload contract (viewer_bundle.jl `terrain_payload`):
//   terrain.site               { lat_deg, lon_deg, height_m, name }
//   terrain.reference_radius_m
//   terrain.grids[]            { name, rows, cols, lat_min, lat_max, lon_min, lon_max, heights (base64 Float32, row-major north to south) }, finest first
//   terrain.tiles              { scheme: "quadtree", root: { lat_min, lat_max, lon_min, lon_max }, tile_px, max_level,
//                                nodes: [ { level, x, y, url (data: JPEG), m_per_px } ],
//                                and optionally resolution { finest_m_per_px, feature_scale_m, ... },
//                                attribution [] and source, which the panel reports }
//     A node at (level, x, y) covers lon_min + (lon_max - lon_min) * x / 2^level
//     eastward by one node width, and latitude from lat_max downward the same
//     way, so y = 0 is the northern row (the grids' convention).
//   terrain.imagery[]          legacy nested squares; only the widest is used, as the root tile.
import * as THREE from 'three';
import { decodeFloat32 } from 'viewer/data.js';

const TERRAIN_M_TO_KM = 1e-3;
const TERRAIN_NODE_SEGMENTS = 16;      // quads per node edge; a node's DEM sampling
const TERRAIN_SPLIT_PIXELS = 256;      // split a node projecting wider than this (one tile's texels)
const TERRAIN_BUILD_BUDGET = 6;        // node builds per frame, so a split never stalls the frame
const TERRAIN_EXTRA_LEVELS = 3;        // levels below the deepest tile, for relief the DEM still has
const TERRAIN_NODE_CAP = 1200;         // built nodes kept; the least recently drawn are dropped
const TERRAIN_TEXTURE_CAP = 512;       // decoded tiles kept on the GPU
const TERRAIN_FADE_FRAC = 0.07;        // outer band of the root that blends into the globe
const TERRAIN_GRID_TAPER = 0.08;       // band inside a finer grid over which it blends into the coarser one
const TERRAIN_SKIRT_RELIEF = 0.6;      // skirt depth as a fraction of a node's own relief
const TERRAIN_SKIRT_SIZE = 0.02;       // plus this fraction of the node's width
const TERRAIN_SKIRT_MAX_M = 3000;      // and never deeper than this
const TERRAIN_SKIRT_SINK = 1e-4;       // the skirt's top, below the node's own edge, as a fraction of its width

// Height at (lat, lon) from one grid, bilinear, NaN outside it.
function terrainGridHeight(g, latDeg, lonDeg) {
  let lon = lonDeg;
  while (lon < g.lon_min - 1e-12 && lon + 360 <= g.lon_max + 1e-9) lon += 360;
  while (lon > g.lon_max + 1e-12 && lon - 360 >= g.lon_min - 1e-9) lon -= 360;
  if (latDeg < g.lat_min || latDeg > g.lat_max || lon < g.lon_min || lon > g.lon_max) return NaN;
  return terrainGridSample(g, latDeg, lon);
}

// The same, with the sample point clamped into the grid: the coarsest grid
// extends its edge values outward instead of stepping down to the sphere.
function terrainGridSample(g, latDeg, lonDeg) {
  const rows = g.rows, cols = g.cols;
  const dlat = (g.lat_max - g.lat_min) / rows, dlon = (g.lon_max - g.lon_min) / cols;
  let fr = (g.lat_max - latDeg) / dlat - 0.5, fc = (lonDeg - g.lon_min) / dlon - 0.5;
  fr = Math.min(rows - 1, Math.max(0, fr)); fc = Math.min(cols - 1, Math.max(0, fc));
  const r0 = Math.min(rows - 2, Math.floor(fr)), c0 = Math.min(cols - 2, Math.floor(fc));
  const tr = fr - r0, tc = fc - c0;
  const h = g.data;
  const h00 = h[r0 * cols + c0], h01 = h[r0 * cols + c0 + 1], h10 = h[(r0 + 1) * cols + c0], h11 = h[(r0 + 1) * cols + c0 + 1];
  return (1 - tr) * ((1 - tc) * h00 + tc * h01) + tr * ((1 - tc) * h10 + tc * h11);
}

// 1 well inside the grid, falling to 0 at its edges over `TERRAIN_GRID_TAPER`
// of its span: a fine grid blends into the coarser surface under it instead of
// ending in a cliff.
function terrainGridWeight(g, latDeg, lonDeg) {
  const latSpan = g.lat_max - g.lat_min, lonSpan = g.lon_max - g.lon_min;
  const a = Math.min(latDeg - g.lat_min, g.lat_max - latDeg) / (latSpan * TERRAIN_GRID_TAPER);
  const b = Math.min(lonDeg - g.lon_min, g.lon_max - lonDeg) / (lonSpan * TERRAIN_GRID_TAPER);
  const t = Math.min(1, Math.max(0, Math.min(a, b)));
  return t * t * (3 - 2 * t);
}

// A quadtree over the legacy nested squares: only the widest one can be placed,
// as the root tile. Payloads written since the quadtree landed carry `tiles`.
function terrainTilesFromImagery(imagery) {
  if (!imagery || imagery.length === 0) return null;
  const widest = imagery.slice().sort((a, b) => (b.lat_max - b.lat_min) - (a.lat_max - a.lat_min))[0];
  return {
    scheme: 'quadtree',
    root: { lat_min: widest.lat_min, lat_max: widest.lat_max, lon_min: widest.lon_min, lon_max: widest.lon_max },
    tile_px: 256,
    max_level: 0,
    nodes: [{ level: 0, x: 0, y: 0, url: widest.url, m_per_px: widest.m_per_px }],
  };
}

/**
 * Build the terrain quadtree. `options`:
 *   anisotropy     max anisotropic filtering of the tile textures
 *   splitPixels    projected node width that triggers a split (default 256)
 *   buildBudget    node builds allowed per frame (default 6)
 *   segments       quads per node edge (default 16)
 * The returned object's `update(camera, viewportHeight)` walks the tree for the
 * current view; `setGlobeTexture(map, lonLeftDeg)` hands it the globe's own map
 * so the outermost ring of the covered region fades into it.
 */
export function createTerrain(spec, planet, options = {}) {
  const group = new THREE.Group();
  group.name = 'terrain';
  const emptyStats = { nodes: 0, drawn: 0, triangles: 0, pending: 0, built: 0, updateMs: 0, maxLevel: 0 };
  if (!spec || !spec.grids || spec.grids.length === 0) {
    return {
      group, levels: [], hole: null, site: null, heightAt: () => NaN, stats: emptyStats,
      update() {}, setGlobeTexture() {}, setVisible(v) { group.visible = v; }, modelStatus: 'no terrain',
    };
  }
  const tiles = spec.tiles && spec.tiles.nodes && spec.tiles.nodes.length ? spec.tiles : terrainTilesFromImagery(spec.imagery);
  if (!tiles) {
    return {
      group, levels: [], hole: null, site: null, heightAt: () => NaN, stats: emptyStats,
      update() {}, setGlobeTexture() {}, setVisible(v) { group.visible = v; }, modelStatus: 'no imagery',
    };
  }

  const R = (spec.reference_radius_m || planet.equatorial_radius_m) * TERRAIN_M_TO_KM;
  const DEG_TO_KM = R * Math.PI / 180;
  const splitPixels = options.splitPixels || TERRAIN_SPLIT_PIXELS;
  const buildBudget = options.buildBudget || TERRAIN_BUILD_BUDGET;
  const segments = options.segments || TERRAIN_NODE_SEGMENTS;
  const root = tiles.root;
  const rootSpan = root.lat_max - root.lat_min;
  const rootLonSpan = root.lon_max - root.lon_min;
  const maxTileLevel = Number.isFinite(tiles.max_level) ? tiles.max_level : 0;
  const maxLevel = maxTileLevel + TERRAIN_EXTRA_LEVELS;

  // Grids, coarsest first: each finer one is blended into the surface below it.
  const grids = spec.grids.map((g) => ({ ...g, data: decodeFloat32(g.heights) }))
    .sort((a, b) => (b.lat_max - b.lat_min) - (a.lat_max - a.lat_min));
  const finest = grids[grids.length - 1];

  // The published height: the finest grid that covers the point, NaN off them
  // all. This is the radar-altitude lookup the panel and the dust decals use,
  // so it stays strict.
  const heightAt = (latDeg, lonDeg) => {
    for (let k = grids.length - 1; k >= 0; k--) {
      const h = terrainGridHeight(grids[k], latDeg, lonDeg);
      if (!Number.isNaN(h)) return h;
    }
    return NaN;
  };

  // The surface the tree draws: the coarsest grid extended to its edge values,
  // every finer grid blended in over its own taper band, and the whole thing
  // faded to the sphere across the outer ring of the covered region, so the
  // terrain meets the globe with neither a step nor a seam.
  const fadeOf = (latDeg, lonDeg) => {
    const a = Math.min(latDeg - root.lat_min, root.lat_max - latDeg) / (rootSpan * TERRAIN_FADE_FRAC);
    const b = Math.min(lonDeg - root.lon_min, root.lon_max - lonDeg) / (rootLonSpan * TERRAIN_FADE_FRAC);
    const t = Math.min(1, Math.max(0, Math.min(a, b)));
    return 1 - t * t * (3 - 2 * t);
  };
  const surfaceHeight = (latDeg, lonDeg) => {
    let h = terrainGridSample(grids[0], latDeg, lonDeg);
    for (let k = 1; k < grids.length; k++) {
      const v = terrainGridHeight(grids[k], latDeg, lonDeg);
      if (Number.isNaN(v)) continue;
      h += terrainGridWeight(grids[k], latDeg, lonDeg) * (v - h);
    }
    return h;
  };

  const bodyPoint = (latDeg, lonDeg, hM, out) => {
    const lat = THREE.MathUtils.degToRad(latDeg), lon = THREE.MathUtils.degToRad(lonDeg);
    const r = R + hM * TERRAIN_M_TO_KM;
    return out.set(r * Math.cos(lat) * Math.cos(lon), r * Math.cos(lat) * Math.sin(lon), r * Math.sin(lat));
  };

  // ---- tiles -------------------------------------------------------------
  // Every tile the payload carries, by "level/x/y", decoded on first use and
  // dropped again when too many are resident.
  const tileKey = (level, x, y) => `${level}/${x}/${y}`;
  const tileIndex = new Map();
  for (const n of tiles.nodes) tileIndex.set(tileKey(n.level, n.x, n.y), { ...n, texture: null, ready: false, used: 0 });
  let tilesResident = 0;
  const anisotropy = options.anisotropy || 4;
  function tileTexture(entry, frameIndex) {
    entry.used = frameIndex;
    if (entry.texture) return entry.ready ? entry.texture : null;
    const texture = new THREE.TextureLoader().load(entry.url, () => { entry.ready = true; });
    texture.colorSpace = THREE.SRGBColorSpace;
    texture.anisotropy = anisotropy;
    texture.generateMipmaps = true;
    texture.minFilter = THREE.LinearMipmapLinearFilter;
    texture.wrapS = THREE.ClampToEdgeWrapping;
    texture.wrapT = THREE.ClampToEdgeWrapping;
    entry.texture = texture;
    tilesResident++;
    if (tilesResident > TERRAIN_TEXTURE_CAP) terrainEvictTiles(frameIndex);
    return null;
  }
  function terrainEvictTiles(frameIndex) {
    const resident = [];
    for (const entry of tileIndex.values()) if (entry.texture) resident.push(entry);
    resident.sort((a, b) => a.used - b.used);
    for (let k = 0; k < resident.length - TERRAIN_TEXTURE_CAP * 0.75; k++) {
      const entry = resident[k];
      if (entry.used >= frameIndex - 1) break;
      entry.texture.dispose();
      entry.texture = null; entry.ready = false;
      tilesResident--;
    }
  }

  // The nearest present ancestor (or the node itself) whose tile has decoded,
  // with the UV rectangle of that tile this node covers. Loading is started for
  // the deepest present ancestor, so a node sharpens as its own tile arrives.
  const tileSource = { key: '', texture: null, u0: 0, v0: 0, size: 1 };
  function resolveTile(level, x, y, frameIndex) {
    for (let d = 0; level - d >= 0; d++) {
      const lv = level - d, ax = x >> d, ay = y >> d;
      const entry = tileIndex.get(tileKey(lv, ax, ay));
      if (!entry) continue;
      const texture = tileTexture(entry, frameIndex);
      if (!texture) continue;                    // present but still decoding: keep climbing
      const scale = 1 / (1 << d);
      tileSource.key = tileKey(lv, ax, ay);
      tileSource.texture = texture;
      tileSource.u0 = (x - (ax << d)) * scale;
      tileSource.v0 = 1 - (y - (ay << d)) * scale;   // v = 1 is the tile's northern edge
      tileSource.size = scale;
      return tileSource;
    }
    return null;
  }

  // ---- materials ---------------------------------------------------------
  // Lambert, so the scene's lights and shadows apply as they did to the fixed
  // patches, plus the globe's own map blended in over the outer ring: at the
  // boundary of the covered region the terrain is exactly the globe, so the
  // region has no visible edge. One program serves every node.
  let globeMap = null;
  const globeLonLeft = Number.isFinite(options.globeLonLeft) ? options.globeLonLeft : -180;
  const materials = [];
  function terrainMaterial() {
    const material = new THREE.MeshLambertMaterial({ side: THREE.DoubleSide });
    material.onBeforeCompile = (shader) => {
      shader.uniforms.uGlobeMap = { value: globeMap };
      shader.uniforms.uGlobeMix = { value: globeMap ? 1 : 0 };
      shader.vertexShader = shader.vertexShader
        .replace('#include <common>', '#include <common>\nattribute vec2 aGlobeUv;\nattribute float aFade;\nvarying vec2 vGlobeUv;\nvarying float vFade;')
        .replace('#include <begin_vertex>', '#include <begin_vertex>\nvGlobeUv = aGlobeUv;\nvFade = aFade;');
      shader.fragmentShader = shader.fragmentShader
        .replace('#include <common>', '#include <common>\nuniform sampler2D uGlobeMap;\nuniform float uGlobeMix;\nvarying vec2 vGlobeUv;\nvarying float vFade;')
        // A double-sided material flips the normal on a back face. The ground is
        // only ever seen from above, and the skirt carries the ground's own
        // normal so that it disappears into the ground where it shows: flipping
        // it turned the skirt black, which is what drew a dotted line along
        // every node boundary. Undo the flip (faceDirection is +/-1).
        .replace('#include <normal_fragment_begin>', `#include <normal_fragment_begin>
  normal *= faceDirection;`)
        .replace('#include <map_fragment>', `#include <map_fragment>
  diffuseColor.rgb = mix(diffuseColor.rgb, texture2D(uGlobeMap, vGlobeUv).rgb, vFade * uGlobeMix);`);
      material.userData.shader = shader;
    };
    material.customProgramCacheKey = () => 'terrain-quadtree';
    materials.push(material);
    return material;
  }
  function setGlobeTexture(map) {
    globeMap = map;
    for (const material of materials) {
      if (!material.userData.shader) continue;
      material.userData.shader.uniforms.uGlobeMap.value = map;
      material.userData.shader.uniforms.uGlobeMix.value = map ? 1 : 0;
    }
  }

  // ---- nodes -------------------------------------------------------------
  const nodes = new Map();
  const tmpPoint = new THREE.Vector3();
  const tmpNormal = new THREE.Vector3();
  let builtCount = 0;

  function makeNode(level, x, y) {
    const size = rootSpan / (1 << level);
    const lonSize = rootLonSpan / (1 << level);
    const latMax = root.lat_max - size * y, latMin = latMax - size;
    const lonMin = root.lon_min + lonSize * x, lonMax = lonMin + lonSize;
    const node = {
      level, x, y, latMin, latMax, lonMin, lonMax,
      sizeKm: Math.max(size, lonSize) * DEG_TO_KM,
      center: bodyPoint(0.5 * (latMin + latMax), 0.5 * (lonMin + lonMax), surfaceHeight(0.5 * (latMin + latMax), 0.5 * (lonMin + lonMax)), new THREE.Vector3()),
      radiusKm: 0, mesh: null, children: null, parent: null, texKey: '', used: -1,
    };
    node.radiusKm = 0.75 * node.sizeKm;   // refined from the real heights when the node is built
    nodes.set(tileKey(level, x, y), node);
    return node;
  }

  function nodeAt(level, x, y) {
    return nodes.get(tileKey(level, x, y)) || makeNode(level, x, y);
  }

  // One node's geometry: an n x n patch of the surface sampled at the node's
  // own resolution, authored relative to the node center (so the GPU never
  // handles 1737 km numbers at meter detail), with a skirt hanging from all
  // four edges to hide the crack against a coarser neighbor.
  function buildNode(node, frameIndex) {
    const n = segments;
    const source = resolveTile(node.level, node.x, node.y, frameIndex);
    const positions = [], normals = [], uvs = [], globeUvs = [], fades = [], index = [];
    const center = node.center;
    const dLat = (node.latMax - node.latMin) / n, dLon = (node.lonMax - node.lonMin) / n;
    const latOf = (i) => node.latMax - dLat * i;
    const lonOf = (j) => node.lonMin + dLon * j;
    const u0 = source ? source.u0 : 0, v0 = source ? source.v0 : 1, us = source ? source.size : 1;
    // Heights on a one-vertex halo around the node, so a normal is the central
    // difference of the same surface on both sides of every edge. Normals
    // averaged from a node's own faces alone lean inward along its border, and
    // the neighbor's lean the other way, which draws a shading line along every
    // boundary of the tree; the halo makes the two agree by construction.
    const stride = n + 3;
    const halo = new Float64Array(stride * stride);
    for (let i = -1; i <= n + 1; i++) {
      const lat = latOf(i);
      for (let j = -1; j <= n + 1; j++) {
        const lon = lonOf(j);
        halo[(i + 1) * stride + (j + 1)] = surfaceHeight(lat, lon) * (1 - fadeOf(lat, lon));
      }
    }
    const heightOf = (i, j) => halo[(i + 1) * stride + (j + 1)];
    const latMeters = dLat * DEG_TO_KM * 1000;
    // The surface normal from the height gradient: up - dh/dNorth * north - dh/dEast * east.
    const normalAt = (i, j, lat, lon, out) => {
      const la = THREE.MathUtils.degToRad(lat), lo = THREE.MathUtils.degToRad(lon);
      const cla = Math.cos(la), sla = Math.sin(la), clo = Math.cos(lo), slo = Math.sin(lo);
      const lonMeters = Math.max(dLon * DEG_TO_KM * 1000 * cla, 1e-6);
      const dhdN = (heightOf(i - 1, j) - heightOf(i + 1, j)) / (2 * latMeters);
      const dhdE = (heightOf(i, j + 1) - heightOf(i, j - 1)) / (2 * lonMeters);
      return out.set(
        cla * clo + dhdN * sla * clo + dhdE * slo,
        cla * slo + dhdN * sla * slo - dhdE * clo,
        sla - dhdN * cla,
      ).normalize();
    };
    let hMin = Infinity, hMax = -Infinity;
    for (let i = 0; i <= n; i++) {
      const lat = latOf(i);
      for (let j = 0; j <= n; j++) {
        const lon = lonOf(j);
        const h = heightOf(i, j);
        if (h < hMin) hMin = h;
        if (h > hMax) hMax = h;
        bodyPoint(lat, lon, h, tmpPoint).sub(center);
        positions.push(tmpPoint.x, tmpPoint.y, tmpPoint.z);
        normalAt(i, j, lat, lon, tmpNormal);
        normals.push(tmpNormal.x, tmpNormal.y, tmpNormal.z);
        uvs.push(u0 + us * j / n, v0 - us * i / n);
        globeUvs.push((lon - globeLonLeft) / 360, (lat + 90) / 180);
        fades.push(fadeOf(lat, lon));
      }
    }
    for (let i = 0; i < n; i++) {
      for (let j = 0; j < n; j++) {
        const a = i * (n + 1) + j, b = a + 1, c = a + (n + 1), d = c + 1;
        index.push(a, c, b, b, c, d);
      }
    }
    // Skirts: deep enough to cover the height a coarser neighbor may sit at.
    //
    // Two neighboring nodes are separate meshes placed at their own centers, so
    // their shared edge lands a hair apart once the positions are rounded to
    // float32, and the pixels whose centers fall in that sub-pixel gap show
    // what is behind the ground: the skirt. That is harmless as long as the
    // skirt looks like the ground, and it did not. The skirt hung from the
    // surface ring, so `computeVertexNormals` had it cancel against itself (it
    // is indexed in both windings) and left the wall with no normal at all;
    // giving it the ground's normal is not enough either, because a
    // double-sided material flips the normal on a back face and the wall is
    // seen from both. So the skirt owns its vertices and carries the ground's
    // own normal, the material undoes the flip (see terrainMaterial), and the
    // wall is sunk under the ground edge by a fraction of the node's width -
    // under a pixel, since a node is drawn at most `splitPixels` wide - so that
    // it cannot win a depth comparison against the neighbor's surface either.
    // The dotted line along every node boundary was this.
    const skirtM = Math.min(TERRAIN_SKIRT_MAX_M, TERRAIN_SKIRT_RELIEF * Math.max(hMax - hMin, 1) + TERRAIN_SKIRT_SIZE * node.sizeKm * 1000);
    const sinkM = TERRAIN_SKIRT_SINK * node.sizeKm * 1000;
    // One edge sample: the top of the wall then its bottom, so the pair is
    // `k` and `k + 1`. Returns the index of the top.
    const pushSkirt = (i, j) => {
      const lat = latOf(i), lon = lonOf(j), h = heightOf(i, j);
      const fade = fadeOf(lat, lon);
      normalAt(i, j, lat, lon, tmpNormal);
      for (let k = 0; k < 2; k++) {
        bodyPoint(lat, lon, h - (k === 0 ? sinkM : skirtM), tmpPoint).sub(center);
        positions.push(tmpPoint.x, tmpPoint.y, tmpPoint.z);
        normals.push(tmpNormal.x, tmpNormal.y, tmpNormal.z);
        uvs.push(u0 + us * j / n, v0 - us * i / n);
        globeUvs.push((lon - globeLonLeft) / 360, (lat + 90) / 180);
        fades.push(fade);
      }
      return positions.length / 3 - 2;
    };
    const skirtEdge = (i0, j0, i1, j1) => {
      const steps = Math.max(Math.abs(i1 - i0), Math.abs(j1 - j0));
      const di = Math.sign(i1 - i0), dj = Math.sign(j1 - j0);
      let prev = pushSkirt(i0, j0);
      for (let m = 1; m <= steps; m++) {
        const cur = pushSkirt(i0 + di * m, j0 + dj * m);
        index.push(prev, cur, prev + 1, cur, cur + 1, prev + 1);
        prev = cur;
      }
    };
    skirtEdge(0, 0, 0, n); skirtEdge(n, 0, n, n); skirtEdge(0, 0, n, 0); skirtEdge(0, n, n, n);

    const geometry = new THREE.BufferGeometry();
    geometry.setAttribute('position', new THREE.Float32BufferAttribute(positions, 3));
    geometry.setAttribute('normal', new THREE.Float32BufferAttribute(normals, 3));
    geometry.setAttribute('uv', new THREE.Float32BufferAttribute(uvs, 2));
    geometry.setAttribute('aGlobeUv', new THREE.Float32BufferAttribute(globeUvs, 2));
    geometry.setAttribute('aFade', new THREE.Float32BufferAttribute(fades, 1));
    geometry.setIndex(index);
    geometry.computeBoundingSphere();

    const mesh = new THREE.Mesh(geometry, terrainMaterial());
    mesh.material.map = source ? source.texture : null;
    mesh.position.copy(center);
    mesh.name = `terrain-${node.level}-${node.x}-${node.y}`;
    // The tree spans the whole visible horizon while the sun's shadow camera is
    // sized to the vehicle: terrain casting into that camera puts its own
    // clamped shadow edge across the ground for kilometers, so the nodes only
    // receive. The lander still casts its shadow onto them.
    mesh.castShadow = false;
    mesh.receiveShadow = true;
    mesh.visible = false;
    mesh.frustumCulled = false;   // the tree culls; three's own test would use the node-local bounds
    node.mesh = mesh;
    node.texKey = source ? source.key : '';
    node.radiusKm = geometry.boundingSphere ? geometry.boundingSphere.radius : node.radiusKm;
    node.heightSpanM = hMax - hMin;
    group.add(mesh);
    builtCount++;
    return node;
  }

  // A drawn node follows its tile as it decodes: when a sharper ancestor (or
  // its own tile) has arrived, the UV rectangle is rewritten in place.
  function refreshTile(node, frameIndex) {
    const source = resolveTile(node.level, node.x, node.y, frameIndex);
    if (!source || source.key === node.texKey) return;
    const n = segments;
    const uv = node.mesh.geometry.getAttribute('uv');
    for (let i = 0; i <= n; i++) {
      for (let j = 0; j <= n; j++) uv.setXY(i * (n + 1) + j, source.u0 + source.size * j / n, source.v0 - source.size * i / n);
    }
    // the skirt vertices repeat the edge rows in the order they were pushed,
    // two of them (the top of the wall and its bottom) per edge sample
    let k = (n + 1) * (n + 1);
    const setSkirt = (i, j) => {
      const u = source.u0 + source.size * j / n, v = source.v0 - source.size * i / n;
      uv.setXY(k++, u, v); uv.setXY(k++, u, v);
    };
    const edge = (i0, j0, i1, j1) => {
      const steps = Math.max(Math.abs(i1 - i0), Math.abs(j1 - j0));
      const di = Math.sign(i1 - i0), dj = Math.sign(j1 - j0);
      for (let m = 0; m <= steps; m++) setSkirt(i0 + di * m, j0 + dj * m);
    };
    edge(0, 0, 0, n); edge(n, 0, n, n); edge(0, 0, n, 0); edge(0, n, n, n);
    uv.needsUpdate = true;
    node.mesh.material.map = source.texture;
    node.mesh.material.needsUpdate = true;
    node.texKey = source.key;
  }

  function childrenOf(node) {
    if (!node.children) {
      node.children = [
        nodeAt(node.level + 1, 2 * node.x, 2 * node.y),
        nodeAt(node.level + 1, 2 * node.x + 1, 2 * node.y),
        nodeAt(node.level + 1, 2 * node.x, 2 * node.y + 1),
        nodeAt(node.level + 1, 2 * node.x + 1, 2 * node.y + 1),
      ];
      for (const child of node.children) child.parent = node;
    }
    return node.children;
  }

  // ---- the per-frame walk -------------------------------------------------
  const rootNode = makeNode(0, 0, 0);
  buildNode(rootNode, 0);

  const terrainFrustum = new THREE.Frustum();
  const terrainMatrix = new THREE.Matrix4();
  const terrainInverse = new THREE.Matrix4();
  const terrainSphere = new THREE.Sphere();
  const camLocal = new THREE.Vector3();
  const camDir = new THREE.Vector3();
  const drawn = [];
  const queue = [];
  let frameIndex = 0;
  let horizonPlane = -Infinity, camRadius = 0, pixelScale = 1000;
  const stats = { nodes: 0, drawn: 0, triangles: 0, pending: 0, built: 1, updateMs: 0, maxLevel: 0 };

  function visible(node) {
    terrainSphere.center.copy(node.center);
    terrainSphere.radius = node.radiusKm * 1.05;
    if (!terrainFrustum.intersectsSphere(terrainSphere)) return false;
    // Below the horizon: the sphere of the covered region hides it from the camera.
    if (horizonPlane > -Infinity && node.center.dot(camDir) + node.radiusKm < horizonPlane) return false;
    return true;
  }

  function projectedPixels(node) {
    const d = Math.max(camLocal.distanceTo(node.center) - node.radiusKm, 1e-6);
    return node.sizeKm / d * pixelScale;
  }

  function show(node) {
    node.used = frameIndex;
    node.mesh.visible = true;
    drawn.push(node);
    stats.triangles += node.mesh.geometry.index.count / 3;
    if (node.level > stats.maxLevel) stats.maxLevel = node.level;
    refreshTile(node, frameIndex);
  }

  function walk(node) {
    node.used = frameIndex;
    if (!visible(node)) return;
    if (node.level < maxLevel && projectedPixels(node) > splitPixels) {
      const kids = childrenOf(node);
      let ready = true;
      for (const child of kids) {
        if (child.mesh) continue;
        ready = false;
        queue.push(child);
      }
      if (ready) {
        for (const child of kids) walk(child);
        return;
      }
    }
    show(node);
  }

  // Built nodes are kept for reuse; the least recently drawn go when there are
  // too many, so a long descent does not accumulate every node it passed.
  function evictNodes() {
    if (nodes.size <= TERRAIN_NODE_CAP) return;
    const built = [];
    for (const node of nodes.values()) if (node.mesh) built.push(node);
    if (built.length <= TERRAIN_NODE_CAP) return;
    built.sort((a, b) => a.used - b.used);
    for (let k = 0; k < built.length - TERRAIN_NODE_CAP * 0.8; k++) {
      const node = built[k];
      if (node === rootNode || node.used >= frameIndex - 1) break;
      group.remove(node.mesh);
      node.mesh.geometry.dispose();
      const at = materials.indexOf(node.mesh.material);
      if (at >= 0) materials.splice(at, 1);
      node.mesh.material.dispose();
      node.mesh = null;
      node.texKey = '';
      if (node.parent) node.parent.children = null;
      nodes.delete(tileKey(node.level, node.x, node.y));
    }
  }

  function update(camera, viewportHeight) {
    if (!group.visible) return;
    const t0 = (typeof performance !== 'undefined' ? performance.now() : 0);
    frameIndex++;
    terrainMatrix.multiplyMatrices(camera.projectionMatrix, camera.matrixWorldInverse).multiply(group.matrixWorld);
    terrainFrustum.setFromProjectionMatrix(terrainMatrix);
    terrainInverse.copy(group.matrixWorld).invert();
    camera.getWorldPosition(camLocal).applyMatrix4(terrainInverse);
    camRadius = camLocal.length();
    camDir.copy(camLocal).divideScalar(camRadius || 1);
    // The horizon plane of a sphere just under the lowest ground: points behind
    // it cannot be seen. Skip it when the camera is at or below that sphere.
    const rFloor = R - 5;
    horizonPlane = camRadius > rFloor * 1.0000001 ? (rFloor * rFloor) / camRadius : -Infinity;
    const fov = THREE.MathUtils.degToRad(camera.fov);
    pixelScale = (viewportHeight || 900) / (2 * Math.tan(0.5 * fov));

    for (const node of drawn) if (node.mesh) node.mesh.visible = false;
    drawn.length = 0;
    queue.length = 0;
    stats.triangles = 0;
    stats.maxLevel = 0;
    walk(rootNode);
    // A few node builds per frame: the camera moves smoothly, so the tree
    // catches up within a few frames and the frame itself never stalls.
    if (queue.length) {
      queue.sort((a, b) => projectedPixels(b) - projectedPixels(a));
      for (let k = 0; k < Math.min(buildBudget, queue.length); k++) buildNode(queue[k], frameIndex);
    }
    evictNodes();
    stats.nodes = nodes.size;
    stats.drawn = drawn.length;
    stats.pending = queue.length;
    stats.built = builtCount;
    stats.updateMs = (typeof performance !== 'undefined' ? performance.now() : 0) - t0;
  }

  // The globe is cut open under the covered region, a hair inside it so the
  // sphere and the terrain's faded outer ring overlap rather than leave a gap.
  const inset = 0.0025 * rootSpan;
  const hole = {
    lat_min: root.lat_min + inset, lat_max: root.lat_max - inset,
    lon_min: root.lon_min + inset, lon_max: root.lon_max - inset,
  };

  // Site marker: a thin ring on the ground.
  let siteMarker = null;
  if (spec.site) {
    const h = Number.isFinite(spec.site.height_m) ? spec.site.height_m : heightAt(spec.site.lat_deg, spec.site.lon_deg) || 0;
    const center = bodyPoint(spec.site.lat_deg, spec.site.lon_deg, h + 0.5, new THREE.Vector3());
    const ring = new THREE.Mesh(new THREE.RingGeometry(0.004, 0.005, 48), new THREE.MeshBasicMaterial({ color: 0xf2b950, side: THREE.DoubleSide, transparent: true, opacity: 0.8, depthWrite: false }));
    ring.position.copy(center);
    ring.lookAt(center.clone().multiplyScalar(2));
    ring.renderOrder = 20;
    ring.name = 'site-marker';
    ring.userData.uiOnly = true;
    group.add(ring);
    siteMarker = ring;
  }

  // The ground the vehicle lands on: the deepest tile that covers the site.
  // The lighting module sets the exposure from it (it reads the texture off
  // `levels[0].mesh.material.map`), and with a corridor-wide root the widest
  // tile is no longer that ground, so the site's own tile is reported instead.
  const siteTileEntry = (() => {
    if (!spec.site) return tileIndex.get(tileKey(0, 0, 0)) || null;
    let found = null;
    for (let level = 0; level <= maxTileLevel; level++) {
      const n = 1 << level;
      const x = Math.floor((spec.site.lon_deg - root.lon_min) / rootLonSpan * n);
      const y = Math.floor((root.lat_max - spec.site.lat_deg) / rootSpan * n);
      const entry = tileIndex.get(tileKey(level, x, y));
      if (entry) found = entry;
    }
    return found || tileIndex.get(tileKey(0, 0, 0)) || null;
  })();
  const siteGround = { material: { get map() { return siteTileEntry ? tileTexture(siteTileEntry, frameIndex) : null; } } };

  const finestTile = tiles.nodes.reduce((best, n) => (n.m_per_px > 0 && (!best || n.m_per_px < best) ? n.m_per_px : best), 0);
  return {
    group,
    // A quadtree has no fixed levels; this reports the ground under the site
    // for the lighting module's exposure and tells the rest of the page that
    // the run has terrain at all.
    levels: [{ mesh: siteGround, m_per_px: siteTileEntry ? siteTileEntry.m_per_px : 0 }],
    rootMesh: rootNode.mesh,
    hole,
    site: spec.site || null,
    referenceRadiusKm: R,
    heightAt,
    surfaceHeight,
    stats,
    root,
    maxLevel,
    siteMarker,
    update,
    setGlobeTexture,
    setVisible(v) { group.visible = v; },
    // The imagery's required credit, which the payload carries so a page that ships
    // archive data says where it came from.
    get attribution() { return (tiles.attribution || []).join(' \u00b7 '); },
    get modelStatus() {
      // the finest tiles are sampled finer than their source resolves, so the panel reports the
      // feature scale the payload states beside the sampling rather than the sampling alone
      const resolves = tiles.resolution && tiles.resolution.feature_scale_m;
      return `quadtree to L${maxTileLevel} (+${TERRAIN_EXTRA_LEVELS} relief), ${tiles.nodes.length} tiles, finest ${finestTile ? finestTile.toFixed(2) : '–'} m/px${resolves ? ` sampling of ~${resolves.toFixed(1)} m detail` : ''}, ${grids.length} grids (${finest.rows}x${finest.cols})`;
    },
  };
}
