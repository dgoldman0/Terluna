import { test, after } from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import path from 'node:path';
import crypto from 'node:crypto';
import { fileURLToPath } from 'node:url';
import * as T from 'three';
import L from '../../world/landscape.js';
import C from '../../engine/core.js';
import { TerrainSystem, triangle } from '../../engine/terrain.js';
import E from '../../engine/ecology.js';
import SurfaceWater from '../../engine/surface-water.js';
const __dirname = path.dirname(fileURLToPath(import.meta.url));
const root = path.join(__dirname, '..', '..'),
  report = {
    schema: 'open-moon-landscape-numeric/1',
    baseCommit: '4c7affae780d3ff7ef5581e88aac7ae32ada7505',
    checks: {},
  };
const close = (a, b, e = 1e-8) => assert.ok(Math.abs(a - b) <= e, `${a} != ${b}; tolerance ${e}`);
const scene = new T.Scene(),
  material = new T.MeshBasicMaterial();
let terrain, layout;
const terrainSystem = () => (terrain ||= new TerrainSystem(T, scene, material));
const ecology = () => (layout ||= E.plan());
const hash = (b) => crypto.createHash('sha256').update(b).digest('hex');
function renderedHeight(system, x, z) {
  for (const l of system.levels) {
    const b = l.bounds;
    if (x < b.minX || x >= b.maxX || z < b.minZ || z >= b.maxZ) continue;
    const i = Math.floor((x - b.minX) / l.step),
      j = Math.floor((z - b.minZ) / l.step),
      n = system.cells + 1,
      k = j * n + i,
      p = l.geo.attributes.position.array;
    return triangle(
      p[k * 3 + 1],
      p[(k + 1) * 3 + 1],
      p[(k + n) * 3 + 1],
      p[(k + n + 1) * 3 + 1],
      (x - b.minX) / l.step - i,
      (z - b.minZ) / l.step - j,
    );
  }
  throw Error('Point outside terrain');
}
function topology(system) {
  const ids = new Map(),
    vertices = [],
    edges = new Map();
  let triangles = 0,
    maxWeldError = 0,
    minArea = Infinity;
  const vertex = (p, k) => {
    const x = p[k * 3],
      y = p[k * 3 + 1],
      z = p[k * 3 + 2],
      key = (Math.round(x * 4) + 1048576) * 2097153 + (Math.round(z * 4) + 1048576);
    let id = ids.get(key);
    if (id === undefined) {
      id = vertices.length;
      ids.set(key, id);
      vertices.push([x, y, z]);
    } else {
      const error = Math.abs(vertices[id][1] - y);
      maxWeldError = Math.max(maxWeldError, error);
      assert.ok(error < 0.00004, `LOD height separation at ${x},${z}: ${error}`);
    }
    return id;
  };
  for (const l of system.levels) {
    const p = l.geo.attributes.position.array,
      idx = l.geo.index.array;
    for (let q = 0; q < l.geo.drawRange.count; q += 3) {
      const a = idx[q],
        b = idx[q + 1],
        c = idx[q + 2];
      assert.ok(a < p.length / 3 && b < p.length / 3 && c < p.length / 3);
      const area =
        (p[b * 3 + 2] - p[a * 3 + 2]) * (p[c * 3] - p[a * 3]) -
        (p[b * 3] - p[a * 3]) * (p[c * 3 + 2] - p[a * 3 + 2]);
      assert.ok(area > 0);
      minArea = Math.min(minArea, area);
      triangles++;
      const v = [vertex(p, a), vertex(p, b), vertex(p, c)];
      for (let k = 0; k < 3; k++) {
        let a = v[k],
          b = v[(k + 1) % 3];
        if (a > b) [a, b] = [b, a];
        const key = a * 1000000 + b;
        edges.set(key, (edges.get(key) || 0) + 1);
      }
    }
  }
  const outer = system.levels.at(-1).bounds;
  let boundary = 0;
  for (const [key, n] of edges) {
    assert.ok(n === 1 || n === 2, `Non-manifold edge multiplicity ${n}`);
    if (n === 1) {
      boundary++;
      const a = vertices[Math.floor(key / 1000000)],
        b = vertices[key % 1000000];
      assert.ok(
        (a[0] === outer.minX && b[0] === outer.minX) ||
          (a[0] === outer.maxX && b[0] === outer.maxX) ||
          (a[2] === outer.minZ && b[2] === outer.minZ) ||
          (a[2] === outer.maxZ && b[2] === outer.maxZ),
        `Interior crack from ${a} to ${b}`,
      );
    }
  }
  assert.equal(boundary, system.cells * 4);
  return {
    triangles,
    boundaryEdges: boundary,
    weldedVertices: vertices.length,
    maxWeldError_m: maxWeldError,
    minTwiceTriangleArea_m2: minArea,
  };
}
test('Authoritative height initializes once and preserves the exact M1 reference eye height', () => {
  const h = L.height(0, 9);
  L.initialize();
  close(h, L.height(0, 9), 0);
  close(h, 1.3458966316214787, 1e-12);
  close(C.groundHeight(0, 9) + 1.7, 3.045873320927338, 1e-12);
});
test('Landscape inputs reject nonfinite coordinates', () => {
  assert.throws(() => L.height(NaN, 0));
  assert.throws(() => L.height(0, Infinity));
});
test('Actual zero-contour roots replace the coastal guide, including extended headlands', () => {
  let maximum = 0;
  for (let x = -500; x <= 500; x += 5) {
    const z = L.shorelineAtX(x);
    maximum = Math.max(maximum, Math.abs(L.height(x, z)));
  }
  assert.ok(maximum < 1e-7);
  report.checks.shorelineMaxElevationResidual_m = maximum;
});
test('Elevation and its gradient stay continuous across the field-grid perimeter', () => {
  let max = 0;
  for (const z of [-128, -127, 0, 50, 372, 384])
    for (const x of [-256, -244, -30, 0, 244, 256]) {
      const eps = 1e-5,
        d = Math.abs(L.height(x + eps, z) - L.height(x - eps, z));
      max = Math.max(max, d);
      assert.ok(d < 0.0001);
    }
  report.checks.boundaryHeightDifference_m = max;
});
test('Material weights are normalized and nonnegative over sea, coast and upland', () => {
  for (let z = -120; z <= 160; z += 7)
    for (let x = -140; x <= 140; x += 7) {
      const f = L.sample(x, z);
      close(
        f.weights.reduce((a, b) => a + b),
        1,
        1e-6,
      );
      assert.ok(f.weights.every((v) => v >= 0 && v <= 1));
      assert.ok(f.soilDepth > 0 && f.soilDepth <= 1.200001);
      assert.ok(f.moisture >= 0 && f.moisture <= 1);
    }
});
test('CPU material placement reads the same raster samples uploaded to the GPU', () => {
  const bytes = L.materialFieldRGBA();
  for (let z = -100; z < 150; z += 13.3)
    for (let x = -100; x < 130; x += 12.7) {
      const f = L.sample(x, z);
      for (let k = 0; k < 4; k++) close(f.weights[k], L.gridChannel(bytes, x, z, k), 2e-7);
    }
});
test('Drainage processes every cell and exports exactly its contributing area', () => {
  const d = L.data.hydrology;
  assert.equal(d.processed, L.GRID.n ** 2);
  close(d.outletArea, d.totalArea, 1e-8);
  report.checks.drainage = {
    cells: d.processed,
    totalArea_m2: d.totalArea,
    outletArea_m2: d.outletArea,
  };
});
test('Receivers descend the filled surface and their graph is acyclic', () => {
  const d = L.data.hydrology,
    order = new Int32Array(d.order.length);
  d.order.forEach((id, i) => (order[id] = i));
  for (let i = 0; i < d.receiver.length; i++) {
    const r = d.receiver[i];
    if (r < 0) continue;
    assert.ok(d.filled[r] <= d.filled[i] + 1e-12);
    assert.ok(order[r] > order[i]);
    assert.ok(d.accumulation[r] >= d.accumulation[i]);
  }
});
test('Priority flood fills an isolated bowl to its explicit spill elevation', () => {
  const n = 7,
    h = new Float64Array(n * n).fill(3);
  h[3 * n + 3] = 1;
  h[3 * n + 2] = 1.5;
  const d = L.drainageGrid(h, n, 2);
  close(d.filled[3 * n + 3], 3);
  close(d.outletArea, d.totalArea);
  assert.equal(d.processed, n * n);
});
test('Drainage rejects malformed grids', () => {
  assert.throws(() => L.drainageGrid(new Float64Array(5), 3, 1));
  assert.throws(() => L.drainageGrid(new Float64Array(9), 3, 0));
});
test('Eleven moving grids retain finite vertices and upward unit normals', () => {
  const t = terrainSystem();
  assert.equal(t.levelCount, 11);
  assert.equal(t.stats.outerHalfExtent_m, 16384);
  assert.equal(t.stats.nearSpacing_m, 0.25);
  for (const l of t.levels) {
    const p = l.geo.attributes.position,
      n = l.geo.attributes.normal;
    for (let i = 0; i < p.count; i++) {
      assert.ok([p.getX(i), p.getY(i), p.getZ(i)].every(Number.isFinite));
      close(Math.hypot(n.getX(i), n.getY(i), n.getZ(i)), 1, 2e-6);
      assert.ok(n.getY(i) > 0);
    }
  }
  report.checks.initialTerrain = t.stats;
});
test('Every initial terrain edge is paired except the true outer perimeter', () => {
  report.checks.initialTopology = topology(terrainSystem());
});
test('Moving through three distinct grid-snapping states preserves the closed topology', () => {
  report.checks.movingTopology = [];
  for (const [x, z] of [
    [0.51, 9.51],
    [3.51, 12.01],
    [-38.1, 47.7],
  ]) {
    terrainSystem().update(x, z);
    report.checks.movingTopology.push({ observer: [x, z], ...topology(terrainSystem()) });
  }
});
test('The fine terrain follows the observer and matches walking height within twelve millimetres', () => {
  const t = terrainSystem(),
    r = C.rng(554);
  let max = 0;
  for (const [cx, cz] of [
    [0, 9],
    [20, 33],
    [-30, 43],
    [-56, 47],
  ]) {
    t.update(cx, cz);
    assert.ok(Math.abs(t.levels[0].bounds.cx - cx) < 0.5);
    assert.ok(Math.abs(t.levels[0].bounds.cz - cz) < 0.5);
    for (let i = 0; i < 180; i++) {
      const x = cx + (r() - 0.5) * 18,
        z = cz + (r() - 0.5) * 18;
      max = Math.max(max, Math.abs(renderedHeight(t, x, z) - C.groundHeight(x, z)));
    }
  }
  assert.ok(max < 0.012, `Collision mesh error ${max} m`);
  report.checks.nearCollisionMaxError_m = max;
});
test('Crossing a snap leaves fixed near-observer samples unchanged', () => {
  const t = terrainSystem(),
    points = Array.from({ length: 20 }, (_, i) => [0.2 + i * 0.1, 9.7 + i * 0.03]);
  t.update(0.49, 9.49);
  const before = points.map((p) => renderedHeight(t, ...p));
  t.update(0.51, 9.51);
  points.forEach((p, i) => close(renderedHeight(t, ...p), before[i], 1e-9));
});
test('An unchanged observer issues no geometry update or new height evaluation', () => {
  const t = terrainSystem();
  t.update(4, 12);
  const version = t.stats.updates,
    count = t.stats.heightEvaluations;
  assert.equal(t.update(4, 12), false);
  assert.equal(t.stats.updates, version);
  assert.equal(t.stats.heightEvaluations, count);
});
test('Terrain caches and topology templates remain bounded during traversal', () => {
  const local = new TerrainSystem(T, new T.Scene(), material, 'moon', { cells: 32, levels: 5 });
  for (let i = 0; i < 70; i++) local.update(-40 + i * 0.77, 45 + i * 0.11);
  for (const l of local.levels) {
    assert.ok(l.cache.size <= 33 * 33 * 5);
    assert.ok(l.topologies.size <= 4);
  }
  report.checks.traversalCacheSamples = local.stats.cachedSamples;
  local.dispose();
});
test('Moving-grid options reject incompatible topology and nonfinite observers', () => {
  assert.throws(() => new TerrainSystem(T, new T.Scene(), material, 'moon', { cells: 17 }));
  assert.throws(() => terrainSystem().update(Infinity, 0));
});
test('Earth comparison changes only the spherical sag in the terrain geometry', () => {
  const t = terrainSystem();
  t.setWorld('earth');
  t.update(0, 9);
  const l = t.levels[0],
    p = l.geo.attributes.position;
  for (let i = 100; i < p.count; i += 317) {
    const x = p.getX(i),
      z = p.getZ(i);
    if (Math.abs(x) < 10 && Math.abs(z - 9) < 10)
      close(p.getY(i), C.groundHeight(x, z, 'earth'), 2e-6);
  }
  t.setWorld('moon');
  t.update(0, 9);
});
test('Plant, rock and understory identities are reproducible after canopy accumulation', () => {
  const a = ecology(),
    b = E.plan();
  assert.deepEqual(a, b);
  report.checks.ecology = {
    trees: a.trees.length,
    rocks: a.rocks.length,
    gravel: a.gravel.length,
    tufts: a.tufts.length,
    understory: a.understory.length,
    managedTrees: a.trees.filter((t) => t.managed).map((t) => t.id),
  };
});
test('Habitat placement excludes water, the shelter and the walking corridor', () => {
  for (const t of ecology().trees) {
    assert.ok(L.height(t.x, t.z) > (t.managed ? 0.7 : 1.2));
    assert.equal(C.roofMask(t.x, t.z, 4), 0);
    assert.ok(C.pathDistance(t.x, t.z) > 3.7);
    assert.ok(t.habitat.soilDepth > 0.16);
  }
  for (const r of ecology().rocks) assert.ok(C.pathDistance(r.x, r.z) >= 2 + r.s * 1.2);
});
test('The planted canopy changes the shared environmental field and understory light', () => {
  const trees = ecology().trees;
  assert.ok(trees.length > 10);
  const t = trees.find((t) => !t.managed),
    f = L.sample(t.x, t.z);
  assert.ok(f.canopy > 0.5);
  assert.ok(f.lightAvailability < 0.65);
  assert.ok(L.environmentFieldRGBA().some((v, i) => i % 4 === 1 && v > 0.5));
});
test('Tree families and both age classes are represented', () => {
  assert.deepEqual([...new Set(ecology().trees.map((t) => t.family))].sort(), [0, 1, 2]);
  assert.deepEqual([...new Set(ecology().trees.map((t) => t.ageClass))].sort(), [
    'juvenile',
    'mature',
  ]);
});
test('Branch hierarchy is connected at real trunk/limb endpoints and radii taper', () => {
  let branches = 0,
    leaves = 0;
  for (const tree of ecology().trees) {
    const s = E.skeleton(tree),
      nodes = new Set([JSON.stringify([tree.x, tree.y, tree.z])]);
    for (const b of s.segments) {
      assert.ok(nodes.has(JSON.stringify(b.a)), `Disconnected branch in ${tree.id}`);
      assert.ok(b.r0 >= b.r1 && b.r1 > 0);
      assert.ok(b.a.concat(b.b).every(Number.isFinite));
      nodes.add(JSON.stringify(b.b));
    }
    branches += s.segments.length;
    leaves += s.leaves.length;
  }
  report.checks.treeArchitecture = { branchSegments: branches, leaves };
});
test('Every leaf has a narrow geometric petiole and valid finite lamina vertices', () => {
  for (let family = 0; family < 3; family++) {
    const g = E.leafGeometry(T, family),
      p = g.attributes.position;
    assert.ok(p.count >= 14);
    assert.equal(Math.min(...Array.from({ length: p.count }, (_, i) => p.getY(i))), 0);
    assert.ok(p.array.every(Number.isFinite));
    g.dispose();
  }
});
test('Surface-water grids reject invalid dimensions and steps', () => {
  assert.throws(() => new SurfaceWater({ n: 2 }));
  assert.throws(() => new SurfaceWater({ cell: 0 }));
  const w = new SurfaceWater({ n: 5 });
  assert.throws(() => w.advance(-1, C.weatherAt(0)));
});
test('Clear weather has dry exposed surfaces and conserved antecedent soil water', () => {
  const w = new SurfaceWater();
  assert.ok(w.ledger.initialStorage > 0);
  w.seek(14400, 'clear');
  for (const key of ['leaf', 'film', 'pond']) assert.ok(w[key].every((v) => v === 0));
  assert.ok(w.soil.some((v) => v > 0));
  assert.ok(w.ledger.relativeResidual < 1e-11);
});
test('Antecedent soil water is initialized from the shared habitat-moisture field', () => {
  const w = new SurfaceWater({ n: 17, minX: -40, minZ: -12 });
  for (let j = 0; j < w.n; j++)
    for (let i = 0; i < w.n; i++) {
      const k = j * w.n + i,
        f = L.sample(w.minX + i * w.cell, w.minZ + j * w.cell);
      close(w.initialSoil[k], f.elevation > 0 ? w.capacity[k] * f.moisture * 0.65 : 0, 1e-10);
    }
});
test('Rainfall, four stores and explicit exports conserve water across the full episode', () => {
  const w = new SurfaceWater();
  w.seek(14400, 'episode');
  const r = w.ledger;
  assert.ok(r.input > 0 && r.storage > 0 && r.roofExport > 0 && r.boundaryExport > 0);
  assert.ok(r.relativeResidual < 1e-11);
  for (const key of ['leaf', 'film', 'soil', 'pond'])
    assert.ok(w[key].every((v) => Number.isFinite(v) && v >= -1e-12));
  for (let i = 0; i < w.soil.length; i++) {
    assert.ok(w.soil[i] <= w.capacity[i] + 1e-9);
    assert.ok(w.leaf[i] <= 0.35 * w.cover[i] + 1e-9);
  }
  report.checks.episodeWater = r;
});
test('Surface-water replay is deterministic after a time reset', () => {
  const w = new SurfaceWater({ n: 33, minX: -64, minZ: -24 });
  w.seek(8000, 'episode');
  const stores = ['leaf', 'film', 'soil', 'pond'].map((k) => Array.from(w[k])),
    ledger = w.ledger;
  w.seek(0, 'clear');
  w.seek(8000, 'episode');
  assert.deepEqual(
    stores,
    ['leaf', 'film', 'soil', 'pond'].map((k) => Array.from(w[k])),
  );
  assert.deepEqual(ledger, w.ledger);
});
test('Forward seeks on the same integration boundaries equal direct replay', () => {
  const a = new SurfaceWater({ n: 17 }),
    b = new SurfaceWater({ n: 17 });
  a.seek(4000, 'episode').seek(8000, 'episode');
  b.seek(8000, 'episode');
  assert.deepEqual(a.ledger, b.ledger);
  assert.deepEqual(a.textureData(), b.textureData());
});
test('Soil and pond stores persist after rain; post-episode drying conserves mass', () => {
  const w = new SurfaceWater({ n: 33, minX: -64, minZ: -16 });
  w.seek(14400, 'episode');
  const initial = w.ledger.storage;
  assert.ok(w.soil.some((v) => v > 0));
  w.advance(3600, C.weatherAt(14400));
  assert.ok(w.ledger.storage < initial);
  assert.ok(w.ledger.relativeResidual < 1e-11);
  assert.equal(w.ledger.elapsedSeconds, 18000);
});
test('Five- and ten-second substeps remain close under the same rain forcing', () => {
  const a = new SurfaceWater({ n: 17, minX: -40, minZ: -12 }),
    b = new SurfaceWater({ n: 17, minX: -40, minZ: -12 });
  const w = { ...C.weatherAt(7000), rain: 8 };
  for (let t = 0; t < 1800; t += 10) a.advance(10, w);
  for (let t = 0; t < 1800; t += 5) b.advance(5, w);
  const relative = Math.abs(a.ledger.storage - b.ledger.storage) / Math.max(1, b.ledger.storage);
  assert.ok(relative < 0.02);
  report.checks.waterStepStorageDifferenceFraction = relative;
});
test('Rendered channels reflect distinct stores and reset restores their initial values', () => {
  const w = new SurfaceWater({ n: 17, minX: -40, minZ: -12 });
  const initial = Array.from(w.textureData());
  w.seek(8500, 'episode');
  const wet = w.textureData();
  assert.ok(wet.some((v) => v > 0));
  w.seek(0, 'clear');
  const dry = w.textureData();
  assert.notEqual(wet, dry);
  assert.deepEqual(Array.from(dry), initial);
});
test('Material asset hashes, encoding and physical scales match their manifest', () => {
  const m = JSON.parse(fs.readFileSync(path.join(root, 'assets/surfaces/manifest.json')));
  assert.equal(m.layers.length, 6);
  for (const [name, f] of Object.entries(m.files)) {
    const b = fs.readFileSync(path.join(root, 'assets/surfaces', name));
    assert.equal(b.length, f.bytes);
    assert.equal(hash(b), f.sha256);
  }
  for (const l of m.layers) {
    assert.ok(l.tile_metres > 0 && l.height_range_metres > 0);
    assert.ok(l.mean_albedo_linear.every((x) => x >= 0 && x <= 1));
  }
});
test('Independent array mipmaps and shared surface stores are wired into PBR shading', () => {
  const m = fs.readFileSync(path.join(root, 'engine/materials.js'), 'utf8');
  assert.ok(m.includes('new T.DataArrayTexture'));
  assert.ok(m.includes('sampler2DArray'));
  assert.ok(m.includes('roughnessFactor=mix(microRough'));
  assert.ok(m.includes('uSurfaceState'));
  const s = fs.readFileSync(path.join(root, 'engine/scene.js'), 'utf8');
  assert.ok(s.includes('image.data = state.surfaceWater.textureData()'));
});
test('Warm walking updates are measured separately from initial construction', () => {
  const t = terrainSystem(),
    times = [];
  for (let i = 0; i < 30; i++) {
    const start = performance.now();
    t.update(i * 0.51, 9 + i * 0.08);
    times.push(performance.now() - start);
  }
  times.sort((a, b) => a - b);
  report.checks.cpuTerrainUpdates = {
    samples: times.length,
    median_ms: times[Math.floor(times.length * 0.5)],
    p95_ms: times[Math.floor(times.length * 0.95)],
    maximum_ms: times.at(-1),
    context: 'Node CPU geometry updates; excludes GPU work and makes no hardware frame-rate claim',
  };
});
after(() => {
  if (terrain) terrain.dispose();
  material.dispose();
  fs.mkdirSync(path.join(root, 'test-results'), { recursive: true });
  fs.writeFileSync(
    path.join(root, 'test-results/landscape-numeric.json'),
    JSON.stringify(report, null, 2) + '\n',
  );
});
