import { test } from 'node:test';
import assert from 'node:assert/strict';
import { OM } from '../../engine/om.js';
import '../../engine/vegetation.js';
import cove from '../../world/cove.js';
import { plan, treeAt, TREE_LATTICE } from '../../world/cove-flora.js';
import { ForestField, FLOATS_PER_TREE, TILE_CELLS } from '../../engine/forest.js';
import { PROTOTYPES, prototypeFor, frame } from '../../engine/impostors.js';
// Engine systems read the active world; these tests run in the development cove.
OM.world = cove;
// The far field reads the landscape after the detailed plan has set its canopies.
plan();
const field = (range = 600) =>
  new ForestField(cove.vegetation, (x, z) => cove.landscape.height(x, z), { range });

test('the far field draws every lattice cell except the detailed planting, by the same rule', () => {
  const f = field(),
    { i0, i1, j0, j1 } = TREE_LATTICE.detailed;
  for (const [ti, tj] of [
    [0, 0],
    [-1, 0],
    [1, -1],
    [3, 2],
  ]) {
    const data = f.generate(ti, tj),
      expected = [];
    for (let j = tj * TILE_CELLS; j < (tj + 1) * TILE_CELLS; j++)
      for (let i = ti * TILE_CELLS; i < (ti + 1) * TILE_CELLS; i++) {
        if (i >= i0 && i <= i1 && j >= j0 && j <= j1) continue;
        const t = treeAt(i, j);
        if (t) expected.push(t);
      }
    assert.equal(data.length, expected.length * FLOATS_PER_TREE, `tile ${ti},${tj}`);
    expected.forEach((t, k) => {
      const o = k * FLOATS_PER_TREE;
      assert.equal(data[o], Math.fround(t.x));
      assert.equal(data[o + 2], Math.fround(t.z));
      assert.equal(data[o + 1], Math.fround(cove.landscape.height(t.x, t.z)));
      assert.equal(PROTOTYPES[data[o + 4]].family, t.family);
      assert.equal(PROTOTYPES[data[o + 4]].ageClass, t.ageClass);
      assert.ok(Math.abs(data[o + 3] * PROTOTYPES[data[o + 4]].height - t.height) < 1e-5);
    });
  }
});

test('the detailed plan is the same rule on its own cells', () => {
  const { i0, i1, j0, j1 } = TREE_LATTICE.detailed,
    candidates = [];
  for (let j = j0; j <= j1; j++)
    for (let i = i0; i <= i1; i++) {
      const t = treeAt(i, j);
      if (t) candidates.push(t.id);
    }
  // The plan then removes crowns that overlap more closely than crown spacing allows.
  const planted = plan().trees.filter((t) => !t.managed);
  for (const t of planted) assert.ok(candidates.includes(t.id), t.id);
});

test('tiles are chosen nearest first and generated within the budget', () => {
  const f = field(700);
  assert.equal(f.select(0, 9), true);
  assert.equal(f.select(1, 10), false);
  const d = f.desired.map((t) => t.d);
  assert.deepEqual(
    d,
    [...d].sort((a, b) => a - b),
  );
  assert.ok(d.every((v) => v <= 700));
  const first = f.work(0);
  assert.equal(first.length, 1);
  assert.equal(first[0].key, f.desired[0].key);
  f.work(Infinity);
  assert.equal(f.complete, true);
  const trees = [...f.selected()].reduce((n, a) => n + a.length / FLOATS_PER_TREE, 0);
  assert.ok(trees > 500, `${trees} trees within 700 m`);
});

test('generation is deterministic and moving within a tile keeps the selection', () => {
  const a = field(),
    b = field();
  assert.deepEqual(a.generate(-2, 3), b.generate(-2, 3));
  a.select(0, 9);
  const keys = a.desired.map((t) => t.key).join();
  a.select(20, 30);
  assert.equal(a.desired.map((t) => t.key).join(), keys);
  a.select(0, 9 + TILE_CELLS * TREE_LATTICE.spacing);
  assert.notEqual(a.desired.map((t) => t.key).join(), keys);
});

test('each tree has a prototype of its own family and age class, and the card holds the model', () => {
  for (const family of [0, 1, 2])
    for (const ageClass of ['mature', 'juvenile'])
      for (const seed of [1, 2, 3, 99991]) {
        const p = PROTOTYPES[prototypeFor({ family, ageClass, seed })];
        assert.equal(p.family, family);
        assert.equal(p.ageClass, ageClass);
      }
  for (const spec of PROTOTYPES) {
    const model = OM.Vegetation.skeleton({ x: 0, y: 0, z: 0, ...spec }),
      f = frame(model);
    for (const l of model.leaves) {
      assert.ok(Math.hypot(l.position[0], l.position[2]) <= f.size / 2);
      assert.ok(l.position[1] >= f.bottom && l.position[1] <= f.bottom + f.size);
    }
  }
});

test('distant woodland cover follows the same rule, and the terrain carries it', async () => {
  const L = cove.landscape,
    at = (x, z) => {
      const h = L.height(x, z),
        g = L.gradient(x, z, 1),
        s = Math.hypot(g.x, g.z);
      return cove.vegetation.cover(L.classify(x, z, h, s, 0.24, 0, 0, 0.8), h, s);
    };
  assert.equal(at(0, -600), 0, 'open water carries no woodland');
  const inland = at(0, 4000);
  assert.ok(inland > 0.05 && inland < 0.3, `inland cover ${inland}`);
  const T = await import('three'),
    { TerrainSystem } = await import('../../engine/terrain.js'),
    terrain = new TerrainSystem(T, new T.Scene(), new T.MeshBasicMaterial());
  const outer = terrain.levels.at(-1).geo.attributes.omCanopy.array;
  assert.ok(outer.some((v) => v > 0.05));
  assert.ok(outer.every((v) => v >= 0 && v <= 1));
});
