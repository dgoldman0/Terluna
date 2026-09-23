/* Export a world for the Unreal game: heightmap, material layers, vegetation, stones
 * and sites, from the same world definition the web experience uses.
 *
 * Usage: node bake/game/world.mjs --out <dir> [--size 8129] [--spacing 2]
 *
 * Coordinates. The world is metres, right-handed, y up (three.js). Unreal is
 * centimetres, left-handed, z up. Unreal X = x, Unreal Y = z and Unreal Z = y keeps
 * every shape (swapping two axes also swaps the handedness), so files here stay in
 * world metres and manifest.json records the mapping and the landscape transform.
 *
 * Heights are the uncurved landscape: the Moon's curvature is for the renderer to
 * apply, as the web engine does. Sea level is 0 m.
 */
import fs from 'node:fs';
import path from 'node:path';
import crypto from 'node:crypto';
import { fileURLToPath } from 'node:url';
import { OM } from '../../engine/om.js';
import cove from '../../world/cove.js';
import { TREE_LATTICE } from '../../world/cove-flora.js';

const here = path.dirname(fileURLToPath(import.meta.url)),
  immersion = path.resolve(here, '../..');
const arg = (name, fallback) => {
  const i = process.argv.indexOf(name);
  return i > 0 ? process.argv[i + 1] : fallback;
};
const out = arg('--out');
if (!out) throw new Error('Give the output folder: --out <dir>');
// 8129 vertices = 32 x 32 Unreal landscape components of 254 quads.
const SIZE = +arg('--size', 8129),
  SPACING = +arg('--spacing', 2),
  LAYER_STEP = 2; // material layers are computed every LAYER_STEP samples and interpolated
const world = cove,
  L = world.landscape;
OM.world = world;
L.initialize();
const dir = path.resolve(out, 'world');
fs.mkdirSync(dir, { recursive: true });
const half = ((SIZE - 1) * SPACING) / 2,
  x0 = -half,
  z0 = -half;
const log = (m) => console.log(`${(performance.now() / 1000).toFixed(1)}s ${m}`);

// Heights.
const heights = new Float32Array(SIZE * SIZE);
let hmin = Infinity,
  hmax = -Infinity;
for (let j = 0; j < SIZE; j++)
  for (let i = 0; i < SIZE; i++) {
    const h = L.height(x0 + i * SPACING, z0 + j * SPACING);
    heights[j * SIZE + i] = h;
    if (h < hmin) hmin = h;
    if (h > hmax) hmax = h;
  }
log(`heights ${hmin.toFixed(2)} to ${hmax.toFixed(2)} m`);
// Unreal: height (cm) = location z + (value - 32768) * scale z / 128.
const mid = (hmin + hmax) / 2,
  scaleZ = Math.ceil(((hmax - hmin) * 100 * 128) / 65000),
  r16 = Buffer.alloc(SIZE * SIZE * 2);
for (let k = 0; k < heights.length; k++) {
  const v = Math.round(32768 + ((heights[k] - mid) * 100 * 128) / scaleZ);
  r16.writeUInt16LE(Math.max(0, Math.min(65535, v)), k * 2);
}
fs.writeFileSync(path.join(dir, 'height.r16'), r16);
const quantum = scaleZ / 128 / 100;

// Material layers from the landscape's habitat classification (sand, gravel, rock, soil).
const LAYERS = ['sand', 'gravel', 'rock', 'soil'],
  n = Math.ceil((SIZE - 1) / LAYER_STEP) + 1,
  coarse = LAYERS.map(() => new Float32Array(n * n));
for (let j = 0; j < n; j++)
  for (let i = 0; i < n; i++) {
    const f = L.sample(x0 + i * LAYER_STEP * SPACING, z0 + j * LAYER_STEP * SPACING);
    f.weights.forEach((w, k) => (coarse[k][j * n + i] = w));
  }
LAYERS.forEach((name, k) => {
  const bytes = Buffer.alloc(SIZE * SIZE);
  for (let j = 0; j < SIZE; j++)
    for (let i = 0; i < SIZE; i++) {
      const u = i / LAYER_STEP,
        v = j / LAYER_STEP,
        a = Math.min(n - 2, Math.floor(u)),
        b = Math.min(n - 2, Math.floor(v)),
        tu = u - a,
        tv = v - b,
        c = coarse[k],
        w =
          (c[b * n + a] * (1 - tu) + c[b * n + a + 1] * tu) * (1 - tv) +
          (c[(b + 1) * n + a] * (1 - tu) + c[(b + 1) * n + a + 1] * tu) * tv;
      bytes[j * SIZE + i] = Math.round(255 * Math.max(0, Math.min(1, w)));
    }
  fs.writeFileSync(path.join(dir, `layer_${name}.r8`), bytes);
});
log('material layers');

// Vegetation: the detailed plan, then the same rule on every other lattice cell in the extent.
const plan = world.flora(),
  { spacing: s, originX, originZ, detailed } = TREE_LATTICE,
  rows = ['id,x,z,surface_height_m,height_m,family,age_class,crown_radius_m,seed,detailed'];
const tree = (t, isDetailed) =>
  rows.push(
    [
      t.id,
      t.x,
      t.z,
      L.height(t.x, t.z),
      t.height,
      t.family,
      t.ageClass,
      t.crownRadius,
      t.seed,
      isDetailed ? 1 : 0,
    ]
      .map((v) => (typeof v === 'number' && !Number.isInteger(v) ? +v.toFixed(4) : v))
      .join(','),
  );
for (const t of plan.trees) tree(t, true);
const cells = Math.floor(half / s);
for (let j = -cells; j <= cells; j++)
  for (let i = -cells; i <= cells; i++) {
    if (i >= detailed.i0 && i <= detailed.i1 && j >= detailed.j0 && j <= detailed.j1) continue;
    const t = world.vegetation.treeAt(i, j);
    if (t && Math.abs(t.x - originX) < half && Math.abs(t.z - originZ) < half) tree(t, false);
  }
fs.writeFileSync(path.join(dir, 'trees.csv'), rows.join('\n') + '\n');
log(`${rows.length - 1} trees`);
const stones = ['id,x,z,surface_height_m,size_m,angle_rad,family,stretch'];
for (const r of plan.rocks)
  stones.push([r.id, r.x, r.z, L.height(r.x, r.z), r.s, r.angle, r.family, r.stretch].join(','));
fs.writeFileSync(path.join(dir, 'stones.csv'), stones.join('\n') + '\n');
fs.writeFileSync(
  path.join(dir, 'sites.json'),
  JSON.stringify({ waypoints: world.waypoints, shelter: world.shelter }, null, 1),
);

// Manifest: sources, the coordinate mapping and the landscape transform for import.
const sha = (file) => crypto.createHash('sha256').update(fs.readFileSync(file)).digest('hex');
const sources = [
  'world/cove.js',
  'world/landscape.js',
  'world/cove-flora.js',
  'world/cove-sites.js',
  'bake/game/world.mjs',
].map((f) => ({ file: `immersion/${f}`, sha256: sha(path.join(immersion, f)) }));
const files = fs
  .readdirSync(dir)
  .filter((f) => f !== 'manifest.json')
  .map((f) => ({ file: f, sha256: sha(path.join(dir, f)) }));
const manifest = {
  schema: 'terluna.game.world/1',
  world: { id: world.id, status: world.status },
  evidence:
    'The development cove: an authored placeholder landscape and planting rule for building ' +
    'the engine, not a proposed Open Moon environment.',
  coordinates:
    'Files are in world metres (x east-west, z north-south, y up). Unreal X = x, Y = z, ' +
    'Z = y, in centimetres. Heights are uncurved; sea level is 0 m.',
  landscape: {
    vertices: SIZE,
    spacing_m: SPACING,
    origin_m: { x: x0, z: z0 },
    height_range_m: [hmin, hmax],
    height_quantum_m: quantum,
    heightmap:
      'height.r16 (uint16 little-endian, row j = z0 + j * spacing, column i = x0 + i * spacing)',
    unreal_transform: {
      location_cm: [x0 * 100, z0 * 100, mid * 100],
      scale: [SPACING * 100, SPACING * 100, scaleZ],
    },
    layers: LAYERS.map((name) => ({ name, file: `layer_${name}.r8` })),
  },
  vegetation: { trees: rows.length - 1, file: 'trees.csv', lattice_spacing_m: s },
  sources,
  files,
};
fs.writeFileSync(path.join(dir, 'manifest.json'), JSON.stringify(manifest, null, 1));
log(`wrote ${dir}`);
