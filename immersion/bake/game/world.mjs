/* Export a world for the Unreal game: heightmap, material layers, vegetation, stones
 * and sites, from the same world definition the web experience uses.
 *
 * Usage: node bake/game/world.mjs --out <dir> [--size 8129] [--spacing 2] [--flat]
 *                                          [--far-size 2049] [--far-spacing 32]
 *
 * Coordinates. The world is metres, right-handed, y up (three.js). Unreal is
 * centimetres, left-handed, z up. Unreal X = x, Unreal Y = z and Unreal Z = y keeps
 * every shape (swapping two axes also swaps the handedness), so files here stay in
 * world metres and manifest.json records the mapping and the landscape transform.
 *
 * Heights include the Moon's curvature: the ground of a sphere of the Moon's radius
 * touching the world origin, which is where Unreal's SkyAtmosphere puts its planet
 * by default. Terrain, collision and sky then agree, and the horizon is where it
 * should be. --flat writes the uncurved landscape instead. Sea level is 0 m at the
 * origin and follows the sphere.
 *
 * Beyond the landscape, a coarser far grid (--far-size vertices at --far-spacing metres,
 * centred on the origin) carries the same height field, ground layers and the planting
 * rule's expected canopy cover out past the horizon, for distant terrain drawn without
 * individual trees. The surface textures (assets/surfaces) are copied alongside.
 */
import fs from 'node:fs';
import path from 'node:path';
import crypto from 'node:crypto';
import { fileURLToPath } from 'node:url';
import { OM } from '../../engine/om.js';
import C from '../../engine/core.js';
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
const CURVED = !process.argv.includes('--flat'),
  ground = (x, z) => (CURVED ? C.groundHeight(x, z, 'moon') : L.height(x, z));
const SIZE = +arg('--size', 8129),
  SPACING = +arg('--spacing', 2),
  LAYER_STEP = 2, // material layers are computed every LAYER_STEP samples and interpolated
  FAR_SIZE = +arg('--far-size', 2049),
  FAR_SPACING = +arg('--far-spacing', 32);
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
    const h = ground(x0 + i * SPACING, z0 + j * SPACING);
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
// The same layers at their computed spacing, interleaved: one RGBA texel per LAYER_STEP.
const packed = Buffer.alloc(n * n * 4);
for (let k = 0; k < n * n; k++)
  LAYERS.forEach((_, c) => (packed[k * 4 + c] = Math.round(255 * Math.max(0, Math.min(1, coarse[c][k])))));
fs.writeFileSync(path.join(dir, 'layers.rgba8'), packed);
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

// The far grid: heights (float32 metres, curved like the landscape's), and per vertex the
// ground layers sand, gravel, rock (soil is the rest) and the expected canopy cover.
const farHalf = ((FAR_SIZE - 1) * FAR_SPACING) / 2,
  farHeights = new Float32Array(FAR_SIZE * FAR_SIZE),
  farSurface = Buffer.alloc(FAR_SIZE * FAR_SIZE * 4),
  byte = (v) => Math.round(255 * Math.max(0, Math.min(1, v)));
for (let j = 0; j < FAR_SIZE; j++)
  for (let i = 0; i < FAR_SIZE; i++) {
    const x = -farHalf + i * FAR_SPACING,
      z = -farHalf + j * FAR_SPACING,
      k = j * FAR_SIZE + i,
      f = L.sample(x, z, { canopy: false });
    farHeights[k] = ground(x, z);
    farSurface[k * 4] = byte(f.weights[0]);
    farSurface[k * 4 + 1] = byte(f.weights[1]);
    farSurface[k * 4 + 2] = byte(f.weights[2]);
    farSurface[k * 4 + 3] = byte(world.vegetation.cover(f, f.elevation, f.slope));
  }
fs.writeFileSync(path.join(dir, 'far_height.f32'), Buffer.from(farHeights.buffer));
fs.writeFileSync(path.join(dir, 'far_surface.rgba8'), farSurface);
log(`far grid ${FAR_SIZE} x ${FAR_SIZE} at ${FAR_SPACING} m`);

// The surface textures the ground material samples (original, generated by bake/surfaces.py).
const surfaces = path.join(dir, 'surfaces');
fs.mkdirSync(surfaces, { recursive: true });
for (const f of ['albedo-ao.png', 'normal-roughness-height.png', 'manifest.json'])
  fs.copyFileSync(path.join(immersion, 'assets/surfaces', f), path.join(surfaces, f));

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
      ground(t.x, t.z),
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
  stones.push([r.id, r.x, r.z, ground(r.x, r.z), r.s, r.angle, r.family, r.stretch].join(','));
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
  'assets/surfaces/manifest.json',
].map((f) => ({ file: `immersion/${f}`, sha256: sha(path.join(immersion, f)) }));
const files = fs
  .readdirSync(dir, { recursive: true })
  .filter((f) => f !== 'manifest.json' && fs.statSync(path.join(dir, f)).isFile())
  .map((f) => ({ file: f, sha256: sha(path.join(dir, f)) }));
const manifest = {
  schema: 'terluna.game.world/1',
  world: { id: world.id, status: world.status },
  evidence:
    'The development cove: an authored placeholder landscape and planting rule for building ' +
    'the engine, not a proposed Open Moon environment.',
  coordinates:
    'Files are in world metres (x east-west, z north-south, y up). Unreal X = x, Y = z, ' +
    'Z = y, in centimetres. ' +
    (CURVED
      ? `Heights follow a sphere of radius ${C.worldRadius('moon')} m touching the origin ` +
        '(Unreal SkyAtmosphere planet top at the world origin); sea level is on that sphere.'
      : 'Heights are uncurved; sea level is 0 m.'),
  curvature: CURVED ? { radius_m: C.worldRadius('moon'), touching: 'world origin' } : null,
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
    layers_packed: {
      file: 'layers.rgba8',
      vertices: n,
      spacing_m: LAYER_STEP * SPACING,
      channels: LAYERS,
      note: 'uint8 weights, row j = z0 + j * spacing, the grid the layers are computed on',
    },
    surfaces: {
      folder: 'surfaces',
      layer_for_channel: { sand: 'sand', gravel: 'gravel', rock: 'basalt', soil: 'soil' },
      note: 'Atlas layout, tile sizes and encodings are in surfaces/manifest.json.',
    },
  },
  far: {
    vertices: FAR_SIZE,
    spacing_m: FAR_SPACING,
    origin_m: { x: -farHalf, z: -farHalf },
    heights: 'far_height.f32 (float32 little-endian metres, same curvature as the landscape)',
    surface:
      'far_surface.rgba8 (uint8: sand, gravel, rock weights, soil the rest; alpha the expected ' +
      "canopy cover of the planting rule, the web experience's distant woodland shading)",
    evidence: 'The development cove\'s placeholder height field and rules, continued past the landscape.',
  },
  vegetation: { trees: rows.length - 1, file: 'trees.csv', lattice_spacing_m: s },
  sources,
  files,
};
fs.writeFileSync(path.join(dir, 'manifest.json'), JSON.stringify(manifest, null, 1));
log(`wrote ${dir}`);
