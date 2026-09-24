/* The Unreal game's data exports carry the world and sky unchanged. */
import { test } from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import os from 'node:os';
import path from 'node:path';
import zlib from 'node:zlib';
import { execFileSync } from 'node:child_process';
import { fileURLToPath } from 'node:url';
import { OM } from '../../engine/om.js';
import C from '../../engine/core.js';
import cove from '../../world/cove.js';
import { treeAt } from '../../world/cove-flora.js';

OM.world = cove;
const immersion = path.resolve(path.dirname(fileURLToPath(import.meta.url)), '../..'),
  out = fs.mkdtempSync(path.join(os.tmpdir(), 'game-export-'));

test('the world export keeps heights, the planting rule and the landscape transform', () => {
  execFileSync(
    process.execPath,
    [
      'bake/game/world.mjs',
      ...['--out', out, '--size', '65', '--spacing', '8', '--far-size', '33', '--far-spacing', '64'],
    ],
    {
      cwd: immersion,
      stdio: 'ignore',
    },
  );
  const dir = path.join(out, 'world'),
    m = JSON.parse(fs.readFileSync(path.join(dir, 'manifest.json'))),
    L = m.landscape,
    r16 = fs.readFileSync(path.join(dir, 'height.r16')),
    [locX, locY, locZ] = L.unreal_transform.location_cm,
    scaleZ = L.unreal_transform.scale[2];
  assert.equal(r16.length, 65 * 65 * 2);
  assert.equal(locX, L.origin_m.x * 100);
  assert.equal(locY, L.origin_m.z * 100);
  for (const [i, j] of [
    [0, 0],
    [32, 32],
    [64, 10],
    [7, 51],
  ]) {
    const v = r16.readUInt16LE((j * 65 + i) * 2),
      metres = (locZ + ((v - 32768) * scaleZ) / 128) / 100,
      want = C.groundHeight(L.origin_m.x + i * 8, L.origin_m.z + j * 8, 'moon');
    assert.ok(Math.abs(metres - want) <= L.height_quantum_m, `${i},${j}: ${metres} vs ${want}`);
  }
  for (const layer of L.layers) assert.equal(fs.statSync(path.join(dir, layer.file)).size, 65 * 65);
  const packed = L.layers_packed;
  assert.equal(packed.vertices, 33);
  assert.equal(fs.statSync(path.join(dir, packed.file)).size, 33 * 33 * 4);
  // The packed layers are the full-resolution layers at every LAYER_STEP-th sample.
  const sand = fs.readFileSync(path.join(dir, 'layer_sand.r8')),
    rgba = fs.readFileSync(path.join(dir, packed.file));
  for (const [i, j] of [
    [0, 0],
    [16, 16],
    [32, 5],
  ])
    assert.equal(rgba[(j * 33 + i) * 4], sand[j * 2 * 65 + i * 2]);
  const far = m.far,
    farHeights = new Float32Array(new Uint8Array(fs.readFileSync(path.join(dir, 'far_height.f32'))).buffer);
  assert.equal(farHeights.length, 33 * 33);
  assert.equal(fs.statSync(path.join(dir, 'far_surface.rgba8')).size, 33 * 33 * 4);
  for (const [i, j] of [
    [0, 0],
    [16, 16],
    [30, 3],
  ]) {
    const want = C.groundHeight(far.origin_m.x + i * 64, far.origin_m.z + j * 64, 'moon');
    assert.ok(Math.abs(farHeights[j * 33 + i] - want) < 1e-3 * Math.max(1, Math.abs(want)), `far ${i},${j}`);
  }
  for (const f of ['albedo-ao.png', 'normal-roughness-height.png', 'manifest.json'])
    assert.ok(fs.statSync(path.join(dir, 'surfaces', f)).size > 0, f);
  const rows = fs.readFileSync(path.join(dir, 'trees.csv'), 'utf8').trim().split('\n').slice(1);
  assert.equal(rows.length, m.vegetation.trees);
  for (const row of rows.filter((r) => r.endsWith(',0')).slice(0, 50)) {
    const [id, x, z] = row.split(','),
      [, i, j] = id.split(':').map(Number),
      t = treeAt(i, j);
    assert.ok(t, id);
    assert.ok(Math.abs(t.x - +x) < 1e-3 && Math.abs(t.z - +z) < 1e-3, id);
  }
  assert.ok(
    rows.some((r) => r.endsWith(',1')),
    'the detailed plan is included',
  );
});

test('the sky export holds the clear-sky atlas frame for frame', () => {
  execFileSync('python3', ['bake/game/sky.py', '--out', out], { cwd: immersion, stdio: 'ignore' });
  const dir = path.join(out, 'sky'),
    layout = JSON.parse(fs.readFileSync(path.join(dir, 'atlas.json'))),
    source = JSON.parse(fs.readFileSync(path.join(immersion, 'assets/sky/atmosphere.json'))),
    w = source.worlds.moon,
    raw = zlib.gunzipSync(Buffer.from(w.data, 'base64')),
    half = (v) => {
      const s = v & 0x8000 ? -1 : 1,
        e = (v >> 10) & 31,
        f = v & 1023;
      return s * (e ? 2 ** (e - 15) * (1 + f / 1024) : (2 ** -14 * f) / 1024);
    };
  const L = layout.worlds.moon,
    [width] = L.size_px,
    [tw, th] = L.tile_px,
    floats = new Float32Array(
      new Uint8Array(fs.readFileSync(path.join(dir, 'atlas_moon.f32'))).buffer,
    );
  for (const [f, b, a, c] of [
    [173, 40, 0, 0],
    [173, 20, 16, 1],
    [100, 25, 3, 2],
    [0, 30, 31, 0],
  ]) {
    const want = half(raw.readUInt16LE(((f * th + b) * tw + a) * 3 * 2 + c * 2)) * w.scale[f],
      row = Math.floor(f / 14) * th + (th - 1 - b),
      col = (f % 14) * tw + a;
    assert.ok(
      Math.abs(floats[(row * width + col) * 3 + c] - want) <= 1e-6 * Math.abs(want) + 1e-12,
    );
  }
  const engine = JSON.parse(fs.readFileSync(path.join(dir, 'engine_atmosphere.json')));
  assert.equal(engine.schema, 'terluna.illumination.engine-atmosphere/1');
  const manifest = JSON.parse(fs.readFileSync(path.join(dir, 'manifest.json')));
  assert.ok(manifest.sources.some((s) => s.file === 'immersion/assets/sky/atmosphere.json'));
});
