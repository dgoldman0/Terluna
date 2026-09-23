import { test } from 'node:test';
import assert from 'node:assert/strict';
global.OM = {};
import L from '../../world/landscape.js';
import C from '../../engine/core.js';
import E from '../../engine/vegetation.js';
import W from '../../engine/surface-water.js';
import { PondModel, triangleVolume, clipTriangle } from '../../engine/ponds.js';
import { OM } from '../../engine/om.js';
import cove from '../../world/cove.js';
// Engine systems read the active world; these tests run in the development cove.
OM.world = cove;
const close = (a, b, e = 1e-8) => assert.ok(Math.abs(a - b) < e, `${a} != ${b}`);
let w, p;
function setup() {
  if (!w) {
    cove.flora();
    w = new W();
    p = new PondModel(w);
  }
  return { w, p };
}
test('Triangle hypsometry is exact for flat, one-wet-corner and fully submerged cases', () => {
  close(triangleVolume(0, [1, 2, 3], 0.5), 0);
  close(triangleVolume(2, [1, 1, 1], 0.5), 0.5);
  close(triangleVolume(1, [0, 2, 2], 0.5), 1 / 24);
  close(triangleVolume(3, [0, 1, 2], 0.5), 1);
});
test('Clipped triangle integration equals the analytic partial volume', () => {
  const points = [
      [0, 0, 0],
      [0, 1, 1],
      [1, 0, 3],
    ],
    sorted = [0, 1, 3];
  for (const h of [0.05, 0.7, 1, 1.1, 2.8, 3.5]) {
    const poly = clipTriangle(points, h);
    let volume = 0;
    for (let i = 1; i < poly.length - 1; i++) {
      const a = poly[0],
        b = poly[i],
        c = poly[i + 1],
        area = Math.abs((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])) / 2;
      volume += (area * (a[2] + b[2] + c[2])) / 3;
    }
    close(volume, triangleVolume(h, sorted, 0.5));
  }
});
test('A dry state produces no pond geometry or invented water', () => {
  const { w, p } = setup();
  assert.equal(p.summary.columnPondVolume_m3, 0);
  assert.equal(p.meshData().positions.length, 0);
});
test('Standing water, mobile storage and above-spill excess close the reconstruction budget', () => {
  const { w, p } = setup();
  w.seek(8100, 'episode');
  const before = JSON.stringify(w.ledger);
  p.update();
  assert.equal(JSON.stringify(w.ledger), before);
  const s = p.summary;
  assert.ok(s.activeBasins > 0);
  assert.ok(s.reconstructedVolume_m3 > 0);
  assert.ok(Math.abs(s.accountingResidual_m3) < 1e-7);
  assert.ok(s.maxBasinSolveError_m3 < 1e-7);
});
test('Rendered pond triangles integrate to exactly the represented volume', () => {
  const { p } = setup(),
    d = p.meshData();
  let v = 0;
  for (let i = 0; i < d.depths.length; i += 3) {
    const a = i * 3,
      b = a + 3,
      c = a + 6;
    const area =
      Math.abs(
        (d.positions[b] - d.positions[a]) * (d.positions[c + 2] - d.positions[a + 2]) -
          (d.positions[b + 2] - d.positions[a + 2]) * (d.positions[c] - d.positions[a]),
      ) / 2;
    v += (area * (d.depths[i] + d.depths[i + 1] + d.depths[i + 2])) / 3;
  }
  close(v, p.summary.reconstructedVolume_m3, 2e-5);
});
test('Connected basin surfaces share one elevation and stay below their spill level', () => {
  const { p } = setup();
  for (const g of p.groups) {
    if (g.volume > 0) {
      assert.ok(g.level <= g.spill + 1e-10);
      assert.ok(g.level >= g.min);
      for (const t of g.triangles)
        for (const point of clipTriangle(t.points, g.level)) assert.ok(point[2] >= -1e-10);
    }
  }
});
test('World curvature changes render positions while pond volumes remain identical', () => {
  const { p } = setup(),
    a = p.meshData('moon'),
    b = p.meshData('earth');
  assert.deepEqual(a.depths, b.depths);
  for (let i = 0; i < a.positions.length; i += 111) {
    const x = a.positions[i],
      z = a.positions[i + 2];
    close(
      b.positions[i + 1] - a.positions[i + 1],
      C.curvatureSag(x, z, 'moon') - C.curvatureSag(x, z, 'earth'),
      2e-6,
    );
  }
});
test('Reset clears pond geometry and restores the untouched antecedent soil budget', () => {
  const { w, p } = setup();
  w.seek(0, 'clear');
  p.update();
  assert.equal(p.meshData().positions.length, 0);
  assert.equal(p.summary.columnPondVolume_m3, 0);
  assert.ok(w.ledger.relativeResidual < 1e-12);
});
test('Stationary water state skips pool solving', () => {
  const { p } = setup();
  assert.equal(p.update(), false);
});
test('Coastal forcing has a nonnegative bounded bathymetric envelope', () => {
  for (let z = -60; z < 12; z += 1)
    for (let x = -12; x < 12; x += 2) {
      const e = C.coastalEnvelope(x, z, 6);
      assert.ok(e >= 0 && e <= 1);
    }
});
test('Land vertices carry zero wave amplitude within the local field', () => {
  for (const p of [
    [0, 9],
    [20, 33],
    [-56, 47],
  ])
    assert.equal(C.coastalEnvelope(...p), 0);
});
test('Envelope bounds the sum of wave amplitudes below 45% of sampled depth', () => {
  const data = L.initialize();
  for (let z = -40; z < -1; z += 0.7) {
    const h = Math.max(0, -L.gridSample(data.elevation, 0, z)),
      a = C.WATER_BANDS.reduce(
        (s, k) => s + 0.095 * Math.pow(0.095 / k, 1.12) * (0.4 + 6 * 0.18),
        0,
      );
    assert.ok(a * C.coastalEnvelope(0, z, 6) <= 0.45 * h + 1e-10);
  }
});
