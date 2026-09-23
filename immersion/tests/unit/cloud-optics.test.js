// The viewer's real-time RGB cloud-lighting approximation, checked against analytic limits.
import test from 'node:test';
import assert from 'node:assert/strict';
import O from '../../engine/cloud-optics.js';
import C from '../../engine/core.js';
import fs from 'node:fs';
import { loadColumns, getColumn } from '../../engine/columns.js';
import { OM } from '../../engine/om.js';
import cove from '../../world/cove.js';
// Engine systems read the active world; these tests run in the development cove.
OM.world = cove;
// Baked by `npm run bake:columns` (run automatically before `npm test`).
loadColumns(
  JSON.parse(fs.readFileSync(new URL('../../assets/columns/columns.json', import.meta.url))),
);
const near = (a, b, e = 1e-7) => assert.ok(Math.abs(a - b) < e, `${a} vs ${b} exceeds ${e}`);
test('Scene weather scenarios leave the existing rainfall ledger unchanged', () => {
  const before = cove.weather.ledgerAt(14400);
  getColumn('moon', 'convection');
  getColumn('moon', 'fog');
  assert.deepEqual(cove.weather.ledgerAt(14400), before);
});
test('Baked columns interpolate their rows and alias the no-ozone Moon', () => {
  const fog = getColumn('moon', 'fog'),
    alias = getColumn('moon_no_ozone', 'fog');
  assert.equal(alias.summary.cloudTop_m, fog.summary.cloudTop_m);
  const [a, b] = fog.rows;
  near(fog.sample((a.z + b.z) / 2).T, (a.T + b.T) / 2, 1e-9);
  assert.equal(fog.texture.data.length, fog.texture.n * 4);
});
test('Optical vertical molecular column converges to its analytic exponential integral', () => {
  for (const w of ['earth', 'moon']) {
    const p = O.profile(w),
      a = O.columns(w, 1000, 1, 1024),
      exact = p.density * p.H * (Math.exp(-1000 / p.H) - Math.exp(-p.top / p.H));
    assert.ok(Math.abs(a.molecular / exact - 1) < 0.00001);
  }
});
test('No-ozone atmosphere removes the selected ozone attenuation only', () => {
  const a = O.transmission('moon', 10000, 1),
    b = O.transmission('moon_no_ozone', 10000, 1);
  assert.ok(b.every((v, i) => v >= a[i]));
  near(O.columns('moon_no_ozone', 10000, 1).ozone, 0);
});
test('Sunlight attenuation decreases with height at zenith', () => {
  const a = O.transmission('moon', 0, 1),
    b = O.transmission('moon', 50000, 1);
  assert.ok(b.every((v, i) => v > a[i] && v <= 1));
});
test('Planet occultation distinguishes surface night from an illuminated upper cloud', () => {
  assert.deepEqual(O.transmission('moon', 0, -0.1), [0, 0, 0]);
  assert.ok(O.transmission('moon', 60000, -0.1).some((v) => v > 0));
  assert.deepEqual(O.transmission('earth', 10000, -0.1), [0, 0, 0]);
});
test('RGB ray quadrature converges across zenith and low-Sun cloud cases', () => {
  for (const [h, mu] of [
    [0, 1],
    [10000, 0.3],
    [50000, 0.05],
    [60000, -0.1],
  ]) {
    const a = O.transmission('moon', h, mu, 96),
      b = O.transmission('moon', h, mu, 384);
    for (let i = 0; i < 3; i++) near(a[i], b[i], 0.003);
  }
});
test('Light lookup keeps the exact CPU samples and finite nonnegative transmission', () => {
  const a = O.buildTable('moon', 65000);
  assert.ok(a.data.every((x) => Number.isFinite(x) && x >= 0 && x <= 1));
  const i = 105,
    j = 40,
    v = (2 * i) / (a.width - 1) - 1,
    t = O.transmission('moon', (65000 * j) / (a.height - 1), Math.sign(v) * v * v);
  t.forEach((v, k) => near(a.data[(j * a.width + i) * 4 + k], v, 1e-7));
});
test('Longer eye paths and the lunar column attenuate cloud contrast more strongly', () => {
  const m = O.segmentColumn('moon', 2, 0.2, 50000),
    e = O.segmentColumn('earth', 2, 0.2, 50000);
  assert.ok(m > e);
  assert.ok(O.segmentColumn('moon', 2, 0.2, 100000) > m);
});
test('Invalid optical rays are rejected', () => {
  for (const args of [
    ['moon', -1, 1],
    ['moon', 0, 2],
    ['x', 0, 1],
    ['moon', 0, NaN],
  ])
    assert.throws(() => O.transmission(...args));
});
