import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import { SkyBodies } from '../../engine/sky-bodies.js';
import C from '../../engine/core.js';
import { OM } from '../../engine/om.js';
import cove from '../../world/cove.js';
// Engine systems read the active world; these tests run in the development cove.
OM.world = cove;

// Baked by `npm run bake:sky` from illumination/ephemeris.py.
const product = JSON.parse(
  fs.readFileSync(new URL('../../assets/sky/site-sky.json', import.meta.url)),
);
const sky = new SkyBodies(product);
const close = (a, b, tol, what) =>
  a.forEach((v, i) => assert.ok(Math.abs(v - b[i]) <= tol, `${what}[${i}]: ${v} vs ${b[i]}`));

test('The experience evaluates the published ephemeris exactly at every golden sample', () => {
  for (const g of product.samples) {
    const s = sky.state(g.t_s);
    close(s.sun, g.sun_enu, 1e-9, 'sun');
    close(s.earth, g.earth_enu, 1e-9, 'earth');
    close(s.celestialToEnu.flat(), g.celestial_to_enu.flat(), 1e-9, 'celestial');
    assert.ok(Math.abs(s.earthFraction - g.earth_illuminated_fraction) < 1e-12);
    assert.ok(
      Math.abs(s.earthlightRatio / g.earthlight_ratio - 1) < 1e-9 || g.earthlight_ratio < 1e-12,
    );
    assert.ok(Math.abs(s.earthRotation - g.earth_rotation_rad) < 1e-9);
  }
});

test("The ephemeris Sun matches the scene's solar clock", () => {
  for (const phase of [0, 0.1, 0.25, 0.4, 0.5, 0.77, 0.99]) {
    const a = sky.scene(phase * sky.period).sun,
      b = C.sunAt(phase);
    close(a, [b.x, b.y, b.z], 1e-9, 'sun at phase ' + phase);
  }
});

test('Scene-frame celestial matrix is a proper rotation', () => {
  const m = sky.scene(12345678).celestialToScene;
  for (let i = 0; i < 3; i++)
    for (let j = 0; j < 3; j++) {
      const d = m[i][0] * m[j][0] + m[i][1] * m[j][1] + m[i][2] * m[j][2];
      assert.ok(Math.abs(d - (i === j ? 1 : 0)) < 1e-12);
    }
  const det =
    m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]) -
    m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]) +
    m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
  assert.ok(Math.abs(det - 1) < 1e-12);
});
