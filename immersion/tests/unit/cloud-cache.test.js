/* Renderer-independent cache orchestration checks. Pixel checks are separate. */
import test from 'node:test';
import assert from 'node:assert/strict';
import * as T from 'three';
import { OM } from '../../engine/om.js';
import '../../engine/atmosphere.js';
import CloudRenderer from '../../engine/cloud-renderer.js';
import fs from 'node:fs';
import { loadColumns } from '../../engine/columns.js';
import cove from '../../world/cove.js';
// Engine systems read the active world; these tests run in the development cove.
OM.world = cove;
// Baked by `npm run bake:columns` (run automatically before `npm test`).
loadColumns(
  JSON.parse(fs.readFileSync(new URL('../../assets/columns/columns.json', import.meta.url))),
);
function fixture() {
  const a = new OM.Atmosphere(T, {
    worlds: {
      moon: { direct: [[10, 10, 10]] },
      earth: { direct: [[10, 10, 10]] },
      moon_no_ozone: { direct: [[10, 10, 10]] },
    },
  });
  const renderer = {
    toneMapping: T.ACESFilmicToneMapping,
    toneMappingExposure: 0.24,
    target: null,
    draws: [],
    getRenderTarget() {
      return this.target;
    },
    setRenderTarget(t) {
      this.target = t;
    },
    render() {
      this.draws.push(this.target?.scissor?.toArray() || []);
    },
  };
  a.columnClouds.initRenderer(renderer, '');
  return { a, c: a.columnClouds, r: renderer };
}
test('Cloud bake tiles cover the entire angular cache once', () => {
  const { c, r } = fixture();
  c.select('fair');
  c.update('moon');
  c.setQuality('economy');
  assert.equal(c.prepare(r), true);
  assert.equal(r.draws.length, 10);
  const tiles = r.draws.slice(0, -1);
  assert.equal(
    tiles.reduce((n, a) => n + a[2] * a[3], 0),
    768 * 384,
  );
  assert.equal(c.bakeCount, 1);
  assert.equal(r.target, null);
  assert.equal(r.toneMapping, T.ACESFilmicToneMapping);
  c.dispose();
});
test('Stable observer and sub-cadence wind reuse the existing transfer textures', () => {
  const { a, c, r } = fixture();
  c.select('convection');
  c.update('moon');
  c.prepare(r);
  const count = r.draws.length;
  a.uniforms.uTime.value = 1.5;
  a.uniforms.uObserver.value.set(10, 10);
  assert.equal(c.prepare(r), false);
  assert.equal(r.draws.length, count);
  c.dispose();
});
test('Wind cadence, observer cell, Sun direction and quality invalidate the cloud bake', () => {
  const { a, c, r } = fixture();
  c.select('fair');
  c.update('moon');
  c.prepare(r);
  a.uniforms.uTime.value = 2;
  assert.equal(c.prepare(r), true);
  a.uniforms.uObserver.value.x = 35;
  assert.equal(c.prepare(r), true);
  a.uniforms.uSun.value.set(0.7, 0.7, 0);
  assert.equal(c.prepare(r), true);
  c.setQuality('balanced');
  assert.equal(c.skyTarget.width, 1536);
  assert.equal(c.prepare(r), true);
  assert.equal(c.bakeCount, 5);
  c.dispose();
});
test('Reference selection clears diagnosed cloud rendering and cache flags', () => {
  const { a, c, r } = fixture();
  c.select('fog');
  c.update('moon');
  c.prepare(r);
  assert.ok(a.uniforms.uColumnCacheReady.value);
  c.select('reference');
  c.update('moon');
  assert.equal(a.uniforms.uColumnMode.value, 0);
  assert.equal(a.uniforms.uColumnCacheReady.value, 0);
  assert.equal(a.uniforms.uColumnShadowReady.value, 0);
  assert.equal(c.snapshot(), null);
  c.dispose();
});
test('No-ozone export reports the actual selected transport profile', () => {
  const { c } = fixture();
  c.select('fair');
  c.update('moon_no_ozone');
  assert.equal(c.snapshot().optics.ozoneEnabled, false);
  assert.equal(c.snapshot().optics.ozone, 0);
  c.dispose();
});
test('A failed render restores renderer state and leaves the bake invalid', () => {
  const { a, c, r } = fixture();
  c.select('fair');
  c.update('moon');
  r.render = () => {
    throw Error('synthetic draw failure');
  };
  assert.throws(() => c.prepare(r), /synthetic/);
  assert.equal(c.renderKey, null);
  assert.equal(a.uniforms.uColumnCacheReady.value, 0);
  assert.equal(r.target, null);
  assert.equal(r.toneMappingExposure, 0.24);
  c.dispose();
});

test('Periodic cloud noise lattice is deterministic, bounded and nonconstant', () => {
  const { noiseLattice } = CloudRenderer,
    a = noiseLattice(),
    b = noiseLattice();
  assert.equal(a.length, 64 ** 3);
  assert.deepEqual(a, b);
  assert.equal(new Set(a).size, 256);
  assert.throws(() => noiseLattice(128));
});
