import test from 'node:test';
import assert from 'node:assert/strict';
import { adapt, targetExposure, REFERENCE_EXPOSURE } from '../../engine/exposure.js';

test('Clear lunar noon keeps the calibrated reference exposure', () => {
  assert.ok(Math.abs(targetExposure(94200) - REFERENCE_EXPOSURE) < 1e-12);
  assert.ok(targetExposure(2) > targetExposure(2000));
  assert.equal(targetExposure(NaN), REFERENCE_EXPOSURE);
});

test('Adaptation depends on elapsed time, not on frame count', () => {
  const coarse = adapt(2350, 13.6, 0.5);
  let fine = 2350;
  for (let i = 0; i < 50; i++) fine = adapt(fine, 13.6, 0.01);
  assert.ok(Math.abs(Math.log(coarse / fine)) < 1e-9);
});

test('A scene that brightens is adapted to within a second', () => {
  // From midnight earthlight to bright twilight: the white-out case.
  const after = adapt(2350, 13.6, 1);
  assert.ok(after / 13.6 < 1.5, `still ${after.toFixed(1)} after one second`);
  // Darkening is slower, as for eyes.
  const dark = adapt(13.6, 2350, 1);
  assert.ok(dark / 13.6 < 20 && dark > 13.6);
});

test('Invalid current exposure snaps to the target', () => {
  assert.equal(adapt(NaN, 5, 0.016), 5);
  assert.equal(adapt(0, 5, 0.016), 5);
});
