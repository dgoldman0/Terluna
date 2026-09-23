import test from 'node:test';
import assert from 'node:assert/strict';
import fs from 'node:fs';
import path from 'node:path';
import { fileURLToPath } from 'node:url';
import C from '../../engine/core.js';
import { OM } from '../../engine/om.js';
import cove from '../../world/cove.js';
// Engine systems read the active world; these tests run in the development cove.
OM.world = cove;
const __dirname = path.dirname(fileURLToPath(import.meta.url));
const close = (a, b, tol = 1e-9) => assert.ok(Math.abs(a - b) < tol, `${a} != ${b}`);
test('Sun follows a continuous normalized full-cycle path', () => {
  for (let i = 0; i <= 1000; i++) {
    const s = C.sunAt(i / 1000);
    close(s.x * s.x + s.y * s.y + s.z * s.z, 1);
  }
});
test('Noon, sunset, midnight, sunrise have expected positions', () => {
  close(C.sunAt(0).elevation, 90);
  close(C.sunAt(0.25).elevation, 0);
  close(C.sunAt(0.5).elevation, -90);
  close(C.sunAt(0.75).elevation, 0);
});
test('Sun phase wraps both directions', () => {
  close(C.sunAt(-0.25).x, -1);
  close(C.sunAt(1.25).x, 1);
});
test('Solar elevation inverses preserve dawn and dusk branches', () => {
  for (let a = -90; a <= 90; a += 5) {
    for (const rising of [false, true])
      close(C.sunAt(C.phaseForElevation(a, rising)).elevation, a, 1e-6);
  }
});
test('The lunar clock has the expected period', () => close(C.PERIOD / C.DAY, 29.53059));
test('Calibrated vertical field of view uses actual viewing geometry', () =>
  close(C.calibratedFov(30, 65), (2 * Math.atan(15 / 65) * 180) / Math.PI));
test('Calibration rejects invalid input', () => assert.throws(() => C.calibratedFov(0, 50)));
test('Cloud optical depth has the correct dimensions and value', () =>
  close(C.cloudTau(0.1, 1000, 10), 15));
test('Cloud optical depth responds correctly to droplet radius', () =>
  close(C.cloudTau(0.1, 1000, 20), 7.5));
test('Cloud data rejects negative water content', () =>
  assert.throws(() => C.cloudTau(-1, 1000, 10)));
test('Weather episode is bounded and finite', () => {
  for (let t = -100; t <= 15000; t += 13) {
    const w = cove.weather.at(t);
    assert.ok(w.rain >= 0 && w.rain <= 10);
    assert.ok(w.coverage >= 0 && w.coverage <= 1);
    assert.ok(w.humidity > 0 && w.humidity <= 1);
    assert.ok(Object.values(w).every(Number.isFinite));
  }
});
test('Episode forcing is continuous at keyframe boundaries', () => {
  for (const a of cove.weather.EPISODE.slice(1, -1)) {
    const x = cove.weather.at(a[0] - 0.001),
      y = cove.weather.at(a[0] + 0.001);
    for (const k of Object.keys(x)) close(x[k], y[k], 1e-4);
  }
});
test('Clear forcing does not inherit a rain episode', () =>
  close(cove.weather.at(8000, 'clear').rain, 0));
test('Fog has prescribed low visibility and no rain', () => {
  const w = cove.weather.at(8000, 'fog');
  close(w.visibility, 110);
  close(w.rain, 0);
});
test('Reservoir zero-rate branch preserves inflow', () => {
  const r = C.reservoirStep(0.2, 0.01, 0, 0, 2, 10);
  close(r.storage, 0.3);
  close(r.residual, 0);
});
test('Reservoir overflow is explicitly recorded', () => {
  const r = C.reservoirStep(0.8, 0.1, 0, 0, 1, 10);
  close(r.storage, 1);
  close(r.overflow, 0.8);
  close(r.residual, 0);
});
test('Reservoir exponential drying matches its analytic solution', () => {
  const r = C.reservoirStep(1, 0, 0.01, 0.02, 2, 10);
  close(r.storage, Math.exp(-0.3));
  close(r.residual, 0);
});
test('Capped reservoir conserves water across boundary contact', () => {
  const r = C.reservoirStep(0.1, 0.3, 0.01, 0.02, 0.5, 15);
  close(r.storage, 0.5);
  assert.ok(r.overflow > 0);
  close(r.residual, 0);
});
test('Boundary-start reservoir can drain below capacity', () => {
  const r = C.reservoirStep(1, 0, 0.01, 0.02, 1, 10);
  assert.ok(r.storage < 1);
  close(r.residual, 0);
});
test('Invalid water data fails immediately', () =>
  assert.throws(() => C.reservoirStep(3, 1, 0, 0, 1, 10)));
test('Reservoir subdivision is exact for constant forcing', () => {
  let s = 0.2;
  for (let i = 0; i < 100; i++) s = C.reservoirStep(s, 0.014, 0.001, 0.002, 1, 2).storage;
  close(s, C.reservoirStep(0.2, 0.014, 0.001, 0.002, 1, 200).storage, 1e-12);
});
test('Water ledger conserves exposure plus canopy column water', () => {
  for (let t = 0; t <= 14400; t += 300) {
    const l = cove.weather.ledgerAt(t);
    assert.ok(Math.abs(l.residual) < 1e-9);
    assert.ok(l.exposed >= 0 && l.exposed <= 2);
    assert.ok(l.canopy >= 0 && l.canopy <= 1.4);
    assert.ok(l.leaf >= 0 && l.leaf <= 0.35);
  }
});
test('Sheltered ground remains dry through the rain', () =>
  close(cove.weather.ledgerAt(14400).sheltered, 0));
test('Clearing sky leaves persistent surface water', () => {
  const l = cove.weather.ledgerAt(14400);
  assert.equal(cove.weather.at(14400).rain, 0);
  assert.ok(l.exposed > 0.05 && l.canopy > 0.05);
});
test('Replay produces identical weather-water state', () =>
  assert.deepEqual(cove.weather.ledgerAt(8210), cove.weather.ledgerAt(8210)));
test('Two-second and ten-second water steps are close', () => {
  const a = cove.weather.ledgerAt(9800, 'episode', 2),
    b = cove.weather.ledgerAt(9800, 'episode', 10);
  close(a.exposed, b.exposed, 0.001);
  close(a.canopy, b.canopy, 0.001);
  close(a.leaf, b.leaf, 0.001);
});
test('Spherical-drop drag predicts a lower lunar terminal speed', () => {
  const a = C.dropletTerminalSpeed(0.0007, 1.62, 1.45),
    b = C.dropletTerminalSpeed(0.0007, 9.80665, 1.2);
  assert.ok(a > 0 && b > a && b < 10);
});
test('Wave frequency follows square-root gravity scaling', () =>
  close(C.waveOmega(0.2, 12, 1.62) / C.waveOmega(0.2, 12, 9.80665), Math.sqrt(1.62 / 9.80665)));
test('Wave gradients agree with centered finite differences', () => {
  const x = 3,
    z = 7,
    t = 24,
    w = C.waveAt(x, z, t),
    e = 1e-4;
  close(w.nx, -(C.waveAt(x + e, z, t).y - C.waveAt(x - e, z, t).y) / (2 * e), 1e-7);
  close(w.nz, -(C.waveAt(x, z + e, t).y - C.waveAt(x, z - e, t).y) / (2 * e), 1e-7);
});
test('Terrain is finite around every waypoint', () => {
  for (const p of cove.waypoints) {
    const y = C.groundHeight(p.x, p.z);
    assert.ok(Number.isFinite(y) && y > 0);
  }
});
test('Roof footprint separates dry inside from rain outside', () => {
  assert.equal(cove.roofMask(20, 33), 1);
  assert.equal(cove.roofMask(20, 20), 0);
});
test('Path centerline connects the requested places', () => {
  for (const p of cove.waypoints) assert.ok(cove.pathDistance(p.x, p.z) < 5);
});
test('Procedural seeds are deterministic', () => {
  const a = C.rng(123),
    b = C.rng(123);
  for (let i = 0; i < 100; i++) close(a(), b());
});
test('The atlas stores all three inherited profiles and the full angular range', () => {
  const d = JSON.parse(fs.readFileSync(path.join(__dirname, '../../assets/sky/atmosphere.json')));
  assert.equal(Object.keys(d.worlds).length, 3);
  for (const w of Object.values(d.worlds)) {
    assert.equal(w.suns[0], -90);
    assert.equal(w.suns.at(-1), 90);
    assert.equal(w.width, 33);
    assert.equal(w.height, 41);
    assert.ok(w.quantization_error_relative_to_frame_max < 0.00025);
  }
});
test('Cloud displacement integrates a constant wind', () =>
  close(cove.weather.windDistance(3600, 'clear'), cove.weather.at(0, 'clear').wind * 3600, 1e-8));
test('Cloud displacement integrates varying wind reproducibly', () => {
  close(cove.weather.windDistance(8100), cove.weather.windDistance(8100), 1e-10);
  close(
    cove.weather.windDistance(8100, 'episode', 10),
    cove.weather.windDistance(8100, 'episode', 2),
    0.02,
  );
  assert.ok(cove.weather.windDistance(8200) > cove.weather.windDistance(8100));
});
