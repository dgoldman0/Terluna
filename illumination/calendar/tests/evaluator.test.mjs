import test from "node:test";
import assert from "node:assert/strict";
import { readFileSync } from "node:fs";
import {
  geometry,
  ttOffset,
  julianDate,
  diskIntegral,
  crossing,
} from "../evaluator.mjs";

const astronomy = JSON.parse(
  readFileSync(new URL("../results/astronomy.json", import.meta.url)),
);
test("browser geometry agrees with the independent Python export across the supported dates", () => {
  for (const c of astronomy.golden_geometry) {
    const actual = geometry(c.jd_tt, c.longitude, c.latitude, astronomy);
    for (const [key, expected] of Object.entries(c.values)) {
      const delta = key.endsWith("azimuth_deg")
        ? Math.abs(((actual[key] - expected + 540) % 360) - 180)
        : Math.abs(actual[key] - expected);
      assert.ok(
        delta <= Math.max(1e-8, Math.abs(expected) * 2e-12),
        `${key}: ${actual[key]} vs ${expected}`,
      );
    }
  }
});
test("historical leap seconds and explicit TT do not depend on computer time zone", () => {
  assert.equal(ttOffset(Date.parse("2000-01-01T00:00:00Z")), 64.184);
  assert.equal(ttOffset(Date.parse("2016-12-31T23:59:59Z")), 68.184);
  assert.equal(ttOffset(Date.parse("2017-01-01T00:00:00Z")), 69.184);
  assert.equal(julianDate(Date.parse("2000-01-01T12:00:00Z"), "TT"), 2451545);
  assert.ok(
    Math.abs(
      julianDate(Date.parse("2000-01-01T12:00:00Z")) - 2451545 - 64.184 / 86400,
    ) < 1e-9,
  );
});
test("phase weights preserve total energy and rotate the bright limb at the horizon", () => {
  const constant = () => [1, 2];
  for (const phase of [0, 60, 90, 144, 179]) {
    const value = diskIntegral(constant, 0, 1, phase);
    assert.ok(Math.abs(value[0] - 1) < 1e-12);
    assert.ok(Math.abs(value[1] - 2) < 1e-12);
  }
  const beam = (e) => [Math.max(0, Math.sin((e * Math.PI) / 180)), 0];
  assert.ok(diskIntegral(beam, 0, 1, 90, 0)[0] > 0.008);
  assert.equal(diskIntegral(beam, 0, 1, 90, Math.PI)[0], 0);
  assert.deepEqual(diskIntegral(constant, 0, 1, 180), [0, 0]);
});
test("crossing refinement respects the bracket and recognizes missing crossings", () => {
  assert.ok(Math.abs(crossing((ms) => ms / 1000, 0, 10000, 3.2, 1) - 3200) < 1);
  assert.equal(
    crossing((ms) => ms, 0, 10000, -1),
    null,
  );
});
