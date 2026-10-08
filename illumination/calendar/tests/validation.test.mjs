import test from "node:test";
import assert from "node:assert/strict";
import { readFileSync, existsSync } from "node:fs";
import { geometry, lightAt } from "../evaluator.mjs";

const path = (name) => new URL(`../results/${name}.json`, import.meta.url);
const available = existsSync(path("transfer")) && existsSync(path("validation"));

test("JavaScript finite-source light matches independent Python golden cases", { skip: !available }, () => {
  const transfer = JSON.parse(readFileSync(path("transfer")));
  const validation = JSON.parse(readFileSync(path("validation")));
  const astronomy = JSON.parse(readFileSync(path("astronomy")));
  assert.equal(validation.schema, "terluna.illumination.calendar-validation/1");
  assert.ok(validation.golden_light.length >= 40);
  const compare = (actual, expected, name, relativeTolerance = 2e-10) => {
    for (const key of ["solar_direct", "solar_diffuse", "solar", "earth_direct", "earth_diffuse", "earth", "total"]) {
      const error = Math.abs(actual[key] - expected[key]);
      assert.ok(error <= Math.max(1e-9, Math.abs(expected[key]) * relativeTolerance), `${name} ${key}: ${actual[key]} vs ${expected[key]}`);
    }
    assert.equal(actual.eclipse, expected.eclipse, name);
    assert.equal(actual.earthIncluded, expected.earthIncluded, name);
  };
  for (const c of validation.golden_light) {
    compare(lightAt(c.geometry, transfer, c.includeEarth, c.order), c.expected, c.name);
    if (c.kind === "native_ephemeris") {
      const g = geometry(c.jd_tt, c.longitude, c.latitude, astronomy);
      // Independent polynomial evaluation can differ by a few nanodegrees
      // centuries from J2000; this tolerance remains far below model accuracy.
      compare(lightAt(g, transfer, c.includeEarth, c.order), c.expected, `${c.name} end-to-end`, 2e-8);
    }
  }
});
