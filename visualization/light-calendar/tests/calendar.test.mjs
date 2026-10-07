import test from "node:test";
import assert from "node:assert/strict";
import { readFileSync } from "node:fs";
import { calculate, findCrossings } from "../dist/worker.mjs";
import { stateAt } from "../dist/model/evaluator.mjs";

const read = (name) =>
  JSON.parse(
    readFileSync(new URL(`../dist/data/${name}.json`, import.meta.url)),
  );
const astronomy = read("astronomy"),
  transfer = read("transfer");
const settings = {
  year: 2030,
  month: 2,
  longitude: -73.875,
  latitude: 28.875,
  earth: true,
  scale: "TT",
};

test("a grazing dip that fits between samples produces both crossings", () => {
  const times = [0, 1800000, 3600000, 5400000];
  const fn = (t) => ((t - 2100000) / 60000) ** 2 - 1;
  const roots = findCrossings(times, times.map(fn), fn);
  assert.equal(roots.length, 2);
  assert.ok(Math.abs(roots[0] - 2040000) < 500);
  assert.ok(Math.abs(roots[1] - 2160000) < 500);
});

test("calendar events straddle the displayed geometric or brightness threshold", () => {
  const result = calculate(settings, astronomy, transfer);
  assert.equal(result.daily.length, 31);
  assert.ok(result.events.length > 4);
  for (const event of result.events) {
    const state = (delta) =>
      stateAt(
        event.ms + delta,
        settings.longitude,
        settings.latitude,
        astronomy,
        transfer,
        true,
        "TT",
      );
    const value = (s) =>
      event.key === "light"
        ? s.total - event.level
        : s[`${event.key}_elevation_deg`] + s[`${event.key}_radius_deg`];
    assert.ok(value(state(-1000)) * value(state(1000)) <= 0, event.label);
  }
  for (const day of result.daily) {
    assert.ok(day.sunHours >= 0 && day.sunHours <= 24);
    assert.ok(day.min >= 0 && day.max >= day.min);
  }
});

test("Earthlight toggle preserves the celestial geometry and changes only the light contribution", () => {
  const on = calculate(
    { ...settings, year: 2038, month: 1 },
    astronomy,
    transfer,
  );
  const off = calculate(
    { ...settings, year: 2038, month: 1, earth: false },
    astronomy,
    transfer,
  );
  assert.ok(on.series.some((r) => r.earth > 0.01));
  for (let i = 0; i < on.series.length; i++) {
    assert.equal(off.series[i].earth, 0);
    assert.equal(off.series[i].total, off.series[i].solar);
    assert.equal(
      off.series[i].sun_elevation_deg,
      on.series[i].sun_elevation_deg,
    );
    assert.equal(
      off.series[i].earth_elevation_deg,
      on.series[i].earth_elevation_deg,
    );
    assert.ok(on.series[i].total >= off.series[i].total);
  }
});

test("leap-day and end-of-range calendars stay within the supported dates", () => {
  assert.equal(
    calculate(
      { ...settings, year: 2000, month: 1, earth: false },
      astronomy,
      transfer,
    ).daily.length,
    29,
  );
  const final = calculate(
    { ...settings, year: 2500, month: 11, earth: false },
    astronomy,
    transfer,
  );
  assert.equal(final.daily.length, 31);
  assert.ok(final.series.every((r) => r.ms < Date.UTC(2501, 0, 1)));
  assert.throws(
    () => calculate({ ...settings, year: 2501 }, astronomy, transfer),
    /2000 and 2500/,
  );
});
