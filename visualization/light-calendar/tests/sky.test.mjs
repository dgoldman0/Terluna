import test from "node:test";
import assert from "node:assert/strict";
import { readFileSync, existsSync } from "node:fs";
import { stateAt } from "../dist/model/evaluator.mjs";
import { skyFrame, subObserverPoint, apply, localVector, lonLat } from "../dist/sky-frame.mjs";
import { span, setProducts } from "../dist/explorer-worker.mjs";
import { sentence, brightness, conditionOf, fromNow } from "../dist/words.mjs";

const read = (path) => JSON.parse(readFileSync(new URL(path, import.meta.url)));
const astronomy = read("../dist/data/astronomy.json"),
  transfer = read("../dist/data/transfer.json");
const horizons = read("../../sea-appearance/nubium-earth-inputs.json").horizons;
const solvedSky = read("../../../illumination/sky/results/solved_sky.json");
const control = {
  suns: solvedSky.worlds.earth_control.samples.map((s) => s.sun_deg),
  horizontal_lux: solvedSky.worlds.earth_control.samples.map((s) => s.total_horizontal_lux),
};
const angle = (a, b) => Math.abs(((a - b + 540) % 360) - 180);

test("Earth's face and its sunlit point match the saved JPL Horizons record", () => {
  const ms = (horizons.jd_tt - 2440587.5) * 86400000;
  const [lon, lat] = horizons.observer_moon_lon_lat_height_km;
  const s = stateAt(ms, lon, lat, astronomy, transfer, true, "TT");
  const frame = skyFrame(ms, "TT", lon, lat, astronomy);
  const [subLon, subLat] = subObserverPoint(frame, s.earth_elevation_deg, s.earth_azimuth_deg);
  assert.ok(angle(subLon, horizons.subobserver_lon_east_deg) < 0.2, `sub-observer longitude ${subLon}`);
  assert.ok(Math.abs(subLat - horizons.subobserver_lat_deg) < 0.2, `sub-observer latitude ${subLat}`);
  const sun = apply(frame.localToEarth, localVector(s.sun_elevation_deg, s.sun_azimuth_deg));
  const [solarLon, solarLat] = lonLat(sun);
  assert.ok(angle(solarLon, horizons.subsolar_lon_east_deg) < 0.2, `sub-solar longitude ${solarLon}`);
  assert.ok(Math.abs(solarLat - horizons.subsolar_lat_deg) < 0.2, `sub-solar latitude ${solarLat}`);
});

test("Earth crosses the Sun over the Sea of Showers at 2026's two lunar eclipses", () => {
  setProducts(astronomy, transfer);
  const result = span(
    { longitude: -17.875, latitude: 15.875, earth: true, scale: "UTC" },
    Date.UTC(2026, 0, 1),
    Date.UTC(2027, 0, 1),
    2 * 3600000,
    4,
  );
  const days = result.events
    .filter((e) => e.type === "eclipse-begins")
    .map((e) => new Date(e.ms).toISOString().slice(0, 10));
  assert.deepEqual(days, ["2026-03-03", "2026-08-28"]);
  const rises = result.events.filter((e) => e.type === "sunrise").length;
  assert.ok(rises >= 12 && rises <= 13, `${rises} sunrises in a year`);
});

test("the sentence follows evening and morning twilight", () => {
  const now = Date.UTC(2026, 9, 7);
  const s = { sun_elevation_deg: -20, sun_azimuth_deg: 270, earth_elevation_deg: 60, earth_azimuth_deg: 120 };
  const evening = sentence({ state: s, condition: 1, rising: false, now, lastSunset: now - 2 * 86400000, nextSunrise: now + 12 * 86400000 });
  assert.equal(evening, "The Sun set 2 days ago, but the deep air still glows with its light. Sunrise is 12 days away.");
  const morning = sentence({ state: s, condition: 1, rising: true, now, lastSunset: now - 14 * 86400000, nextSunrise: now + 6 * 3600000 });
  assert.match(morning, /^Morning twilight\..*rises in 6 hours\.$/);
  assert.equal(fromNow(-90 * 60000), "2 hours ago");
});

test("brightness comparisons use the Earth control sky and full moonlight", () => {
  assert.equal(brightness(10, 1, control), "As bright as Earth just after sunset");
  assert.equal(brightness(1, 3, control), "As bright as late dusk on Earth");
  assert.equal(brightness(4.5, 2, control), "About 15 times brighter than a full Moon on Earth");
  assert.equal(brightness(80000, 0, control), "As bright as a sunny day on Earth");
  const base = { sun_radius_deg: 0.27, solar: 1, earth: 0 };
  assert.equal(conditionOf({ ...base, sun_elevation_deg: -0.2, total: 4000 }), 0);
  assert.equal(conditionOf({ ...base, sun_elevation_deg: -20, total: 2.99 }), 1);
  assert.equal(conditionOf({ ...base, sun_elevation_deg: -50, total: 2.5, solar: 0.5, earth: 2 }), 2);
  assert.equal(conditionOf({ ...base, sun_elevation_deg: -50, total: 0.4 }), 3);
  assert.equal(conditionOf({ ...base, sun_elevation_deg: -80, total: 0.01 }), 4);
});

const places = new URL("../dist/assets/places.json", import.meta.url);
test("the six places carry their settings", { skip: !existsSync(places) && "page assets not built" }, () => {
  const list = JSON.parse(readFileSync(places)).places;
  assert.equal(list.length, 6);
  const byId = Object.fromEntries(list.map((p) => [p.id, p]));
  assert.equal(byId.tranquility.setting, "sea");
  assert.match(byId.tranquility.line, /467 m of sea/);
  assert.equal(byId.descartes.setting, "land");
  for (const p of list) assert.ok(Math.abs(p.latitude) <= 90 && Math.abs(p.longitude) <= 180, p.id);
});
