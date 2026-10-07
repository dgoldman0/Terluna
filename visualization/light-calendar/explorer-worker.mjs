// Calculations for the explorer page: spans of samples with their events,
// month calendars and a year of hourly light. All light comes from the shared
// evaluator; this module samples it and finds crossings and turning points.
import {
  stateAt,
  geometry,
  julianDate,
  crossing,
  DAY_MS,
} from "./model/evaluator.mjs";
import { calculate, findCrossings } from "./worker.mjs";

export const PRACTICAL_DUSK_LUX = 2.98;
const HOUR = 3600000;
const minMs = Date.UTC(2000, 0, 1),
  maxMs = Date.UTC(2501, 0, 1) - 1;
let astronomy, transfer;

export function regimeCode(s) {
  if (s.sun_elevation_deg + s.sun_radius_deg > 0) return 0; // day
  if (s.earth > s.solar && s.total >= 0.05) return 2; // Earthlit night
  if (s.total >= PRACTICAL_DUSK_LUX) return 1; // twilight
  if (s.total >= 0.05) return 3; // night with a twilight glow
  return 4; // dark night
}

const clampMs = (ms) => Math.min(maxMs, Math.max(minMs, ms));
function state(ms, settings, order) {
  return stateAt(
    clampMs(ms),
    settings.longitude,
    settings.latitude,
    astronomy,
    transfer,
    settings.earth,
    settings.scale,
    order,
  );
}
function geo(ms, settings) {
  return geometry(
    julianDate(clampMs(ms), settings.scale),
    settings.longitude,
    settings.latitude,
    astronomy,
  );
}

// Golden-section refinement of a sampled turning point of fn.
function turningPoint(fn, lo, hi, sign) {
  const r = (Math.sqrt(5) - 1) / 2;
  let x = hi - r * (hi - lo),
    y = lo + r * (hi - lo),
    fx = sign * fn(x),
    fy = sign * fn(y);
  while (hi - lo > 60000) {
    if (fx > fy) {
      hi = y;
      y = x;
      fy = fx;
      x = hi - r * (hi - lo);
      fx = sign * fn(x);
    } else {
      lo = x;
      x = y;
      fx = fy;
      y = lo + r * (hi - lo);
      fy = sign * fn(y);
    }
  }
  return (lo + hi) / 2;
}

export function span(settings, start, end, step = HOUR, order = 8) {
  const times = [];
  for (let t = start; t <= end; t += step) times.push(clampMs(t));
  const states = times.map((t) => state(t, settings, order));
  const events = [];
  const add = (ms, type, extra = {}) => events.push({ ms, type, ...extra });
  const rising = (fn, ms) =>
    fn(Math.min(ms + 1000, maxMs)) > fn(Math.max(ms - 1000, minMs));

  const sunFn = (ms) => {
    const g = geo(ms, settings);
    return g.sun_elevation_deg + g.sun_radius_deg;
  };
  for (const ms of findCrossings(
    times,
    states.map((s) => s.sun_elevation_deg + s.sun_radius_deg),
    sunFn,
  ))
    add(ms, rising(sunFn, ms) ? "sunrise" : "sunset");

  const earthFn = (ms) => {
    const g = geo(ms, settings);
    return g.earth_elevation_deg + g.earth_radius_deg;
  };
  for (const ms of findCrossings(
    times,
    states.map((s) => s.earth_elevation_deg + s.earth_radius_deg),
    earthFn,
  ))
    add(ms, rising(earthFn, ms) ? "earthrise" : "earthset");

  const lightFn = (ms) => state(ms, settings, order).total - PRACTICAL_DUSK_LUX;
  for (const ms of findCrossings(
    times,
    states.map((s) => s.total - PRACTICAL_DUSK_LUX),
    lightFn,
  ))
    add(ms, rising(lightFn, ms) ? "dawn" : "dusk");

  const litFn = (ms) => geo(ms, settings).earth_lit_fraction;
  for (let i = 1; i < states.length - 1; i++) {
    const a = states[i - 1].earth_lit_fraction,
      b = states[i].earth_lit_fraction,
      c = states[i + 1].earth_lit_fraction;
    if (b >= a && b > c) {
      const ms = turningPoint(litFn, times[i - 1], times[i + 1], 1);
      add(ms, "full-earth", { lit: litFn(ms) });
    } else if (b <= a && b < c) {
      const ms = turningPoint(litFn, times[i - 1], times[i + 1], -1);
      add(ms, "new-earth", { lit: litFn(ms) });
    }
  }

  const overlapFn = (ms) => {
    const g = geo(ms, settings);
    return g.earth_radius_deg + g.sun_radius_deg - g.source_separation_deg;
  };
  const overlap = times.map(overlapFn);
  for (let i = 0; i < times.length - 1; i++) {
    if (overlap[i] * overlap[i + 1] < 0) {
      const ms = crossing(overlapFn, times[i], times[i + 1]);
      if (ms !== null)
        add(ms, overlap[i + 1] > 0 ? "eclipse-begins" : "eclipse-ends");
    }
  }
  // An eclipse shorter than the step can hide between samples; refine the
  // closest approach where separation dips without a sampled overlap.
  for (let i = 1; i < times.length - 1; i++) {
    if (
      overlap[i] >= overlap[i - 1] &&
      overlap[i] >= overlap[i + 1] &&
      overlap[i] < 0 &&
      overlap[i] > -2
    ) {
      const peak = turningPoint(overlapFn, times[i - 1], times[i + 1], 1);
      if (overlapFn(peak) > 0) {
        add(crossing(overlapFn, times[i - 1], peak), "eclipse-begins");
        add(crossing(overlapFn, peak, times[i + 1]), "eclipse-ends");
      }
    }
  }
  events.sort((a, b) => a.ms - b.ms);
  const samples = {
    ms: times,
    total: states.map((s) => s.total),
    solar: states.map((s) => s.solar),
    earth: states.map((s) => s.earth),
    sunEl: states.map((s) => s.sun_elevation_deg),
    earthEl: states.map((s) => s.earth_elevation_deg),
    lit: states.map((s) => s.earth_lit_fraction),
    regime: states.map(regimeCode),
  };
  return { start, end, step, samples, events };
}

// The Earth calendar month with the explorer's day markers.
export function month(settings) {
  const result = calculate(settings, astronomy, transfer);
  const extra = span(
    settings,
    result.start - 2 * HOUR,
    result.end + 2 * HOUR,
    HOUR,
    6,
  );
  const events = extra.events.filter(
    (e) => e.ms >= result.start && e.ms < result.end,
  );
  for (const day of result.daily) {
    day.events = events.filter((e) => e.ms >= day.ms && e.ms < day.ms + DAY_MS);
    day.regimes = result.series
      .filter((r) => r.ms >= day.ms && r.ms < day.ms + DAY_MS)
      .filter((_, i) => i % 2 === 0)
      .map(regimeCode);
  }
  return { ...result, markers: events };
}

// Hourly light for one Earth year.
export function year(settings) {
  const start = Date.UTC(settings.year, 0, 1),
    end = Date.UTC(settings.year + 1, 0, 1);
  const days = Math.round((end - start) / DAY_MS);
  const total = new Float32Array(days * 24),
    regime = new Uint8Array(days * 24);
  for (let d = 0; d < days; d++)
    for (let h = 0; h < 24; h++) {
      const s = state(start + d * DAY_MS + (h + 0.5) * HOUR, settings, 4);
      total[d * 24 + h] = s.total;
      regime[d * 24 + h] = regimeCode(s);
    }
  return { start, days, total, regime };
}

export function setProducts(a, t) {
  astronomy = a;
  transfer = t;
}

if (typeof self !== "undefined" && typeof self.postMessage === "function")
  self.onmessage = ({ data }) => {
    try {
      if (data.type === "init") {
        setProducts(data.astronomy, data.transfer);
        return;
      }
      let result;
      if (data.type === "span")
        result = span(data.settings, data.start, data.end, data.step, data.order);
      else if (data.type === "month") result = month(data.settings);
      else if (data.type === "year") result = year(data.settings);
      else throw Error(`Unknown request ${data.type}`);
      self.postMessage({ id: data.id, type: data.type, result });
    } catch (error) {
      self.postMessage({ id: data.id, type: data.type, error: error.message });
    }
  };
