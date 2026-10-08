import {
  stateAt,
  geometry,
  julianDate,
  crossing,
  DAY_MS,
} from "./model/evaluator.mjs";
let astronomy, transfer;
const minMs = Date.UTC(2000, 0, 1),
  maxMs = Date.UTC(2501, 0, 1) - 1;
function evaluate(ms, settings, order = 12) {
  return stateAt(
    ms,
    settings.longitude,
    settings.latitude,
    astronomy,
    transfer,
    settings.earth,
    settings.scale,
    order,
  );
}
function geo(ms, s) {
  return geometry(julianDate(ms, s.scale), s.longitude, s.latitude, astronomy);
}

export function findCrossings(times, values, fn) {
  const roots = [];
  for (let i = 0; i < times.length - 1; i++) {
    if (values[i] === 0) roots.push(times[i]);
    if (values[i] * values[i + 1] < 0)
      roots.push(crossing(fn, times[i], times[i + 1]));
  }
  if (values.at(-1) === 0) roots.push(times.at(-1));
  // A grazing dip can cross twice between adjacent samples. Refine sampled
  // extrema before deciding that a same-sign interval contains no event.
  for (let i = 1; i < times.length - 1; i++) {
    const left = values[i - 1],
      middle = values[i],
      right = values[i + 1];
    if (
      (middle - left) * (right - middle) >= 0 ||
      left * middle <= 0 ||
      middle * right <= 0
    )
      continue;
    const sign = middle > 0 ? 1 : -1;
    let lo = times[i - 1],
      hi = times[i + 1];
    const ratio = (Math.sqrt(5) - 1) / 2;
    let x = hi - ratio * (hi - lo),
      y = lo + ratio * (hi - lo),
      fx = sign * fn(x),
      fy = sign * fn(y);
    while (hi - lo > 500) {
      if (fx < fy) {
        hi = y;
        y = x;
        fy = fx;
        x = hi - ratio * (hi - lo);
        fx = sign * fn(x);
      } else {
        lo = x;
        x = y;
        fx = fy;
        y = lo + ratio * (hi - lo);
        fy = sign * fn(y);
      }
    }
    const extremum = (lo + hi) / 2;
    if (sign * fn(extremum) < 0) {
      roots.push(
        crossing(fn, times[i - 1], extremum),
        crossing(fn, extremum, times[i + 1]),
      );
    }
  }
  return roots
    .filter((r) => r !== null)
    .sort((a, b) => a - b)
    .filter((r, i, all) => i === 0 || r - all[i - 1] > 1000)
    .filter(
      (r) =>
        fn(Math.max(times[0], r - 1000)) *
          fn(Math.min(times.at(-1), r + 1000)) <=
        0,
    );
}

export function calculate(settings, a, t) {
  if (
    !Number.isInteger(settings.year) ||
    settings.year < 2000 ||
    settings.year > 2500 ||
    !Number.isInteger(settings.month) ||
    settings.month < 0 ||
    settings.month > 11
  )
    throw Error("Calendar dates must lie between 2000 and 2500");
  astronomy = a;
  transfer = t;
  const start = Date.UTC(settings.year, settings.month, 1),
    end = Date.UTC(settings.year, settings.month + 1, 1),
    step = DAY_MS / 48;
  const series = [];
  for (let ms = start; ms <= end; ms += step)
    series.push(evaluate(Math.min(ms, maxMs), settings));
  const daily = [];
  for (let begin = start; begin < end; begin += DAY_MS) {
    const rows = series.filter((r) => r.ms >= begin && r.ms < begin + DAY_MS);
    const slices = [begin, begin + DAY_MS];
    const visible = (ms) => {
      const g = geo(Math.min(ms, maxMs), settings);
      return g.sun_elevation_deg + g.sun_radius_deg;
    };
    const times = Array.from({ length: 49 }, (_, k) => begin + k * step);
    slices.push(...findCrossings(times, times.map(visible), visible));
    slices.sort((a, b) => a - b);
    let up = 0;
    for (let i = 0; i < slices.length - 1; i++)
      if (visible((slices[i] + slices[i + 1]) / 2) > 0)
        up += slices[i + 1] - slices[i];
    daily.push({
      ms: begin,
      min: Math.min(...rows.map((r) => r.total)),
      max: Math.max(...rows.map((r) => r.total)),
      earth: rows[0].earth,
      sunHours: up / 3600000,
      strip: rows.filter((_, i) => i % 2 === 0).map((r) => r.total),
      eclipse: rows.some((r) => r.eclipse),
    });
  }
  // Event labels are geometric disk-edge crossings and source-sum lux thresholds.
  const defs = [
    {
      key: "sun",
      labelUp: "Sunrise begins",
      labelDown: "Sunset ends",
      value: (g) => g.sun_elevation_deg + g.sun_radius_deg,
    },
    {
      key: "earth",
      labelUp: "Earthrise begins",
      labelDown: "Earthset ends",
      value: (g) => g.earth_elevation_deg + g.earth_radius_deg,
    },
    ...[100, 10, 2.98, 0.1].map((level) => ({
      key: "light",
      level,
      labelUp: `Above ${level} lux`,
      labelDown: `Below ${level} lux`,
      value: (g) => g.total - level,
    })),
  ];
  const events = [];
  for (const def of defs) {
    const fn =
      def.key === "light"
        ? (ms) => def.value(evaluate(ms, settings))
        : (ms) => def.value(geo(ms, settings));
    for (const ms of findCrossings(
      series.map((r) => r.ms),
      series.map(def.value),
      fn,
    ))
      events.push({
        ms,
        key: def.key,
        level: def.level ?? null,
        label:
          fn(Math.min(ms + 1000, maxMs)) > fn(Math.max(ms - 1000, minMs))
            ? def.labelUp
            : def.labelDown,
      });
  }
  events.sort((a, b) => a.ms - b.ms);
  return { start, end, series, daily, events };
}

if (typeof self !== "undefined" && typeof self.postMessage === "function")
  self.onmessage = ({ data }) => {
    try {
      if (data.type === "init") {
        astronomy = data.astronomy;
        transfer = data.transfer;
        return;
      }
      const result = calculate(data.settings, astronomy, transfer);
      self.postMessage({ id: data.id, result });
    } catch (error) {
      self.postMessage({ id: data.id, error: error.message });
    }
  };
