// Portable reader for versioned calendar products. Physics/provenance live in illumination.
const rad = Math.PI / 180,
  deg = 180 / Math.PI,
  day = 86400000;
const clamp = (x, a, b) => Math.max(a, Math.min(b, x));
const dot = (a, b) => a.reduce((s, v, i) => s + v * b[i], 0);
const norm = (v) => Math.hypot(...v);
const add = (a, b) => a.map((v, i) => v + b[i]);
const sub = (a, b) => a.map((v, i) => v - b[i]);
const mul = (a, k) => a.map((v) => v * k);
const unit = (v) => mul(v, 1 / norm(v));
const cross = (a, b) => [
  a[1] * b[2] - a[2] * b[1],
  a[2] * b[0] - a[0] * b[2],
  a[0] * b[1] - a[1] * b[0],
];
const direction = (lon, lat) => [
  Math.cos(lat) * Math.cos(lon),
  Math.cos(lat) * Math.sin(lon),
  Math.sin(lat),
];

export function bracket(a, x) {
  if (x <= a[0]) return [0, 0];
  if (x >= a.at(-1)) return [a.length - 2, 1];
  let lo = 0,
    hi = a.length - 1;
  while (hi - lo > 1) {
    const m = (lo + hi) >> 1;
    if (a[m] <= x) lo = m;
    else hi = m;
  }
  return [lo, (x - a[lo]) / (a[lo + 1] - a[lo])];
}
export function interpolate(a, v, x) {
  const [i, w] = bracket(a, x);
  return v[i] * (1 - w) + v[i + 1] * w;
}

export function ttOffset(ms) {
  const changes = [
    ["2000-01-01", 64.184],
    ["2006-01-01", 65.184],
    ["2009-01-01", 66.184],
    ["2012-07-01", 67.184],
    ["2015-07-01", 68.184],
    ["2017-01-01", 69.184],
  ];
  let result = 64.184;
  for (const [date, value] of changes)
    if (ms >= Date.parse(date + "T00:00:00Z")) result = value;
  return result;
}
export function julianDate(ms, scale = "UTC") {
  return ms / day + 2440587.5 + (scale === "TT" ? 0 : ttOffset(ms) / 86400);
}

export function ephemeris(jd, data) {
  const t = (jd - 2451545) / 36525;
  const polynomial = (a, b, c, d, e) =>
    rad * (a + t * (b + t * (c + t * (d + t * e))));
  const L = polynomial(
      218.3164477,
      481267.88123421,
      -0.0015786,
      1 / 538841,
      -1 / 65194000,
    ),
    D = polynomial(
      297.8501921,
      445267.1114034,
      -0.0018819,
      1 / 545868,
      -1 / 113065000,
    ),
    M = polynomial(357.5291092, 35999.0502909, -0.0001536, 1 / 24490000, 0),
    Mp = polynomial(
      134.9633964,
      477198.8675055,
      0.0087414,
      1 / 69699,
      -1 / 14712000,
    ),
    F = polynomial(
      93.272095,
      483202.0175233,
      -0.0036539,
      -1 / 3526000,
      1 / 863310000,
    ),
    node = polynomial(
      125.0445479,
      -1934.1362891,
      0.0020754,
      1 / 467441,
      -1 / 60616000,
    );
  const eccentricity = 1 - 0.002516 * t - 0.0000074 * t * t;
  const series = (table, column, trig) =>
    table.reduce(
      (s, r) =>
        s +
        r[column] *
          eccentricity ** Math.abs(r[1]) *
          trig(r[0] * D + r[1] * M + r[2] * Mp + r[3] * F),
      0,
    );
  const a1 = rad * (119.75 + 131.849 * t),
    a2 = rad * (53.09 + 479264.29 * t),
    a3 = rad * (313.45 + 481266.484 * t);
  const lon =
    L +
    rad *
      1e-6 *
      (series(data.longitude_distance, 4, Math.sin) +
        3958 * Math.sin(a1) +
        1962 * Math.sin(L - F) +
        318 * Math.sin(a2));
  const lat =
    rad *
    1e-6 *
    (series(data.latitude, 4, Math.sin) -
      2235 * Math.sin(L) +
      382 * Math.sin(a3) +
      175 * Math.sin(a1 - F) +
      175 * Math.sin(a1 + F) +
      127 * Math.sin(L - Mp) -
      115 * Math.sin(L + Mp));
  const distance =
    385000.56 + series(data.longitude_distance, 5, Math.cos) / 1000;
  const solarMean = rad * (280.46646 + 36000.76983 * t + 0.0003032 * t * t),
    anomaly = rad * (357.52911 + 35999.05029 * t - 0.0001537 * t * t);
  const ecc = 0.016708634 - 0.000042037 * t - 0.0000001267 * t * t;
  const centre =
    rad *
    ((1.914602 - 0.004817 * t - 0.000014 * t * t) * Math.sin(anomaly) +
      (0.019993 - 0.000101 * t) * Math.sin(2 * anomaly) +
      0.000289 * Math.sin(3 * anomaly));
  const solarLon = solarMean + centre,
    solarDistance =
      (((1.000001018 * (1 - ecc * ecc)) /
        (1 + ecc * Math.cos(anomaly + centre))) *
        data.au_m) /
      1000;
  const selenographic = (lo, la) => {
    const w = lo - node,
      i = data.equator_to_ecliptic;
    const a = Math.atan2(
      Math.sin(w) * Math.cos(la) * Math.cos(i) - Math.sin(la) * Math.sin(i),
      Math.cos(w) * Math.cos(la),
    );
    return direction(
      a - F,
      Math.asin(
        clamp(
          -Math.sin(w) * Math.cos(la) * Math.sin(i) -
            Math.sin(la) * Math.cos(i),
          -1,
          1,
        ),
      ),
    );
  };
  const relative = sub(
    mul(direction(lon, lat), distance),
    mul(direction(solarLon, 0), solarDistance),
  );
  const ds = norm(relative);
  return {
    earth: selenographic(lon, lat),
    earthDistance: distance * 1000,
    sun: selenographic(
      Math.atan2(relative[1], relative[0]),
      Math.asin(relative[2] / ds),
    ),
    sunDistance: ds * 1000,
  };
}

export function geometry(jd, longitude, latitude, data) {
  if (
    !Number.isFinite(jd) ||
    latitude < -90 ||
    latitude > 90 ||
    longitude < -180 ||
    longitude > 180 ||
    !Number.isFinite(latitude + longitude)
  )
    throw Error("Invalid date or coordinates");
  const g = ephemeris(jd, data),
    lon = longitude * rad,
    lat = latitude * rad;
  const frame = [
    [-Math.sin(lon), Math.cos(lon), 0],
    [
      -Math.sin(lat) * Math.cos(lon),
      -Math.sin(lat) * Math.sin(lon),
      Math.cos(lat),
    ],
    direction(lon, lat),
  ];
  const local = (v) => frame.map((row) => dot(row, v));
  const observer = mul(frame[2], data.moon_radius_m),
    earthVector = mul(g.earth, g.earthDistance),
    sunVector = mul(g.sun, g.sunDistance);
  const earth = local(sub(earthVector, observer)),
    sun = local(sub(sunVector, observer)),
    de = norm(earth),
    ds = norm(sun),
    eu = unit(earth),
    su = unit(sun);
  const earthToSun = sub(sunVector, earthVector),
    des = norm(earthToSun),
    se = local(unit(earthToSun));
  const phase = Math.acos(clamp(dot(se, mul(eu, -1)), -1, 1));
  let vertical = sub([0, 0, 1], mul(eu, eu[2]));
  vertical = norm(vertical) > 1e-12 ? unit(vertical) : [0, 1, 0];
  const horizontal = cross(vertical, eu);
  const out = {
    earth_distance_m: de,
    sun_distance_m: ds,
    earth_phase_deg: phase * deg,
    earth_lit_fraction: (1 + Math.cos(phase)) / 2,
    earth_bright_limb_rad: Math.atan2(dot(se, horizontal), dot(se, vertical)),
    earth_radius_deg: Math.asin(data.earth_radius_m / de) * deg,
    sun_radius_deg: Math.asin(data.sun_radius_m / ds) * deg,
    sun_distance_factor: (data.au_m / ds) ** 2,
    earth_sunlight_factor: (data.au_m / des) ** 2,
    source_separation_deg: Math.acos(clamp(dot(eu, su), -1, 1)) * deg,
  };
  for (const [name, v] of [
    ["sun", su],
    ["earth", eu],
  ]) {
    out[name + "_elevation_deg"] = Math.asin(clamp(v[2], -1, 1)) * deg;
    out[name + "_azimuth_deg"] = (Math.atan2(v[0], v[1]) * deg + 360) % 360;
  }
  return out;
}

const gaussCache = new Map();
function gauss(n) {
  if (gaussCache.has(n)) return gaussCache.get(n);
  const nodes = [],
    weights = [];
  for (let i = 0; i < n; i++) {
    let x = Math.cos((Math.PI * (i + 0.75)) / (n + 0.5)),
      derivative;
    for (let j = 0; j < 30; j++) {
      let a = 1,
        b = x;
      for (let k = 2; k <= n; k++) {
        const c = ((2 * k - 1) * x * b - (k - 1) * a) / k;
        a = b;
        b = c;
      }
      derivative = (n * (x * b - a)) / (x * x - 1);
      const delta = b / derivative;
      x -= delta;
      if (Math.abs(delta) < 1e-15) break;
    }
    nodes.push(x);
    weights.push(2 / ((1 - x * x) * derivative * derivative));
  }
  const result = { nodes, weights };
  gaussCache.set(n, result);
  return result;
}
export function diskIntegral(
  callback,
  elevation,
  radius,
  phase = null,
  beta = 0,
  order = 12,
) {
  if (phase !== null && phase >= 180) return [0, 0];
  const { nodes, weights } = gauss(order),
    alpha = (phase ?? 0) * rad,
    se = Math.sin(elevation * rad),
    ce = Math.cos(elevation * rad),
    sr = Math.sin(radius * rad);
  let a = 0,
    b = 0,
    total = 0;
  for (let i = 0; i < order; i++) {
    const y = nodes[i],
      half = Math.sqrt(1 - y * y),
      lower = phase === null ? -half : -Math.cos(alpha) * half,
      width = half - lower;
    for (let j = 0; j < order; j++) {
      const x = lower + ((nodes[j] + 1) * width) / 2,
        z = Math.sqrt(Math.max(0, 1 - x * x - y * y)),
        directionCos = Math.sqrt(1 - sr * sr * (x * x + y * y));
      let w = (weights[i] * weights[j] * width) / 2 / directionCos;
      if (phase !== null)
        w *= Math.max(0, x * Math.sin(alpha) + z * Math.cos(alpha));
      const vertical = x * Math.cos(beta) + y * Math.sin(beta),
        e =
          Math.asin(clamp(se * directionCos + ce * sr * vertical, -1, 1)) * deg;
      const v = callback(e);
      a += w * v[0];
      b += w * v[1];
      total += w;
    }
  }
  return total > 0 ? [a / total, b / total] : [0, 0];
}

export function lightAt(g, transfer, includeEarth = true, order = 12) {
  if (transfer.schema !== "terluna.illumination.calendar-transfer/1")
    throw Error("Unsupported lighting data");
  const da = transfer.direct_elevation_deg,
    fa = transfer.diffuse_elevation_deg;
  const s = diskIntegral(
    (e) => [
      interpolate(da, transfer.solar.direct_lux, e),
      interpolate(fa, transfer.solar.diffuse_lux, e),
    ],
    g.sun_elevation_deg,
    g.sun_radius_deg,
    null,
    0,
    order,
  ).map((v) => v * g.sun_distance_factor);
  let earth = [0, 0];
  if (includeEarth && g.earth_phase_deg < 180) {
    const [p, w] = bracket(transfer.phase_deg, g.earth_phase_deg);
    const sample = (rows, angles, e) =>
      (1 - w) * interpolate(angles, rows[p], e) +
      w * interpolate(angles, rows[p + 1], e);
    earth = diskIntegral(
      (e) => [
        sample(transfer.earth.direct_lux, da, e),
        sample(transfer.earth.diffuse_lux, fa, e),
      ],
      g.earth_elevation_deg,
      g.earth_radius_deg,
      g.earth_phase_deg,
      g.earth_bright_limb_rad,
      order,
    );
    const R = transfer.constants.earth_radius_m,
      reference = transfer.constants.earth_reference_distance_m;
    const omega = (d) =>
      (2 * Math.PI * (R / d) ** 2) / (1 + Math.sqrt(1 - (R / d) ** 2));
    earth = earth.map(
      (v) =>
        ((v * omega(g.earth_distance_m)) / omega(reference)) *
        g.earth_sunlight_factor,
    );
  }
  const eclipse =
    g.source_separation_deg < g.earth_radius_deg + g.sun_radius_deg;
  return {
    solar_direct: s[0],
    solar_diffuse: s[1],
    solar: s[0] + s[1],
    earth_direct: earth[0],
    earth_diffuse: earth[1],
    earth: earth[0] + earth[1],
    total: s[0] + s[1] + earth[0] + earth[1],
    eclipse,
    earthIncluded: includeEarth,
  };
}

export function stateAt(
  ms,
  longitude,
  latitude,
  astronomy,
  transfer,
  includeEarth = true,
  scale = "UTC",
  order = 12,
) {
  const g = geometry(julianDate(ms, scale), longitude, latitude, astronomy);
  return { ms, ...g, ...lightAt(g, transfer, includeEarth, order) };
}

export function crossing(fn, start, end, target = 0, toleranceMs = 500) {
  let a = fn(start) - target,
    b = fn(end) - target;
  if (a === 0) return start;
  if (b === 0) return end;
  if (a * b > 0) return null;
  while (end - start > toleranceMs) {
    const mid = (start + end) / 2,
      c = fn(mid) - target;
    if (c === 0) return mid;
    if (c * a > 0) {
      start = mid;
      a = c;
    } else {
      end = mid;
      b = c;
    }
  }
  return (start + end) / 2;
}

export const DAY_MS = day;
