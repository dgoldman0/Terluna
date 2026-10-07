import { stateAt, ephemeris, julianDate, DAY_MS } from "./model/evaluator.mjs";
const $ = (id) => document.getElementById(id),
  svgNS = "http://www.w3.org/2000/svg";
const minMs = Date.UTC(2000, 0, 1),
  maxMs = Date.UTC(2501, 0, 1) - 60000;
const clamp = (n, a, b) => Math.min(b, Math.max(a, n));
let astronomy,
  transfer,
  validation,
  worker,
  request = 0,
  calendar,
  selected,
  drawFrame = 0;
const query = new URLSearchParams(location.search),
  requested = Date.parse(query.get("at") ?? "");
let ms = Number.isFinite(requested)
  ? clamp(requested, minMs, maxMs)
  : Date.now();
let settings = {
  longitude: query.has("lon") ? Number(query.get("lon")) : -16.125,
  latitude: query.has("lat") ? Number(query.get("lat")) : -23.125,
  earth: query.get("earth") !== "0",
  scale: query.get("scale") === "TT" ? "TT" : "UTC",
  year: new Date(ms).getUTCFullYear(),
  month: new Date(ms).getUTCMonth(),
};
if (!Number.isFinite(settings.longitude) || Math.abs(settings.longitude) > 180)
  settings.longitude = -16.125;
if (!Number.isFinite(settings.latitude) || Math.abs(settings.latitude) > 90)
  settings.latitude = -23.125;
const formatDate = (value, options = {}) =>
  new Intl.DateTimeFormat("en-US", { timeZone: "UTC", ...options }).format(
    new Date(value),
  );
const number = (v, d = 0) =>
  v.toLocaleString("en-US", {
    maximumFractionDigits: d,
    minimumFractionDigits: 0,
  });
function lux(v) {
  if (v === 0) return "0 lux";
  if (v < 0.001) return `${v.toExponential(1)} lux`;
  if (v < 1) return `${number(v, 3)} lux`;
  if (v < 100) return `${number(v, 2)} lux`;
  return `${number(v, 0)} lux`;
}
function bearing(a) {
  return [
    "N",
    "NNE",
    "NE",
    "ENE",
    "E",
    "ESE",
    "SE",
    "SSE",
    "S",
    "SSW",
    "SW",
    "WSW",
    "W",
    "WNW",
    "NW",
    "NNW",
  ][Math.round(a / 22.5) % 16];
}
function el(tag, attrs = {}, text) {
  const n = document.createElementNS(svgNS, tag);
  for (const [k, v] of Object.entries(attrs)) n.setAttribute(k, String(v));
  if (text !== undefined) n.textContent = text;
  return n;
}
function notice(message) {
  $("toast").textContent = message;
  $("toast").hidden = false;
  setTimeout(() => ($("toast").hidden = true), 3000);
}
function dateInput() {
  const iso = new Date(ms).toISOString();
  $("date").value = iso.slice(0, 10);
  $("time").value = iso.slice(11, 16);
}
function updateQuery() {
  const p = new URLSearchParams({
    lat: settings.latitude,
    lon: settings.longitude,
    at: new Date(ms).toISOString(),
    earth: settings.earth ? "1" : "0",
    scale: settings.scale,
  });
  history.replaceState(null, "", `${location.pathname}?${p}${location.hash}`);
}
function pickPreset() {
  let site = astronomy.sites.find(
    (s) =>
      Math.abs(s.latitude - settings.latitude) < 1e-5 &&
      Math.abs(s.longitude - settings.longitude) < 1e-5,
  );
  $("place").value = site?.id ?? "custom";
  $("place-description").textContent =
    site?.description ??
    "A custom point on the Open Moon. Negative coordinates indicate south and west.";
  return site;
}
function sync() {
  dateInput();
  $("latitude").value = settings.latitude;
  $("longitude").value = settings.longitude;
  $("earth-toggle").checked = settings.earth;
  $("timescale").value = settings.scale;
  $("month").value =
    `${settings.year}-${String(settings.month + 1).padStart(2, "0")}`;
  pickPreset();
  updateQuery();
}
function regime(s) {
  if (s.sun_elevation_deg + s.sun_radius_deg > 0) return "Sun above horizon";
  if (s.earth > s.solar) return "Earthlit night";
  if (s.total >= 100) return "Bright twilight";
  if (s.total >= 2.98) return "Long twilight";
  if (s.total >= 0.1) return "Faint twilight";
  return "Low light";
}
function lightDescription(s) {
  if (s.eclipse) return "Unobscured reference during an eclipse";
  if (!settings.earth) return "Sunlight and scattered sunlight";
  if (s.earth > s.solar)
    return `${number((100 * s.earth) / s.total, 0)}% of this light comes from Earth`;
  if (s.solar_diffuse > s.solar_direct)
    return "Scattered sunlight fills the sky";
  return "Direct sunlight with a bright diffuse sky";
}

function updateSnapshot() {
  selected = stateAt(
    ms,
    settings.longitude,
    settings.latitude,
    astronomy,
    transfer,
    settings.earth,
    settings.scale,
  );
  const site = pickPreset();
  $("moment-label").textContent =
    `${formatDate(ms, { month: "short", day: "numeric", year: "numeric" })} · ${formatDate(ms, { hour: "2-digit", minute: "2-digit", hour12: false })} ${settings.scale}`;
  $("snapshot-title").textContent =
    site?.name ??
    `${number(Math.abs(settings.latitude), 2)}° ${settings.latitude < 0 ? "S" : "N"}, ${number(Math.abs(settings.longitude), 2)}° ${settings.longitude < 0 ? "W" : "E"}`;
  $("regime").textContent = regime(selected);
  $("regime").style.color =
    selected.earth > selected.solar ? "var(--earth)" : "var(--sun)";
  const value =
    selected.total >= 10000 ? selected.total / 1000 : selected.total;
  $("total-lux").textContent =
    value > 0 && value < 0.001
      ? value.toExponential(1)
      : number(value, value < 1 ? 3 : value < 100 ? 1 : 0);
  $("lux-unit").textContent = selected.total >= 10000 ? "klux" : "lux";
  $("light-description").textContent = lightDescription(selected);
  $("sun-lux").textContent = lux(selected.solar);
  $("earth-lux").textContent = settings.earth ? lux(selected.earth) : "Off";
  for (const [id, key] of [
    ["sun-direct", "solar_direct"],
    ["sun-diffuse", "solar_diffuse"],
    ["earth-direct", "earth_direct"],
    ["earth-diffuse", "earth_diffuse"],
  ])
    $(id).textContent = lux(selected[key]);
  const f = selected.total > 0 ? selected.solar / selected.total : 1;
  $("solar-bar").style.width = `${f * 100}%`;
  $("earth-bar").style.width = `${(1 - f) * 100}%`;
  for (const name of ["sun", "earth"]) {
    const e = selected[name + "_elevation_deg"],
      a = selected[name + "_azimuth_deg"];
    $(name + "-position").textContent =
      `${number(Math.abs(e), 1)}° ${e >= 0 ? "above" : "below"}`;
    $(name + "-bearing").textContent =
      `${bearing(a)} · ${number(a, 1)}° azimuth`;
  }
  $("earth-fraction").textContent =
    `${number(selected.earth_lit_fraction * 100, 1)}%`;
  $("earth-diameter").textContent =
    `${number(selected.earth_radius_deg * 2, 2)}° apparent diameter`;
  $("eclipse-notice").hidden = !selected.eclipse;
  $("earth-legend").style.opacity = settings.earth ? "1" : ".35";
  $("event-scale").textContent = settings.scale;
  $("daily-scale").textContent = settings.scale;
  $("scrubber").setAttribute("aria-valuetext", $("moment-label").textContent);
  drawEarth();
  drawSky();
  drawMap();
  drawMarker();
  renderEvents();
}

function drawEarth() {
  const c = $("earth-phase"),
    ctx = c.getContext("2d"),
    image = ctx.createImageData(80, 80),
    a = (selected.earth_phase_deg * Math.PI) / 180,
    beta = selected.earth_bright_limb_rad;
  for (let y = 0; y < 80; y++)
    for (let x = 0; x < 80; x++) {
      const xx = (x - 39.5) / 36,
        yy = (39.5 - y) / 36,
        r2 = xx * xx + yy * yy;
      if (r2 > 1) continue;
      const z = Math.sqrt(1 - r2),
        // The model's horizontal tangent is screen-left (up × body direction).
        toward = -xx * Math.sin(beta) + yy * Math.cos(beta),
        illum = Math.max(0, toward * Math.sin(a) + z * Math.cos(a));
      const v = 0.06 + 0.94 * Math.sqrt(illum),
        i = 4 * (y * 80 + x);
      image.data.set([100 * v, 200 * v, 185 * v, 255], i);
    }
  ctx.clearRect(0, 0, 80, 80);
  ctx.putImageData(image, 0, 0);
  c.setAttribute(
    "aria-label",
    `Earth is ${number(selected.earth_lit_fraction * 100, 1)} percent illuminated`,
  );
}

// Map: the water mask is resampled at the canvas's own pixel resolution
// (bilinear between 1° cells, 2×2 supersampled) and cached per canvas size;
// only the day/night shading is recomputed when the moment changes.
const mapCache = { field: null, key: "", base: null };
const smooth = (a, b, x) => {
  const t = clamp((x - a) / (b - a), 0, 1);
  return t * t * (3 - 2 * t);
};
function waterField() {
  if (mapCache.field) return mapCache.field;
  const lat = astronomy.map_grid.latitude_centres,
    lon = astronomy.map_grid.longitude_centres,
    rows = lat.length,
    cols = lon.length,
    water = new Float32Array(rows * cols);
  for (const [row, a, b] of astronomy.water_spans)
    water.fill(1, row * cols + a, row * cols + Math.min(b, cols));
  mapCache.field = {
    water,
    rows,
    cols,
    lat0: lat[0],
    dlat: lat[0] - lat[1],
    lon0: lon[0],
    dlon: lon[1] - lon[0],
  };
  return mapCache.field;
}
function waterAt(F, latDeg, lonDeg) {
  const fr = clamp((F.lat0 - latDeg) / F.dlat, 0, F.rows - 1),
    fc = ((((lonDeg - F.lon0) / F.dlon) % F.cols) + F.cols) % F.cols,
    r0 = Math.min(F.rows - 2, Math.floor(fr)),
    tr = fr - r0,
    c0 = Math.floor(fc) % F.cols,
    tc = fc - Math.floor(fc),
    c1 = (c0 + 1) % F.cols,
    w = F.water,
    a = r0 * F.cols,
    b = a + F.cols;
  return (
    (w[a + c0] * (1 - tc) + w[a + c1] * tc) * (1 - tr) +
    (w[b + c0] * (1 - tc) + w[b + c1] * tc) * tr
  );
}
function mapBase(w, h) {
  const key = `${w}x${h}`;
  if (mapCache.key === key) return mapCache.base;
  const F = waterField(),
    base = new Float32Array(w * h * 3),
    sea = [52, 121, 168],
    land = [116, 126, 96],
    coast = [150, 196, 222],
    // Half a pixel of transition in field units keeps coastlines smooth.
    e = clamp((0.5 * 360) / w, 0.08, 0.5);
  for (let y = 0; y < h; y++)
    for (let x = 0; x < w; x++) {
      let t = 0,
        edge = 0;
      for (const [sx, sy] of [
        [0.25, 0.25],
        [0.75, 0.25],
        [0.25, 0.75],
        [0.75, 0.75],
      ]) {
        const v = waterAt(
          F,
          90 - ((y + sy) / h) * 180,
          ((x + sx) / w) * 360 - 180,
        );
        t += smooth(0.5 - e, 0.5 + e, v) / 4;
        edge += Math.max(0, 1 - Math.abs(v - 0.5) / (2 * e)) / 4;
      }
      const i = 3 * (y * w + x),
        k = 0.3 * edge;
      for (let c = 0; c < 3; c++)
        base[i + c] =
          (land[c] * (1 - t) + sea[c] * t) * (1 - k) + coast[c] * k;
    }
  mapCache.key = key;
  mapCache.base = base;
  return base;
}
function drawMap() {
  const c = $("world-map"),
    cssWidth = c.getBoundingClientRect().width,
    dpr = Math.min(window.devicePixelRatio || 1, 2),
    w = cssWidth ? clamp(Math.round(cssWidth * dpr), 360, 1200) : c.width,
    h = Math.round(w / 2),
    s = cssWidth ? w / cssWidth : 1;
  if (c.width !== w || c.height !== h) {
    c.width = w;
    c.height = h;
  }
  const ctx = c.getContext("2d"),
    base = mapBase(w, h),
    image = ctx.createImageData(w, h),
    out = image.data,
    sun = ephemeris(julianDate(ms, settings.scale), astronomy).sun,
    night = [5, 10, 24],
    px = (2 * Math.PI) / w,
    cosLon = new Float32Array(w),
    sinLon = new Float32Array(w);
  for (let x = 0; x < w; x++) {
    const lon = (((x + 0.5) / w) * 360 - 180) * (Math.PI / 180);
    cosLon[x] = Math.cos(lon);
    sinLon[x] = Math.sin(lon);
  }
  for (let y = 0; y < h; y++) {
    const lat = (90 - ((y + 0.5) / h) * 180) * (Math.PI / 180),
      cl = Math.cos(lat),
      sl = Math.sin(lat) * sun[2];
    for (let x = 0; x < w; x++) {
      const mu = cl * (cosLon[x] * sun[0] + sinLon[x] * sun[1]) + sl,
        a =
          (0.25 + 0.45 * Math.min(1, Math.max(0, -mu) * 3)) *
          clamp(0.5 - mu / px, 0, 1),
        i = y * w + x;
      for (let k = 0; k < 3; k++)
        out[4 * i + k] = base[3 * i + k] * (1 - a) + night[k] * a;
      out[4 * i + 3] = 255;
    }
  }
  ctx.putImageData(image, 0, 0);
  ctx.lineWidth = Math.max(1, s * 0.75);
  for (let deg = -150; deg <= 150; deg += 30) {
    const x = ((deg + 180) / 360) * w;
    ctx.strokeStyle = deg === 0 ? "#dce8f040" : "#dce8f01c";
    ctx.beginPath();
    ctx.moveTo(x, 0);
    ctx.lineTo(x, h);
    ctx.stroke();
  }
  for (let deg = -60; deg <= 60; deg += 30) {
    const y = ((90 - deg) / 180) * h;
    ctx.strokeStyle = deg === 0 ? "#dce8f040" : "#dce8f01c";
    ctx.beginPath();
    ctx.moveTo(0, y);
    ctx.lineTo(w, y);
    ctx.stroke();
  }
  const x = ((settings.longitude + 180) / 360) * w,
    y = ((90 - settings.latitude) / 180) * h,
    gap = 11 * s;
  ctx.save();
  ctx.strokeStyle = "#ffffff70";
  ctx.lineWidth = s;
  ctx.setLineDash([3 * s, 3 * s]);
  ctx.beginPath();
  ctx.moveTo(0, y);
  ctx.lineTo(Math.max(0, x - gap), y);
  ctx.moveTo(Math.min(w, x + gap), y);
  ctx.lineTo(w, y);
  ctx.moveTo(x, 0);
  ctx.lineTo(x, Math.max(0, y - gap));
  ctx.moveTo(x, Math.min(h, y + gap));
  ctx.lineTo(x, h);
  ctx.stroke();
  ctx.restore();
  ctx.beginPath();
  ctx.arc(x, y, 7 * s, 0, Math.PI * 2);
  ctx.lineWidth = 5 * s;
  ctx.strokeStyle = "#050a10c0";
  ctx.stroke();
  ctx.lineWidth = 2.2 * s;
  ctx.strokeStyle = "#ffffff";
  ctx.stroke();
  ctx.beginPath();
  ctx.arc(x, y, 2.8 * s, 0, Math.PI * 2);
  ctx.fillStyle = "#ffcf7a";
  ctx.fill();
}

// Sky compass: a view looking up, north at the top, east on the left.
// Labels are placed by testing candidate positions around each marker
// against every drawn element, so no two texts or markers overlap.
function boxOverlap(a, b) {
  const w = Math.min(a.x + a.w, b.x + b.w) - Math.max(a.x, b.x),
    h = Math.min(a.y + a.h, b.y + b.h) - Math.max(a.y, b.y);
  return w > 0 && h > 0 ? w * h : 0;
}
function textBox(node, fallbackWidth, size) {
  let b;
  try {
    b = node.getBBox();
  } catch {
    b = null;
  }
  return b && b.width
    ? { w: b.width, h: b.height }
    : { w: fallbackWidth, h: size * 1.2 };
}
function drawSky() {
  const s = $("sky-plot");
  s.replaceChildren();
  const VW = 320,
    VH = 244,
    cx = 160,
    cy = 122,
    R = 78,
    ring = R + 14,
    cardinal = ring + 18,
    obstacles = [];
  s.setAttribute("viewBox", `0 0 ${VW} ${VH}`);
  s.append(
    el("circle", {
      cx,
      cy,
      r: ring,
      fill: "none",
      stroke: "#5b6e82",
      "stroke-dasharray": "2 4",
    }),
    el("circle", {
      cx,
      cy,
      r: R,
      fill: "#101a26",
      stroke: "#6b8095",
      "stroke-width": 1.2,
    }),
  );
  for (const r of [R / 3, (2 * R) / 3])
    s.append(
      el("circle", {
        cx,
        cy,
        r,
        fill: "none",
        stroke: "#3a4e63",
        "stroke-dasharray": "2 4",
      }),
    );
  s.append(
    el("line", { x1: cx - R, y1: cy, x2: cx + R, y2: cy, stroke: "#2b3d50" }),
    el("line", { x1: cx, y1: cy - R, x2: cx, y2: cy + R, stroke: "#2b3d50" }),
    el("circle", { cx, cy, r: 2.2, fill: "#b4c3d0" }),
  );
  obstacles.push({ x: cx - 4, y: cy - 4, w: 8, h: 8 });
  for (const [text, x, y] of [
    ["N", cx, cy - cardinal],
    ["S", cx, cy + cardinal],
    ["E", cx - cardinal, cy],
    ["W", cx + cardinal, cy],
  ]) {
    const t = el(
      "text",
      {
        x,
        y,
        fill: "#b4c3d0",
        "font-size": 12,
        "font-weight": 600,
        "text-anchor": "middle",
        "dominant-baseline": "central",
      },
      text,
    );
    s.append(t);
    const b = textBox(t, 9, 12);
    obstacles.push({
      x: x - b.w / 2 - 5,
      y: y - b.h / 2 - 4,
      w: b.w + 10,
      h: b.h + 8,
    });
  }
  const bodies = ["sun", "earth"].map((name) => {
    const az = (selected[name + "_azimuth_deg"] * Math.PI) / 180,
      e = selected[name + "_elevation_deg"],
      above = e >= 0,
      r = above ? (R * (90 - e)) / 90 : ring;
    return {
      name,
      label: name === "sun" ? "Sun" : "Earth",
      color: name === "sun" ? "#efb778" : "#73c7bb",
      above,
      x: cx - r * Math.sin(az),
      y: cy - r * Math.cos(az),
      size: above ? 6 : 5.5,
    };
  });
  for (const b of bodies) {
    if (b.above && b.name === "sun")
      s.append(
        el("circle", { cx: b.x, cy: b.y, r: 11, fill: b.color, opacity: 0.18 }),
      );
    s.append(
      el("circle", {
        cx: b.x,
        cy: b.y,
        r: b.size,
        fill: b.above ? b.color : "#172130",
        stroke: b.above ? "#0f1822" : b.color,
        "stroke-width": b.above ? 1.5 : 2,
        "data-body": b.name,
      }),
    );
    const m = b.size + 3;
    obstacles.push({ x: b.x - m, y: b.y - m, w: 2 * m, h: 2 * m });
  }
  const straddles = (box, radius) => {
    const nx = clamp(cx, box.x, box.x + box.w) - cx,
      ny = clamp(cy, box.y, box.y + box.h) - cy,
      far = Math.max(
        ...[
          [box.x, box.y],
          [box.x + box.w, box.y],
          [box.x, box.y + box.h],
          [box.x + box.w, box.y + box.h],
        ].map(([x, y]) => Math.hypot(x - cx, y - cy)),
      );
    return Math.hypot(nx, ny) < radius && far > radius;
  };
  for (const b of bodies) {
    const t = el(
      "text",
      {
        fill: b.color,
        "font-size": 12,
        "font-weight": 600,
        "text-anchor": "middle",
        "dominant-baseline": "central",
        stroke: "#172130",
        "stroke-width": 3,
        "stroke-linejoin": "round",
        "paint-order": "stroke",
      },
      b.label,
    );
    s.append(t);
    const size = textBox(t, b.label.length * 7, 12),
      outward =
        Math.hypot(b.x - cx, b.y - cy) > 1
          ? Math.atan2(b.y - cy, b.x - cx)
          : -Math.PI / 2;
    let best;
    [0, 30, -30, 60, -60, 90, -90, 120, -120, 150, -150, 180].forEach(
      (offset, rank) => {
        for (const [step, gap] of [3, 10].entries()) {
          const a = outward + (offset * Math.PI) / 180,
            ux = Math.cos(a),
            uy = Math.sin(a),
            d =
              b.size +
              gap +
              Math.abs(ux) * (size.w / 2) +
              Math.abs(uy) * (size.h / 2),
            box = {
              x: b.x + ux * d - size.w / 2,
              y: b.y + uy * d - size.h / 2,
              w: size.w,
              h: size.h,
            };
          let score = rank * 0.6 + step * 2;
          for (const o of obstacles) score += 40 * boxOverlap(box, o);
          const inside = boxOverlap(box, { x: 2, y: 2, w: VW - 4, h: VH - 4 });
          score += 60 * (box.w * box.h - inside);
          if (straddles(box, R)) score += 6;
          if (straddles(box, ring)) score += 6;
          if (!best || score < best.score) best = { score, box };
        }
      },
    );
    t.setAttribute("x", best.box.x + best.box.w / 2);
    t.setAttribute("y", best.box.y + best.box.h / 2);
    obstacles.push(best.box);
  }
}

// The month chart is drawn at the container's CSS pixel width, so axis text
// keeps its size on every screen.
const chartSize = { W: 900, H: 280 };
function sizeChart() {
  const W = Math.round($("chart-container").clientWidth);
  if (!W) return false;
  const H = clamp(Math.round(W * 0.32), 210, 300);
  if (W === chartSize.W && H === chartSize.H) return false;
  chartSize.W = W;
  chartSize.H = H;
  $("light-chart").setAttribute("viewBox", `0 0 ${W} ${H}`);
  return true;
}
function plotCoordinates() {
  const { W, H } = chartSize,
    left = W < 520 ? 46 : 56,
    right = 14,
    top = 12,
    bottom = 30,
    plotW = W - left - right,
    plotH = H - top - bottom;
  return {
    W,
    H,
    left,
    right,
    top,
    bottom,
    plotW,
    plotH,
    x: (t) =>
      left + ((t - calendar.start) / (calendar.end - calendar.start)) * plotW,
    y: (v) =>
      top + ((5 - Math.log10(Math.max(v, Number.MIN_VALUE))) / 9) * plotH,
    fraction: (svgX) => clamp((svgX - left) / plotW, 0, 1),
  };
}
function drawChart() {
  if (!calendar) return;
  sizeChart();
  const s = $("light-chart");
  s.replaceChildren();
  const p = plotCoordinates();
  const defs = el("defs"),
    clip = el("clipPath", { id: "plot-area" });
  clip.append(
    el("rect", { x: p.left, y: p.top, width: p.plotW, height: p.plotH }),
  );
  defs.append(clip);
  s.append(defs);
  for (const exponent of [-4, -3, -2, -1, 0, 1, 2, 3, 4, 5]) {
    const v = 10 ** exponent,
      y = p.y(v);
    s.append(
      el("line", {
        x1: p.left,
        y1: y,
        x2: p.W - p.right,
        y2: y,
        stroke: exponent === 0 ? "#4a5d71" : "#2f3f51",
      }),
    );
    const label = exponent >= 3 ? `${v / 1000}k` : String(v);
    s.append(
      el(
        "text",
        {
          x: p.left - 8,
          y: y,
          fill: "#a3b4c4",
          "font-size": 11,
          "text-anchor": "end",
          "dominant-baseline": "central",
        },
        label,
      ),
    );
  }
  const days = Math.round((calendar.end - calendar.start) / DAY_MS),
    step = [5, 10, 15].find((d) => (p.plotW / days) * d >= 48) ?? 15;
  for (let d = 0; d <= days; d += step) {
    const x = p.x(calendar.start + d * DAY_MS);
    s.append(
      el("line", {
        x1: x,
        x2: x,
        y1: p.top + p.plotH,
        y2: p.top + p.plotH + 4,
        stroke: "#4a5d71",
      }),
      el(
        "text",
        {
          x,
          y: p.H - 8,
          fill: "#a3b4c4",
          "font-size": 11,
          "text-anchor": d === 0 ? "start" : x + 22 > p.W ? "end" : "middle",
        },
        formatDate(calendar.start + d * DAY_MS, {
          month: "short",
          day: "numeric",
        }),
      ),
    );
  }
  for (const [key, color, width] of [
    ["solar", "#efb778", 1.9],
    ["earth", "#73c7bb", 1.7],
    ["total", "#eaf0f3", 2.2],
  ]) {
    if (key === "earth" && !settings.earth) continue;
    const path = calendar.series
      .map(
        (r, i) =>
          `${i ? "L" : "M"}${p.x(r.ms).toFixed(2)},${p.y(r[key]).toFixed(2)}`,
      )
      .join(" ");
    s.append(
      el("path", {
        d: path,
        fill: "none",
        stroke: color,
        "stroke-width": width,
        "stroke-linejoin": "round",
        "stroke-linecap": "round",
        "clip-path": "url(#plot-area)",
      }),
    );
  }
  // Eclipse reference remains visibly marked instead of implying calculated shadow transport.
  for (const r of calendar.series)
    if (r.eclipse)
      s.append(
        el("line", {
          x1: p.x(r.ms),
          x2: p.x(r.ms),
          y1: p.top,
          y2: p.top + p.plotH,
          stroke: "#d79bb6",
          "stroke-opacity": 0.25,
          "stroke-width": 1.2,
        }),
      );
  s.append(el("g", { id: "chart-marker", "clip-path": "url(#plot-area)" }));
  s.append(el("g", { id: "chart-hover", "pointer-events": "none" }));
  drawMarker();
  $("range-start").textContent = formatDate(calendar.start, {
    month: "short",
    day: "numeric",
  });
  $("range-end").textContent = formatDate(calendar.end - 1, {
    month: "short",
    day: "numeric",
  });
  $("timeline-title").textContent = formatDate(calendar.start, {
    month: "long",
    year: "numeric",
  });
}
function drawMarker() {
  if (!calendar) return;
  const group = $("chart-marker");
  if (!group) return;
  group.replaceChildren();
  const p = plotCoordinates(),
    x = p.x(ms);
  $("scrubber").value = clamp(
    ((ms - calendar.start) / (calendar.end - calendar.start)) * 1000,
    0,
    1000,
  );
  if (ms < calendar.start || ms > calendar.end) return;
  group.append(
    el("line", {
      x1: x,
      x2: x,
      y1: p.top,
      y2: p.top + p.plotH,
      stroke: "#dbe7ee",
      "stroke-opacity": 0.55,
      "stroke-dasharray": "3 4",
    }),
  );
  group.append(
    el("circle", {
      cx: x,
      cy: p.y(selected.total),
      r: 4.5,
      fill: "#eaf0f3",
      stroke: "#142130",
      "stroke-width": 2,
    }),
  );
}
function chartRow(event) {
  const p = plotCoordinates(),
    rect = $("light-chart").getBoundingClientRect(),
    f = p.fraction(((event.clientX - rect.left) / rect.width) * p.W);
  return { p, rect, f };
}
function showHover(event) {
  if (!calendar) return;
  const { p, rect, f } = chartRow(event),
    i = Math.min(
      calendar.series.length - 1,
      Math.round(f * (calendar.series.length - 1)),
    ),
    r = calendar.series[i],
    tip = $("chart-tooltip"),
    hover = $("chart-hover"),
    hx = p.x(r.ms),
    hy = clamp(p.y(r.total), p.top, p.top + p.plotH);
  if (hover) {
    hover.replaceChildren(
      el("line", {
        x1: hx,
        x2: hx,
        y1: p.top,
        y2: p.top + p.plotH,
        stroke: "#dbe7ee",
        "stroke-opacity": 0.3,
      }),
      el("circle", {
        cx: hx,
        cy: hy,
        r: 4,
        fill: "#10161f",
        stroke: "#eaf0f3",
        "stroke-width": 2,
      }),
    );
  }
  tip.replaceChildren();
  const strong = document.createElement("strong");
  strong.textContent = lux(r.total);
  tip.append(
    strong,
    `${formatDate(r.ms, { month: "short", day: "numeric", hour: "2-digit", minute: "2-digit", hour12: false })} ${settings.scale}`,
  );
  tip.hidden = false;
  // Keep the tooltip beside the point and inside the chart on every side.
  const scale = rect.width / p.W,
    px = hx * scale,
    py = hy * scale,
    tw = tip.offsetWidth,
    th = tip.offsetHeight,
    cw = $("chart-container").clientWidth,
    ch = rect.height;
  let left = px + 14,
    top = py - th - 12;
  if (left + tw > cw) left = px - 14 - tw;
  if (top < 0) top = py + 12;
  tip.style.left = `${clamp(left, 0, Math.max(0, cw - tw))}px`;
  tip.style.top = `${clamp(top, 0, Math.max(0, ch - th))}px`;
}
function hideHover() {
  $("chart-tooltip").hidden = true;
  $("chart-hover")?.replaceChildren();
}

function colour(v) {
  const t = clamp((Math.log10(Math.max(v, 1e-5)) + 3) / 8, 0, 1),
    stops = [
      [24, 32, 53],
      [55, 70, 101],
      [140, 118, 136],
      [223, 158, 101],
      [243, 213, 149],
    ],
    x = t * 4,
    i = Math.min(3, Math.floor(x)),
    f = x - i;
  return `rgb(${stops[i].map((a, k) => Math.round(a * (1 - f) + stops[i + 1][k] * f)).join(",")})`;
}
function renderDaily() {
  const body = $("daily-table");
  body.replaceChildren();
  for (const d of calendar.daily) {
    const tr = document.createElement("tr");
    if (new Date(d.ms).toISOString().slice(0, 10) === $("date").value)
      tr.className = "selected";
    const dateCell = document.createElement("td"),
      button = document.createElement("button");
    button.textContent = formatDate(d.ms, { month: "short", day: "numeric" });
    button.setAttribute(
      "aria-label",
      `View ${formatDate(d.ms, { month: "long", day: "numeric", year: "numeric" })}`,
    );
    dateCell.append(button);
    const small = document.createElement("span");
    small.className = "day-weekday";
    small.textContent = formatDate(d.ms, { weekday: "short" });
    dateCell.append(small);
    tr.append(dateCell);
    const stripCell = document.createElement("td"),
      strip = document.createElement("div");
    strip.className = "day-strip";
    strip.setAttribute("role", "img");
    strip.setAttribute(
      "aria-label",
      `Light ranges from ${lux(d.min)} to ${lux(d.max)}`,
    );
    for (const v of d.strip) {
      const bar = document.createElement("i");
      bar.style.background = colour(v);
      strip.append(bar);
    }
    stripCell.append(strip);
    tr.append(stripCell);
    for (const text of [
      lux(d.min) + (d.eclipse ? " *" : ""),
      lux(d.max),
      `${number(d.sunHours, 1)} h`,
      settings.earth ? lux(d.earth) : "Off",
    ]) {
      const td = document.createElement("td");
      td.textContent = text;
      tr.append(td);
    }
    tr.addEventListener("click", () => setMoment(d.ms + 12 * 3600000, false));
    body.append(tr);
  }
}
function renderEvents() {
  if (!calendar) return;
  const list = $("events-list");
  list.replaceChildren();
  const events = calendar.events.filter((e) => e.ms >= ms).slice(0, 4);
  if (!events.length) {
    const p = document.createElement("p");
    p.className = "empty-events";
    p.textContent =
      "No further crossings in this month. Move to the next month to continue the calendar.";
    list.append(p);
    return;
  }
  for (const e of events) {
    const row = document.createElement("div");
    row.className = "event-row" + (e.key === "earth" ? " earth" : "");
    row.append(document.createElement("i"));
    const description = document.createElement("span");
    description.textContent = e.label;
    if (e.level === 2.98) {
      const small = document.createElement("small");
      small.textContent = "Practical twilight threshold";
      description.append(small);
    }
    row.append(description);
    const time = document.createElement("strong");
    time.textContent = formatDate(e.ms, { month: "short", day: "numeric" });
    const small = document.createElement("small");
    small.textContent = formatDate(e.ms, {
      hour: "2-digit",
      minute: "2-digit",
      hour12: false,
    });
    time.append(small);
    row.append(time);
    list.append(row);
  }
}

function calculate() {
  $("chart-status").textContent = "Calculating…";
  $("chart-container").classList.add("computing");
  $("download").disabled = true;
  worker.postMessage({ id: ++request, settings });
}
function setMoment(value, recompute = true) {
  ms = clamp(value, minMs, maxMs);
  dateInput();
  const d = new Date(ms),
    changed =
      d.getUTCFullYear() !== settings.year ||
      d.getUTCMonth() !== settings.month;
  if (recompute && changed) {
    settings.year = d.getUTCFullYear();
    settings.month = d.getUTCMonth();
    sync();
    calculate();
  } else updateQuery();
  updateSnapshot();
  if (calendar) renderDaily();
}
function changeLocation(lon, lat) {
  settings.longitude = Math.round(lon * 10000) / 10000;
  settings.latitude = Math.round(lat * 10000) / 10000;
  sync();
  updateSnapshot();
  calculate();
}
function updateInputs() {
  if (!$("latitude").checkValidity() || !$("longitude").checkValidity()) return;
  changeLocation(Number($("longitude").value), Number($("latitude").value));
}
function updateTime() {
  if (!$("date").checkValidity() || !$("time").checkValidity()) return;
  const time = $("time").value || "00:00";
  const value = Date.parse(
    `${$("date").value}T${time.length === 5 ? time + ":00" : time}Z`,
  );
  if (Number.isFinite(value) && value >= minMs && value <= maxMs)
    setMoment(value);
}
function monthTo(year, month) {
  const d = new Date(Date.UTC(year, month, 1));
  if (d.getUTCFullYear() < 2000 || d.getUTCFullYear() > 2500) return;
  settings.year = d.getUTCFullYear();
  settings.month = d.getUTCMonth();
  ms = d.valueOf() + 12 * 3600000;
  sync();
  updateSnapshot();
  calculate();
}

function bind() {
  $("latitude").addEventListener("change", updateInputs);
  $("longitude").addEventListener("change", updateInputs);
  $("place").addEventListener("change", () => {
    const site = astronomy.sites.find((s) => s.id === $("place").value);
    if (site) changeLocation(site.longitude, site.latitude);
  });
  $("date").addEventListener("change", updateTime);
  $("time").addEventListener("change", updateTime);
  $("timescale").addEventListener("change", () => {
    settings.scale = $("timescale").value;
    updateQuery();
    updateSnapshot();
    calculate();
  });
  $("earth-toggle").addEventListener("change", () => {
    settings.earth = $("earth-toggle").checked;
    updateQuery();
    updateSnapshot();
    calculate();
  });
  $("now").addEventListener("click", () => {
    settings.scale = "UTC";
    $("timescale").value = "UTC";
    setMoment(Date.now());
    calculate();
  });
  $("month").addEventListener("change", () => {
    if (!$("month").value || !$("month").checkValidity()) return;
    const [y, m] = $("month").value.split("-").map(Number);
    monthTo(y, m - 1);
  });
  $("previous-month").addEventListener("click", () =>
    monthTo(settings.year, settings.month - 1),
  );
  $("next-month").addEventListener("click", () =>
    monthTo(settings.year, settings.month + 1),
  );
  $("scrubber").addEventListener("input", () => {
    if (calendar) {
      cancelAnimationFrame(drawFrame);
      drawFrame = requestAnimationFrame(() =>
        setMoment(
          calendar.start +
            (Number($("scrubber").value) / 1000) *
              (calendar.end - calendar.start - 60000),
          false,
        ),
      );
    }
  });
  $("world-map").addEventListener("click", (event) => {
    const r = $("world-map").getBoundingClientRect();
    changeLocation(
      clamp(((event.clientX - r.left) / r.width) * 360 - 180, -180, 180),
      clamp(90 - ((event.clientY - r.top) / r.height) * 180, -90, 90),
    );
  });
  $("world-map").addEventListener("keydown", (e) => {
    const keys = {
      ArrowLeft: [-1, 0],
      ArrowRight: [1, 0],
      ArrowUp: [0, 1],
      ArrowDown: [0, -1],
    };
    if (keys[e.key]) {
      e.preventDefault();
      const [x, y] = keys[e.key];
      changeLocation(
        clamp(settings.longitude + x, -180, 180),
        clamp(settings.latitude + y, -90, 90),
      );
    }
  });
  $("light-chart").addEventListener("pointermove", showHover);
  $("light-chart").addEventListener("pointerleave", hideHover);
  $("light-chart").addEventListener("click", (e) => {
    if (!calendar) return;
    setMoment(
      calendar.start + chartRow(e).f * (calendar.end - calendar.start - 60000),
      false,
    );
  });
  $("light-chart").addEventListener("keydown", (e) => {
    if (["ArrowLeft", "ArrowRight"].includes(e.key)) {
      e.preventDefault();
      setMoment(ms + (e.key === "ArrowLeft" ? -1 : 1) * 3600000);
    }
  });
  $("share").addEventListener("click", async () => {
    try {
      await navigator.clipboard.writeText(location.href);
      notice("Link copied");
    } catch {
      const input = document.createElement("input");
      input.value = location.href;
      document.body.append(input);
      input.select();
      const copied = document.execCommand("copy");
      input.remove();
      notice(
        copied ? "Link copied" : "Copy the current address to share this view",
      );
    }
  });
  $("download").addEventListener("click", () => {
    if (!calendar) return;
    const header = [
      "date_time",
      "time_scale",
      "latitude_deg",
      "longitude_deg",
      "solar_direct_lux",
      "solar_diffuse_lux",
      "earth_direct_lux",
      "earth_diffuse_lux",
      "total_lux",
      "sun_elevation_deg",
      "earth_elevation_deg",
      "earth_lit_fraction",
      "eclipse_reference",
    ];
    const rows = calendar.series.map((r) => [
      new Date(r.ms).toISOString().replace("Z", ""),
      settings.scale,
      settings.latitude,
      settings.longitude,
      r.solar_direct,
      r.solar_diffuse,
      r.earth_direct,
      r.earth_diffuse,
      r.total,
      r.sun_elevation_deg,
      r.earth_elevation_deg,
      r.earth_lit_fraction,
      r.eclipse,
    ]);
    const url = URL.createObjectURL(
      new Blob([[header, ...rows].map((r) => r.join(",")).join("\n")], {
        type: "text/csv;charset=utf-8",
      }),
    );
    const link = document.createElement("a");
    link.href = url;
    link.download = `terluna-light-${$("month").value}.csv`;
    link.click();
    setTimeout(() => URL.revokeObjectURL(url), 1000);
  });
}

async function start() {
  try {
    [astronomy, transfer, validation] = await Promise.all(
      ["data/astronomy.json", "data/transfer.json", "data/validation.json"].map(
        async (url) => {
          const r = await fetch(url);
          if (!r.ok) throw Error(`Lighting data unavailable (${r.status})`);
          return r.json();
        },
      ),
    );
    if (
      astronomy.schema !== "terluna.illumination.calendar-astronomy/1" ||
      transfer.schema !== "terluna.illumination.calendar-transfer/1" ||
      validation.schema !== "terluna.illumination.calendar-validation/1"
    )
      throw Error("These model files have an unsupported version.");
    for (const site of astronomy.sites) {
      const o = document.createElement("option");
      o.value = site.id;
      o.textContent = site.name;
      $("place").append(o);
    }
    const custom = document.createElement("option");
    custom.value = "custom";
    custom.textContent = "Custom coordinates";
    $("place").append(custom);
    if (!query.has("lon") && !query.has("lat")) {
      const site =
        astronomy.sites.find((s) => s.name.includes("Nubium")) ??
        astronomy.sites[0];
      settings.longitude = site.longitude;
      settings.latitude = site.latitude;
    }
    worker = new Worker(new URL("./worker.mjs", import.meta.url), {
      type: "module",
    });
    worker.postMessage({ type: "init", astronomy, transfer });
    worker.onmessage = ({ data }) => {
      if (data.id !== request) return;
      if (data.error) {
        $("chart-status").textContent = "Calculation failed";
        $("error").hidden = false;
        $("error").textContent = data.error;
        return;
      }
      calendar = data.result;
      $("error").hidden = true;
      drawChart();
      renderDaily();
      renderEvents();
      $("chart-container").classList.remove("computing");
      $("chart-status").textContent = "30-minute samples · drag to explore";
      $("download").disabled = false;
    };
    worker.onerror = () => {
      $("chart-status").textContent = "Calendar calculation unavailable";
      $("error").hidden = false;
      $("error").textContent =
        "The calculation worker could not start. Reload the page or check that all app files are served from the same origin.";
    };
    const orbit = astronomy.ephemeris_validation;
    const comparison =
      validation.ephemeris_to_light.directions_and_geometric_ranges;
    const relative = comparison.components.total.relative.find(
      (row) => row.reference_at_least_lux === 0.1,
    );
    const accuracy = [
      `Checked against JPL Horizons at ${number(orbit.unique_epochs)} unique epochs across ${astronomy.calendar_range.join("–")}. The largest sampled direction differences are ${number(orbit.errors.sun_direction_deg.max_abs, 4)}° for the Sun and ${number(orbit.errors.earth_direction_deg.max_abs, 4)}° for Earth. Tested Sun-centre horizon crossings differ by up to ${number(orbit.horizon_event_abs_error_max_minutes, 1)} minutes.`,
      `The finite-source calculation compares positions and distances at ${number(validation.geometric_range_samples)} epochs across ${validation.sites.length} locations. For reference totals of at least 0.1 lux, the largest sampled change in total light is ${number(relative.max_abs_percent, 2)}%. Eclipse cases are excluded. This measures the effect of the geometry difference within the same atmosphere and source model.`,
      "These are sampled comparisons, not uniform accuracy bounds. Near-grazing thresholds can shift substantially: a check using the earlier solar-only twilight curve found shifts of about an hour. Atmospheric conditions and Earth's changing appearance contribute separate uncertainty.",
    ];
    for (const text of accuracy) {
      const p = document.createElement("p");
      p.textContent = text;
      $("accuracy-text").append(p);
    }
    const ul = document.createElement("ul");
    for (const source of transfer.sources ?? []) {
      const li = document.createElement("li"),
        a = document.createElement("a");
      a.href = source.url;
      a.textContent = source.title;
      li.append(a);
      ul.append(li);
    }
    $("sources").append(ul);
    const provenance = document.createElement("a");
    provenance.href = "data/transfer.json";
    provenance.download = "terluna-light-transfer.json";
    provenance.textContent = "Download the versioned model product";
    $("sources").append(provenance);
    const report = document.createElement("p");
    const reportLink = document.createElement("a");
    reportLink.href = "data/validation.json";
    reportLink.download = "terluna-light-validation.json";
    reportLink.textContent = "Download the numerical validation report";
    report.append(reportLink);
    $("sources").append(report);
    $("build-info").textContent =
      `${transfer.schema} · ${transfer.checks?.scattering_orders ?? "—"} scattering orders`;
    $("load-status").hidden = true;
    $("application").hidden = false;
    $("share").disabled = false;
    bind();
    sync();
    updateSnapshot();
    calculate();
    if ("ResizeObserver" in window) {
      let mapWidth = 0;
      new ResizeObserver(() => {
        if (sizeChart() && calendar) drawChart();
      }).observe($("chart-container"));
      new ResizeObserver(([entry]) => {
        const width = Math.round(entry.contentRect.width);
        if (width && width !== mapWidth && selected) {
          mapWidth = width;
          drawMap();
        }
      }).observe($("world-map"));
    }
  } catch (error) {
    $("load-status").hidden = true;
    $("error").hidden = false;
    $("error").textContent = `The calendar could not load. ${error.message}`;
  }
}
start();
