// Open Moon Skies: the public explorer of the light calendar. The light and
// positions come from the shared calendar evaluator; this file arranges them.
import { stateAt, ephemeris, julianDate, interpolate, DAY_MS } from "./model/evaluator.mjs";
import { skyFrame, subObserverPoint } from "./sky-frame.mjs";
import {
  SkyAtlas,
  SkyView,
  sceneFrom,
  earthImage,
  ToneCurve,
  renderDome,
  prepareStars,
} from "./sky-render.mjs";
import { Globe } from "./globe.mjs";
import {
  REGIMES,
  EVENT_INFO,
  PRACTICAL_DUSK_LUX,
  lux,
  duration,
  shortDuration,
  relative,
  earthPhase,
  earthComparison,
  moonlightRatio,
  altitudeWords,
  direction,
  compass,
  fmtNumber,
} from "./words.mjs";

const $ = (id) => document.getElementById(id);
const HOUR = 3600000;
const SYNODIC_DEG_PER_HOUR = 360 / (29.530589 * 24);
const MOON_FROM_EARTH_DEG = 0.518; // mean apparent diameter of the Moon from Earth
const minMs = Date.UTC(2000, 0, 1),
  maxMs = Date.UTC(2501, 0, 1) - 60000;
const clamp = (x, a, b) => Math.max(a, Math.min(b, x));
const rad = Math.PI / 180;
const settings = { earth: true, scale: "UTC" };

const app = {
  ms: Date.now(),
  place: null,
  span: null,
  month: null,
  year: null,
  viewMonth: null,
  viewYear: null,
  playing: false,
  speed: 6,
  sectors: null,
  state: null,
};

// ---------- Formatting ----------
const fmtDate = (ms, opts) =>
  new Intl.DateTimeFormat("en-GB", { timeZone: "UTC", ...opts }).format(new Date(ms));
const longDate = (ms) =>
  fmtDate(ms, { weekday: "long", day: "numeric", month: "long", year: "numeric" });
const shortDate = (ms) => fmtDate(ms, { day: "numeric", month: "short" });
const clock = (ms) => fmtDate(ms, { hour: "2-digit", minute: "2-digit", hour12: false });
function svgEl(tag, attrs = {}, text) {
  const n = document.createElementNS("http://www.w3.org/2000/svg", tag);
  for (const [k, v] of Object.entries(attrs)) n.setAttribute(k, String(v));
  if (text !== undefined) n.textContent = text;
  return n;
}
function el(tag, attrs = {}, ...children) {
  const n = document.createElement(tag);
  for (const [k, v] of Object.entries(attrs)) {
    if (k === "class") n.className = v;
    else if (k === "style") n.style.cssText = v;
    else if (k.startsWith("on")) n.addEventListener(k.slice(2), v);
    else n.setAttribute(k, v);
  }
  for (const c of children) if (c !== null && c !== undefined) n.append(c);
  return n;
}

// ---------- Icons ----------
const ICONS = {
  sunrise: (c) =>
    `<svg viewBox="0 0 32 32"><path d="M5 22h22" stroke="${c}" stroke-width="2.2" stroke-linecap="round"/><path d="M9.5 22a6.5 6.5 0 0 1 13 0" fill="${c}"/><path d="M16 4v6M12.5 7.5 16 4l3.5 3.5" fill="none" stroke="${c}" stroke-width="2.2" stroke-linecap="round" stroke-linejoin="round"/><path d="M6 15.5l2 1.2M26 15.5l-2 1.2" stroke="${c}" stroke-width="2" stroke-linecap="round"/></svg>`,
  sunset: (c) =>
    `<svg viewBox="0 0 32 32"><path d="M5 22h22" stroke="${c}" stroke-width="2.2" stroke-linecap="round"/><path d="M9.5 22a6.5 6.5 0 0 1 13 0" fill="${c}"/><path d="M16 4v6M12.5 6.5 16 10l3.5-3.5" fill="none" stroke="${c}" stroke-width="2.2" stroke-linecap="round" stroke-linejoin="round"/><path d="M6 15.5l2 1.2M26 15.5l-2 1.2" stroke="${c}" stroke-width="2" stroke-linecap="round"/></svg>`,
  dusk: (c) =>
    `<svg viewBox="0 0 32 32"><path d="M5 24h22" stroke="${c}" stroke-width="2.2" stroke-linecap="round"/><path d="M8 19h16M11 14h10" stroke="${c}" stroke-width="2.2" stroke-linecap="round" opacity=".7"/><path d="M16 4v5.5M12.8 6.6 16 9.5l3.2-2.9" fill="none" stroke="${c}" stroke-width="2.2" stroke-linecap="round" stroke-linejoin="round"/></svg>`,
  dawn: (c) =>
    `<svg viewBox="0 0 32 32"><path d="M5 24h22" stroke="${c}" stroke-width="2.2" stroke-linecap="round"/><path d="M8 19h16M11 14h10" stroke="${c}" stroke-width="2.2" stroke-linecap="round" opacity=".7"/><path d="M16 4v5.5M12.8 7 16 4l3.2 3" fill="none" stroke="${c}" stroke-width="2.2" stroke-linecap="round" stroke-linejoin="round"/></svg>`,
  earthrise: (c) =>
    `<svg viewBox="0 0 32 32"><path d="M5 23h22" stroke="${c}" stroke-width="2.2" stroke-linecap="round"/><circle cx="16" cy="16" r="6" fill="${c}"/><path d="M12 15c2-1 3 1 5 0s2-2 3-1" stroke="#0b1026" stroke-width="1.4" fill="none" opacity=".6"/><path d="M16 3v4M13.5 5.5 16 3l2.5 2.5" fill="none" stroke="${c}" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"/></svg>`,
  earthset: (c) =>
    `<svg viewBox="0 0 32 32"><path d="M5 23h22" stroke="${c}" stroke-width="2.2" stroke-linecap="round"/><circle cx="16" cy="16" r="6" fill="${c}"/><path d="M16 3v4M13.5 4.5 16 7l2.5-2.5" fill="none" stroke="${c}" stroke-width="2" stroke-linecap="round" stroke-linejoin="round"/></svg>`,
  full: (c) =>
    `<svg viewBox="0 0 32 32"><circle cx="16" cy="16" r="9.5" fill="${c}"/><path d="M10 13c3-2 5 2 8 0s4-1 5 1M9 19c2 1 4-1 6 1" stroke="#0b1026" stroke-width="1.5" fill="none" opacity=".45"/></svg>`,
  new: (c) =>
    `<svg viewBox="0 0 32 32"><circle cx="16" cy="16" r="9.5" fill="#0b1026" stroke="${c}" stroke-width="2"/><path d="M16 6.5a9.5 9.5 0 0 1 0 19" fill="none" stroke="${c}" stroke-width="1.2" opacity=".5"/></svg>`,
  eclipse: (c) =>
    `<svg viewBox="0 0 32 32"><circle cx="16" cy="16" r="10" fill="${c}" opacity=".35"/><circle cx="16" cy="16" r="8.4" fill="#0b1026" stroke="${c}" stroke-width="2.2"/></svg>`,
};
const icon = (type) => {
  const info = EVENT_INFO[type];
  return ICONS[info.icon](info.colour);
};

// ---------- Data ----------
async function json(url) {
  const r = await fetch(url);
  if (!r.ok) throw Error(`${url} unavailable (${r.status})`);
  return r.json();
}
async function binary(url) {
  const r = await fetch(url);
  if (!r.ok) throw Error(`${url} unavailable (${r.status})`);
  return r.arrayBuffer();
}
async function imageData(url) {
  const img = new Image();
  img.decoding = "async";
  img.src = url;
  await img.decode();
  const c = document.createElement("canvas");
  c.width = img.naturalWidth;
  c.height = img.naturalHeight;
  const ctx = c.getContext("2d", { willReadFrequently: true });
  ctx.drawImage(img, 0, 0);
  return ctx.getImageData(0, 0, c.width, c.height);
}

let astronomy, transfer, assets, moonAtlas, earthAtlas, stars, earthTexture, moonTexture, heightMap, places, build;
let worker,
  requestId = 0;
const pending = new Map();
function ask(type, payload) {
  const id = ++requestId;
  worker.postMessage({ id, type, ...payload });
  return new Promise((resolve, reject) => pending.set(id, { resolve, reject }));
}

// ---------- Places and the ground ----------
function heightAt(lon, lat) {
  const w = heightMap.width,
    h = heightMap.height;
  const x = Math.min(w - 1, Math.max(0, Math.floor(((((lon + 180) % 360) + 360) % 360) / 360 * w)));
  const y = Math.min(h - 1, Math.max(0, Math.floor(((90 - lat) / 180) * h)));
  const code = heightMap.data[(y * w + x) * 4];
  const s = (code - 128) / 127;
  return Math.sign(s) * s * s * 10000;
}
function destination(lon, lat, azDeg, km) {
  const d = km / 1737.4,
    p = lat * rad,
    t = azDeg * rad;
  const lat2 = Math.asin(Math.sin(p) * Math.cos(d) + Math.cos(p) * Math.sin(d) * Math.cos(t));
  const lon2 =
    lon * rad +
    Math.atan2(Math.sin(t) * Math.sin(d) * Math.cos(p), Math.cos(d) - Math.sin(p) * Math.sin(lat2));
  return [lon2 / rad, lat2 / rad];
}
// Sea or land toward each 5° of azimuth, from the atlas around the place.
function sectorsFor(lon, lat) {
  const s = new Uint8Array(72);
  for (let k = 0; k < 72; k++) {
    let wet = 0;
    for (const km of [6, 14, 28]) {
      const [a, b] = destination(lon, lat, k * 5, km);
      if (heightAt(a, b) < 0) wet++;
    }
    s[k] = wet >= 2 ? 1 : 0;
  }
  const out = s.slice();
  for (let k = 0; k < 72; k++) {
    const l = s[(k + 71) % 72],
      r = s[(k + 1) % 72];
    if (l === r && s[k] !== l) out[k] = l;
  }
  return out;
}
function setting(place) {
  if (place.status === "submerged" && place.depth_m)
    return `Under ${fmtNumber(place.depth_m)} m of ${place.water_body ?? "sea"}`;
  const h = heightAt(place.longitude, place.latitude);
  if (h < 0) return `Sea, about ${fmtNumber(Math.round(-h / 10) * 10)} m deep`;
  return `Land, about ${fmtNumber(Math.round(h / 10) * 10)} m above sea level`;
}
function coord(lat, lon) {
  return `${fmtNumber(Math.abs(lat), 1)}° ${lat < 0 ? "S" : "N"}, ${fmtNumber(Math.abs(lon), 1)}° ${lon < 0 ? "W" : "E"}`;
}
const GROUP_COLOURS = {
  "Study coasts": "#ffc53d",
  "Landing sites and heritage": "#ff5d8f",
  Seas: "#38c6f4",
  "Island summits": "#9be37a",
};
function allPlaces() {
  return places.groups.flat();
}
function nearestNamed(lon, lat) {
  let best = null,
    bd = Infinity;
  for (const p of allPlaces()) {
    const d =
      Math.acos(
        clamp(
          Math.sin(lat * rad) * Math.sin(p.latitude * rad) +
            Math.cos(lat * rad) * Math.cos(p.latitude * rad) * Math.cos((lon - p.longitude) * rad),
          -1,
          1,
        ),
      ) * 1737.4;
    if (d < bd) {
      bd = d;
      best = p;
    }
  }
  return [best, bd];
}
function customPlace(lon, lat) {
  const [near, d] = nearestNamed(lon, lat);
  const h = heightAt(lon, lat);
  return {
    id: "custom",
    custom: true,
    name: d < 120 ? `Near ${near.name}` : coord(lat, lon),
    group: "Your point",
    longitude: Math.round(lon * 1000) / 1000,
    latitude: Math.round(lat * 1000) / 1000,
    description: h < 0 ? "A point on the open sea" : "A point on land",
  };
}

// ---------- Time ----------
function setTime(ms, opts = {}) {
  app.ms = clamp(ms, minMs, maxMs);
  scheduleUpdate(opts);
}
let frameRequested = false,
  lastHeavy = 0,
  heavyForce = false;
function scheduleUpdate(opts = {}) {
  if (opts.heavy) heavyForce = true;
  if (frameRequested) return;
  frameRequested = true;
  requestAnimationFrame(() => {
    frameRequested = false;
    update();
  });
}
function update() {
  const ms = app.ms,
    p = app.place;
  const s = stateAt(ms, p.longitude, p.latitude, astronomy, transfer, settings.earth, settings.scale);
  app.state = s;
  const frame = skyFrame(ms, settings.scale, p.longitude, p.latitude, astronomy);
  app.frame = frame;
  app.scene = sceneFrom(s, frame, moonAtlas, app.sectors);
  skyView.render(app.scene);
  drawSkyLabels();
  drawTimeline();
  const now = performance.now();
  const heavyDue = heavyForce || !app.playing || now - lastHeavy > 350;
  if (heavyDue) {
    heavyForce = false;
    lastHeavy = now;
    updateText();
    updateCards();
    updateGlobeLight();
    updateDial();
    syncInputs();
    ensureData();
  }
}

function ensureData() {
  const ms = app.ms;
  if (!app.span || ms < app.span.start + 16 * DAY_MS || ms > app.span.end - 32 * DAY_MS)
    requestSpan();
  const d = new Date(ms);
  const key = `${d.getUTCFullYear()}-${d.getUTCMonth()}`;
  if (!app.playing || !app.month) {
    if (app.monthKey !== key && !app.monthPending) {
      app.viewMonth = [d.getUTCFullYear(), d.getUTCMonth()];
      requestMonth();
    }
  }
  if (!app.playing && app.viewYear !== d.getUTCFullYear() && !app.yearLocked) {
    app.viewYear = d.getUTCFullYear();
    requestYear();
  }
  queryUpdate();
}

let spanToken = 0;
async function requestSpan() {
  if (app.spanPending) return;
  app.spanPending = true;
  const token = ++spanToken;
  const center = app.ms;
  const p = app.place;
  try {
    const result = await ask("span", {
      settings: { ...settings, longitude: p.longitude, latitude: p.latitude },
      start: Math.max(minMs, center - 32 * DAY_MS),
      end: Math.min(maxMs, center + 62 * DAY_MS),
      step: HOUR,
      order: 8,
    });
    if (token === spanToken && p === app.place) {
      app.span = result;
      scheduleUpdate({ heavy: true });
    }
  } finally {
    app.spanPending = false;
    if (token !== spanToken) requestSpan();
  }
}
async function requestMonth() {
  const [y, m] = app.viewMonth;
  const p = app.place;
  app.monthPending = true;
  try {
    const result = await ask("month", {
      settings: { ...settings, longitude: p.longitude, latitude: p.latitude, year: y, month: m },
    });
    if (p === app.place && app.viewMonth[0] === y && app.viewMonth[1] === m) {
      app.month = result;
      app.monthKey = `${y}-${m}`;
      renderMonth();
    }
  } catch (e) {
    showError(e.message);
  } finally {
    app.monthPending = false;
  }
}
async function requestYear() {
  const y = app.viewYear,
    p = app.place;
  $("year-label").textContent = String(y);
  try {
    const result = await ask("year", {
      settings: { ...settings, longitude: p.longitude, latitude: p.latitude, year: y },
    });
    if (p === app.place && app.viewYear === y) {
      app.year = result;
      renderYear();
    }
  } catch (e) {
    showError(e.message);
  }
}

// ---------- Hero ----------
let skyView;
// Frame the view: the horizon low in the picture and the chosen body inside it,
// to the right of the text on wide screens.
function aimCamera(target) {
  const s = app.state,
    cam = skyView.camera;
  const r = skyView.canvas.getBoundingClientRect();
  const W = r.width,
    H = r.height,
    portrait = W / H < 1.05;
  const earthUp = s.earth_elevation_deg > -0.5;
  const sunUp = s.sun_elevation_deg > -0.3;
  let kind = target;
  if (!kind) kind = earthUp && (!sunUp || s.sun_elevation_deg > 50) ? "earth" : "sun";
  let tEl, tAz;
  if (kind === "earth") {
    tEl = s.earth_elevation_deg;
    tAz = s.earth_azimuth_deg;
  } else {
    tEl = sunUp ? s.sun_elevation_deg : 6;
    tAz = s.sun_azimuth_deg;
  }
  const hy = portrait ? 0.5 : 0.72,
    ty = portrait ? 0.16 : 0.22,
    tx = portrait ? 0.5 : 0.68;
  const span = clamp(tEl, 32, 84);
  const px = clamp(((hy - ty) * H) / span, W / 220, W / 60);
  cam.hfov = W / px;
  cam.pitch = Math.min(((hy - 0.5) * H) / px, 90 - H / 2 / px);
  cam.yaw = tAz - ((tx - 0.5) * W) / px;
}
function drawSkyLabels() {
  const box = $("sky-labels");
  box.replaceChildren();
  const W = skyView.canvas.width,
    H = skyView.canvas.height,
    d = skyView.dpr;
  const cam = skyView.camera,
    s = app.state;
  const cssW = W / d,
    cssH = H / d,
    hY = cam.horizonY(H, W) / d;
  const flat = (el, az) => [Math.cos(el * rad) * Math.sin(az * rad), Math.cos(el * rad) * Math.cos(az * rad), Math.sin(el * rad)];
  if (hY > 0 && hY < cssH - 8)
    for (let k = 0; k < 8; k++) {
      const p = cam.project(flat(0, k * 45), W, H);
      if (!p || p[0] < 24 || p[0] > W - 24) continue;
      box.append(el("span", { class: "sky-label compass", style: `left:${p[0] / d}px;top:${hY}px` }, compass(k * 45)));
    }
  const bodies = [
    ["Sun", s.sun_elevation_deg, s.sun_azimuth_deg, s.sun_radius_deg, "var(--sun)", "sun"],
    ["Earth", s.earth_elevation_deg, s.earth_azimuth_deg, s.earth_radius_deg, "var(--earth)", "earth"],
  ];
  let lens = null;
  const chips = $("sky-chips");
  chips.replaceChildren();
  for (const [name, elv, az, r, colour, key] of bodies) {
    const up = elv + r > 0;
    const sub = up ? `${fmtNumber(elv, 0)}° up` : `${fmtNumber(-elv, 0)}° below the horizon`;
    const p = cam.project(flat(Math.max(elv, 0), az), W, H);
    const x = p ? p[0] / d : null,
      y = p ? p[1] / d : null;
    const offscreen = !p || x < 30 || x > cssW - 30 || y < 70;
    if (offscreen) {
      // A button that turns the view toward the body.
      const rel = ((((az - cam.yaw) % 360) + 540) % 360) - 180;
      const arrow = !p || x < 30 || x > cssW - 30 ? (rel > 0 ? "▶" : "◀") : "▲";
      const b = el("button", { type: "button", class: "sky-chip", "aria-label": `Look toward the ${name}` }, el("i", { class: `dot ${key}` }), `${name} · ${sub}`, el("span", { class: "arrow" }, arrow));
      b.addEventListener("click", () => {
        aimCamera(key);
        scheduleUpdate();
      });
      chips.append(b);
    } else if (!up) {
      if (hY < cssH - 30)
        box.append(el("span", { class: "sky-label", style: `left:${x}px;top:${hY + 18}px` }, el("i", { class: `dot ${key}` }), ` ${name} ↓`, el("small", {}, sub)));
    } else {
      const radius = (r * cam.pixelsPerDegree(W)) / d;
      box.append(el("span", { class: "sky-label", style: `left:${x}px;top:${y + Math.max(9, radius + 6)}px` }, name, el("small", {}, sub)));
      if (key === "earth" && radius < 22) lens = { x, y, radius };
    }
  }
  drawLens(lens, cssW, cssH);
}

// A close-up of Earth beside its true-size disk.
function drawLens(lens, cssW, cssH) {
  const c = $("earth-lens"),
    line = $("lens-line");
  if (!lens) {
    c.hidden = true;
    line.setAttribute("hidden", "");
    return;
  }
  const size = cssW < 700 ? 92 : 124;
  let lx = lens.x + 70,
    ly = lens.y - size - 26;
  if (lx + size > cssW - 16) lx = lens.x - 70 - size;
  if (ly < 70) ly = lens.y + 40;
  c.hidden = false;
  line.removeAttribute("hidden");
  c.style.left = `${lx}px`;
  c.style.top = `${ly}px`;
  c.style.width = c.style.height = `${size}px`;
  const n = Math.round(size * Math.min(window.devicePixelRatio || 1, 2));
  if (c.width !== n) c.width = c.height = n;
  paintEarth(c, { exposure: 4, tone: lensTone });
  const cx = lx + size / 2,
    cy = ly + size / 2;
  const dx = lens.x - cx,
    dy = lens.y - cy,
    dist = Math.hypot(dx, dy);
  line.setAttribute("viewBox", `0 0 ${cssW} ${cssH}`);
  line.replaceChildren(
    svgEl("line", {
      x1: cx + (dx / dist) * (size / 2 + 2),
      y1: cy + (dy / dist) * (size / 2 + 2),
      x2: lens.x - (dx / dist) * (lens.radius + 4),
      y2: lens.y - (dy / dist) * (lens.radius + 4),
      stroke: "rgba(255,255,255,0.6)",
      "stroke-width": 1.2,
      "stroke-dasharray": "3 3",
    }),
    svgEl("text", { x: cx, y: ly + size + 16, "text-anchor": "middle", fill: "rgba(255,255,255,0.85)", "font-size": 11.5, "font-weight": 600 }, `Earth, enlarged ${Math.round(size / (2 * lens.radius))}×`),
  );
}
const lensTone = new ToneCurve(1);
function paintEarth(c, { exposure, tone }) {
  const n = c.width;
  tone.sigma = 0.55;
  tone.saturation = 1;
  tone.night = 0;
  const img = earthImage(app.scene, earthTexture, n, tone, { exposure, noTint: true });
  const ctx = c.getContext("2d");
  ctx.clearRect(0, 0, n, n);
  ctx.fillStyle = "rgba(14,22,46,0.92)";
  ctx.beginPath();
  ctx.arc(n / 2, n / 2, n / 2 - 1, 0, Math.PI * 2);
  ctx.fill();
  const tmp = document.createElement("canvas");
  tmp.width = tmp.height = n;
  tmp.getContext("2d").putImageData(new ImageData(img.data, n, n), 0, 0);
  ctx.drawImage(tmp, 0, 0);
}

function bindSky() {
  const c = skyView.canvas;
  let drag = null;
  c.addEventListener("pointerdown", (e) => {
    c.setPointerCapture(e.pointerId);
    drag = { x: e.clientX, y: e.clientY, yaw: skyView.camera.yaw, pitch: skyView.camera.pitch };
    c.classList.add("dragging");
    $("drag-hint").classList.add("gone");
  });
  c.addEventListener("pointermove", (e) => {
    if (!drag) return;
    const r = c.getBoundingClientRect();
    const cam = skyView.camera;
    const k = cam.hfov / r.width;
    cam.yaw = drag.yaw - (e.clientX - drag.x) * k;
    cam.pitch = clamp(drag.pitch + (e.clientY - drag.y) * k, -20, 90 - (r.height / 2) * k);
    scheduleUpdate();
  });
  const end = () => {
    drag = null;
    c.classList.remove("dragging");
  };
  c.addEventListener("pointerup", end);
  c.addEventListener("pointercancel", end);
  c.addEventListener(
    "wheel",
    (e) => {
      if (!e.ctrlKey && Math.abs(e.deltaY) < 1) return;
      e.preventDefault();
      zoom(Math.exp(e.deltaY * 0.0015));
    },
    { passive: false },
  );
  $("zoom-in").addEventListener("click", () => zoom(0.7));
  $("zoom-out").addEventListener("click", () => zoom(1 / 0.7));
  $("look-earth").addEventListener("click", () => {
    aimCamera("earth");
    scheduleUpdate();
  });
  $("look-sun").addEventListener("click", () => {
    aimCamera("sun");
    scheduleUpdate();
  });
  new ResizeObserver(() => scheduleUpdate({ heavy: true })).observe(c);
}
function zoom(f) {
  const cam = skyView.camera;
  cam.hfov = clamp(cam.hfov * f, 12, 240);
  scheduleUpdate();
}

function regimeOf(s) {
  if (s.sun_elevation_deg + s.sun_radius_deg > 0) return 0;
  if (s.earth > s.solar && s.total >= 0.05) return 2;
  if (s.total >= PRACTICAL_DUSK_LUX) return 1;
  if (s.total >= 0.05) return 3;
  return 4;
}
function eventsAround(type) {
  const ev = app.span?.events ?? [];
  const before = ev.filter((e) => e.type === type && e.ms <= app.ms).at(-1);
  const after = ev.find((e) => e.type === type && e.ms > app.ms);
  return [before, after];
}
function earthReach() {
  const el = app.span?.samples.earthEl;
  if (!el) return null;
  const lo = Math.min(...el),
    hi = Math.max(...el);
  if (lo > 0.5) return "always";
  if (hi < -1) return "never";
  return "sometimes";
}
function isWaxing() {
  const s = app.state,
    p = app.place;
  const later = stateAt(app.ms + 6 * HOUR, p.longitude, p.latitude, astronomy, transfer, false, settings.scale);
  return later.earth_lit_fraction > s.earth_lit_fraction;
}

function updateText() {
  const s = app.state,
    p = app.place;
  const code = regimeOf(s),
    regime = REGIMES[code];
  $("regime-pill").style.setProperty("--regime", regime.colour);
  $("regime-name").textContent = regime.name;
  $("hero-date").textContent = `${longDate(app.ms)} · ${clock(app.ms)} UTC`;
  const [, nextRise] = eventsAround("sunrise");
  const [lastSet, nextSet] = eventsAround("sunset");
  const [, nextDusk] = eventsAround("dusk");
  const phase = earthPhase(s.earth_lit_fraction, isWaxing());
  const reach = earthReach();
  let headline = "",
    sub = regime.blurb;
  const rel = (e) => (e ? relative(e.ms - app.ms) : "beyond the next two months");
  if (code === 0) {
    headline = `The Sun is ${altitudeWords(s.sun_elevation_deg, s.sun_azimuth_deg)}.`;
    sub = nextSet ? `It sets ${rel(nextSet)}. Daytime here lasts about two Earth weeks.` : "The Sun stays up for weeks at a time here.";
  } else if (code === 1) {
    headline = lastSet ? `The Sun set ${relative(lastSet.ms - app.ms)}.` : "The Sun is below the horizon.";
    sub = `Sunlight scattered high in the deep air still lights the ground. ${nextDusk ? `Twilight stays bright enough to work by for another ${duration(nextDusk.ms - app.ms)}.` : ""}`;
  } else if (code === 2) {
    headline = `Earth lights the night: ${phase.toLowerCase().replace(" earth", "")}, ${s.earth_elevation_deg >= 0 ? `${fmtNumber(s.earth_elevation_deg, 0)}° up` : "on the horizon"}.`;
    const ratio = moonlightRatio(s.total);
    sub = `${ratio ? `The ground gets ${ratio}. ` : ""}Sunrise comes ${rel(nextRise)}.`;
  } else if (code === 3) {
    headline = `Night. The Sun is ${fmtNumber(-s.sun_elevation_deg, 0)}° below the horizon.`;
    sub = `A faint twilight glow remains. Sunrise comes ${rel(nextRise)}.`;
  } else {
    headline = "Deep night under the stars.";
    sub = `${reach === "never" ? "Earth never rises over this side of the Moon. " : ""}Sunrise comes ${rel(nextRise)}.`;
  }
  $("headline").textContent = headline;
  $("subline").textContent = sub;
  $("lux-now").textContent = lux(s.total);
  $("lux-compare").textContent = earthComparison(s.total, assets.sky.earth);
  // Next notable event.
  const next = (app.span?.events ?? []).find(
    (e) => e.ms > app.ms && ["sunrise", "sunset", "dusk", "dawn", "full-earth", "new-earth", "eclipse-begins", "earthrise", "earthset"].includes(e.type),
  );
  if (next) {
    $("next-label").textContent = EVENT_INFO[next.type].name;
    $("next-value").textContent = `in ${shortDuration(next.ms - app.ms)}`;
    $("next-note").textContent = `${fmtDate(next.ms, { weekday: "short", day: "numeric", month: "short" })} · ${clock(next.ms)} UTC`;
  }
  const up = s.earth_elevation_deg + s.earth_radius_deg > 0;
  $("earth-stat").textContent = `${fmtNumber(s.earth_lit_fraction * 100, 0)}% lit`;
  $("earth-note").textContent = up
    ? `${phase}, ${fmtNumber(s.earth_elevation_deg, 0)}° up in the ${direction(s.earth_azimuth_deg)}`
    : reach === "never"
      ? `${phase}, always below the horizon here`
      : `${phase}, below the horizon`;
  $("place-name").textContent = p.name;
  $("place-kind").textContent = `${p.group === "Your point" ? coord(p.latitude, p.longitude) : p.group.replace(/s$/, "")} · ${setting(p)}`;
  document.title = `${p.name} · Open Moon Skies`;
}

// ---------- Timeline ----------
const REGIME_RGB = REGIMES.map((r) => {
  const h = r.colour.slice(1);
  return [0, 2, 4].map((i) => parseInt(h.slice(i, i + 2), 16));
});
const TIMELINE_HALF = 15 * DAY_MS;
function sampleAt(t) {
  const sp = app.span;
  if (!sp) return null;
  const i = Math.round((t - sp.start) / sp.step);
  if (i < 0 || i >= sp.samples.ms.length) return null;
  return i;
}
function drawTimeline() {
  const c = $("timeline-canvas");
  const r = c.getBoundingClientRect();
  const d = Math.min(window.devicePixelRatio || 1, 2);
  const W = Math.round(r.width * d),
    H = Math.round(r.height * d);
  if (c.width !== W || c.height !== H) {
    c.width = W;
    c.height = H;
  }
  const ctx = c.getContext("2d");
  ctx.clearRect(0, 0, W, H);
  const half = r.width < 600 ? 8 * DAY_MS : TIMELINE_HALF;
  const t0 = app.ms - half,
    t1 = app.ms + half;
  const x = (t) => ((t - t0) / (t1 - t0)) * W;
  const top = 16 * d,
    bh = 20 * d;
  const sp = app.span;
  if (sp) {
    for (let px = 0; px < W; px++) {
      const t = t0 + (px / W) * (t1 - t0);
      const i = sampleAt(t);
      if (i === null) continue;
      const code = sp.samples.regime[i];
      const L = sp.samples.total[i];
      const k = code === 0 ? 0.75 + 0.25 * clamp(Math.log10(L) / 5, 0, 1) : code === 1 ? 0.55 + 0.4 * clamp((Math.log10(L) - 0.47) / 3.3, 0, 1) : code === 4 ? 0.7 : 0.85;
      const [rr, gg, bb] = REGIME_RGB[code];
      ctx.fillStyle = `rgb(${rr * k},${gg * k},${bb * k})`;
      ctx.fillRect(px, top, 1, bh);
    }
    // Event markers.
    for (const e of sp.events) {
      if (e.ms < t0 || e.ms > t1) continue;
      if (!["sunrise", "sunset", "full-earth", "new-earth", "eclipse-begins", "dusk", "dawn"].includes(e.type)) continue;
      const xx = x(e.ms);
      ctx.fillStyle = EVENT_INFO[e.type].colour;
      ctx.beginPath();
      ctx.arc(xx, top - 7 * d, 3.6 * d, 0, Math.PI * 2);
      ctx.fill();
      ctx.fillRect(xx - 0.5 * d, top - 4 * d, 1 * d, 4 * d);
    }
  }
  // Day ticks and labels.
  ctx.fillStyle = "rgba(255,255,255,0.75)";
  ctx.font = `${600} ${11 * d}px system-ui, sans-serif`;
  ctx.textAlign = "center";
  const firstDay = Math.ceil(t0 / DAY_MS) * DAY_MS;
  const every = r.width < 600 ? 2 : r.width < 1000 ? 3 : 2;
  for (let t = firstDay; t <= t1; t += DAY_MS) {
    const xx = x(t);
    ctx.fillStyle = "rgba(255,255,255,0.35)";
    ctx.fillRect(xx, top + bh, 1 * d, 4 * d);
    const day = new Date(t).getUTCDate();
    if ((Math.round(t / DAY_MS) % every) === 0 || day === 1) {
      ctx.fillStyle = day === 1 ? "#fff" : "rgba(255,255,255,0.7)";
      ctx.fillText(day === 1 ? fmtDate(t, { month: "short" }) : String(day), xx, top + bh + 15 * d);
    }
  }
  // Playhead.
  const cx = W / 2;
  ctx.fillStyle = "#fff";
  ctx.fillRect(cx - 1 * d, top - 12 * d, 2 * d, bh + 16 * d);
  ctx.beginPath();
  ctx.moveTo(cx - 6 * d, top - 14 * d);
  ctx.lineTo(cx + 6 * d, top - 14 * d);
  ctx.lineTo(cx, top - 7 * d);
  ctx.fill();
  c.dataset.t0 = t0;
  c.dataset.t1 = t1;
  c.setAttribute("aria-valuetext", `${longDate(app.ms)} ${clock(app.ms)} UTC`);
}
function bindTimeline() {
  const c = $("timeline-canvas");
  let drag = null;
  c.addEventListener("pointerdown", (e) => {
    c.setPointerCapture(e.pointerId);
    drag = { x: e.clientX, ms: app.ms, moved: false };
    stop();
  });
  c.addEventListener("pointermove", (e) => {
    if (!drag) return;
    const r = c.getBoundingClientRect();
    const span = Number(c.dataset.t1) - Number(c.dataset.t0);
    const dx = e.clientX - drag.x;
    if (Math.abs(dx) > 3) drag.moved = true;
    setTime(drag.ms - (dx / r.width) * span);
  });
  c.addEventListener("pointerup", (e) => {
    if (drag && !drag.moved) {
      const r = c.getBoundingClientRect();
      const t = Number(c.dataset.t0) + ((e.clientX - r.left) / r.width) * (Number(c.dataset.t1) - Number(c.dataset.t0));
      setTime(t, { heavy: true });
    } else scheduleUpdate({ heavy: true });
    drag = null;
  });
  c.addEventListener("keydown", (e) => {
    const step = e.shiftKey ? DAY_MS : HOUR;
    if (e.key === "ArrowLeft") setTime(app.ms - step, { heavy: true });
    else if (e.key === "ArrowRight") setTime(app.ms + step, { heavy: true });
    else return;
    e.preventDefault();
  });
}

// ---------- Playback ----------
let lastTick = 0;
function tick(now) {
  if (!app.playing) return;
  const dt = lastTick ? Math.min(now - lastTick, 100) : 16;
  lastTick = now;
  const next = app.ms + (dt / 1000) * app.speed * HOUR;
  if (next >= maxMs) {
    stop();
    return;
  }
  app.ms = next;
  update();
  requestAnimationFrame(tick);
}
function play() {
  app.playing = true;
  lastTick = 0;
  $("play").classList.add("playing");
  $("play").setAttribute("aria-label", "Pause");
  requestAnimationFrame(tick);
}
function stop() {
  if (!app.playing) return;
  app.playing = false;
  $("play").classList.remove("playing");
  $("play").setAttribute("aria-label", "Play");
  scheduleUpdate({ heavy: true });
}
function bindTransport() {
  $("play").addEventListener("click", () => (app.playing ? stop() : play()));
  $("step-back").addEventListener("click", () => setTime(app.ms - DAY_MS, { heavy: true }));
  $("step-forward").addEventListener("click", () => setTime(app.ms + DAY_MS, { heavy: true }));
  $("now-button").addEventListener("click", () => {
    stop();
    setTime(Date.now(), { heavy: true });
  });
  for (const b of document.querySelectorAll(".speeds button"))
    b.addEventListener("click", () => {
      app.speed = Number(b.dataset.speed);
      for (const o of document.querySelectorAll(".speeds button")) o.setAttribute("aria-checked", String(o === b));
    });
  const fromInputs = () => {
    const date = $("date").value,
      time = $("time").value || "12:00";
    if (!date) return;
    const t = Date.parse(`${date}T${time.length === 5 ? time + ":00" : time}Z`);
    if (Number.isFinite(t)) {
      stop();
      setTime(t, { heavy: true });
    }
  };
  $("date").addEventListener("change", fromInputs);
  $("time").addEventListener("change", fromInputs);
  document.addEventListener("keydown", (e) => {
    if (e.target.closest("input, textarea, select, dialog")) return;
    if (e.key === " " && e.target === document.body) {
      e.preventDefault();
      app.playing ? stop() : play();
    }
  });
}
function syncInputs() {
  const iso = new Date(app.ms).toISOString();
  if (document.activeElement !== $("date")) $("date").value = iso.slice(0, 10);
  if (document.activeElement !== $("time")) $("time").value = iso.slice(11, 16);
}

// ---------- Cards ----------
let earthCardTone = new ToneCurve(1);
function updateCards() {
  const s = app.state;
  $("card-lux").textContent = lux(s.total);
  const f = s.total > 0 ? s.solar / s.total : 1;
  $("mix-sun").style.width = `${f * 100}%`;
  $("mix-earth").style.width = `${(1 - f) * 100}%`;
  $("card-sun-lux").textContent = `${lux(s.solar)} lux`;
  $("card-earth-lux").textContent = `${lux(s.earth)} lux`;
  $("card-compare").textContent = `${earthComparison(s.total, assets.sky.earth)}.`;
  const ratio = moonlightRatio(s.total);
  $("card-moonlight").textContent = ratio ? `That is ${ratio}.` : s.total >= 0.4 ? "" : "A full Moon on Earth gives at most about 0.3 lux.";

  // Sun card.
  const up = s.sun_elevation_deg + s.sun_radius_deg > 0;
  $("sun-where").textContent = up
    ? `${fmtNumber(s.sun_elevation_deg, 0)}° up in the ${direction(s.sun_azimuth_deg)}`
    : `${fmtNumber(-s.sun_elevation_deg, 0)}° below the ${direction(s.sun_azimuth_deg)} horizon`;
  const box = $("sun-countdowns");
  box.replaceChildren();
  const pick = up ? ["sunset", "sunrise"] : ["sunrise", "dusk", "dawn"];
  for (const type of pick) {
    const [, after] = eventsAround(type);
    if (!after) continue;
    const row = el("div", { class: "countdown" });
    row.innerHTML = icon(type);
    row.append(
      el("div", {}, el("strong", {}, `${EVENT_INFO[type].name} in ${shortDuration(after.ms - app.ms)}`), el("span", {}, `${fmtDate(after.ms, { weekday: "short", day: "numeric", month: "short" })}, ${clock(after.ms)} UTC`)),
    );
    box.append(row);
    if (box.children.length >= 2) break;
  }
  const ev = app.span?.events ?? [];
  const rise = ev.find((e) => e.type === "sunrise");
  const set = rise && ev.find((e) => e.type === "sunset" && e.ms > rise.ms);
  const rise2 = set && ev.find((e) => e.type === "sunrise" && e.ms > set.ms);
  $("sun-day-length").textContent =
    rise && set && rise2
      ? `Here the Sun stays up for ${duration(set.ms - rise.ms)}, then stays down for ${duration(rise2.ms - set.ms)}.`
      : "";

  // Earth card.
  const reach = earthReach();
  const phase = earthPhase(s.earth_lit_fraction, isWaxing());
  $("earth-phase").textContent = `${phase} · ${fmtNumber(s.earth_lit_fraction * 100, 0)}% lit`;
  const eUp = s.earth_elevation_deg + s.earth_radius_deg > 0;
  $("earth-where").textContent =
    reach === "never"
      ? "Earth never rises over this part of the Moon; the picture shows the face it turns toward the far side."
      : `${eUp ? `${fmtNumber(s.earth_elevation_deg, 0)}° up in the ${direction(s.earth_azimuth_deg)}` : `${fmtNumber(-s.earth_elevation_deg, 0)}° below the horizon`}. ${reach === "always" ? "From here Earth never sets: it hangs in nearly the same spot all month." : "From here Earth rises and sets as the Moon nods on its orbit."}`;
  const ratio2 = (s.earth_radius_deg * 2) / MOON_FROM_EARTH_DEG;
  const [lonE, latE] = subObserverPoint(app.frame, s.earth_elevation_deg, s.earth_azimuth_deg);
  $("earth-size").textContent = `${fmtNumber(s.earth_radius_deg * 2, 2)}° across, ${fmtNumber(ratio2, 1)} times the width of the Moon seen from Earth. Facing you: ${coord(latE, lonE)}.`;
  drawEarthCard();
  renderUpcoming();
}
function drawEarthCard() {
  const c = $("earth-big");
  paintEarth(c, { exposure: 4, tone: earthCardTone });
  c.setAttribute("aria-label", `Earth, ${fmtNumber(app.state.earth_lit_fraction * 100, 0)} percent lit`);
}
function renderUpcoming() {
  const list = $("upcoming");
  list.replaceChildren();
  const types = ["sunrise", "sunset", "dusk", "dawn", "full-earth", "new-earth", "eclipse-begins", "earthrise", "earthset"];
  const ev = (app.span?.events ?? []).filter((e) => e.ms > app.ms && types.includes(e.type)).slice(0, 6);
  if (!ev.length) list.append(el("li", {}, "Calculating…"));
  for (const e of ev) {
    const li = el("li");
    const b = el("button", { type: "button", "aria-label": `Go to ${EVENT_INFO[e.type].name}, ${longDate(e.ms)}` });
    const ic = el("span");
    ic.innerHTML = icon(e.type);
    b.append(
      ic,
      el("span", { class: "ev-name" }, EVENT_INFO[e.type].name, el("span", { class: "ev-date" }, `${fmtDate(e.ms, { weekday: "short", day: "numeric", month: "short" })} · ${clock(e.ms)}`)),
      el("span", { class: "ev-in" }, `in ${shortDuration(e.ms - app.ms)}`),
    );
    b.addEventListener("click", () => {
      stop();
      setTime(e.ms + (e.type === "sunrise" ? 3 * HOUR : e.type === "sunset" ? -3 * HOUR : 0), { heavy: true });
      window.scrollTo({ top: 0, behavior: "smooth" });
    });
    li.append(b);
    list.append(li);
  }
}

// ---------- Dial ----------
function updateDial() {
  const sp = app.span;
  const svg = $("dial");
  if (!sp) return;
  const ev = sp.events;
  let start = ev.filter((e) => e.type === "sunrise" && e.ms <= app.ms).at(-1)?.ms;
  let end = start && ev.find((e) => e.type === "sunrise" && e.ms > start)?.ms;
  let anchored = true;
  if (!start || !end) {
    start = app.ms - 14.765 * DAY_MS;
    end = app.ms + 14.765 * DAY_MS;
    anchored = false;
  }
  if (app.dialKey === `${start}-${end}-${app.place.id}`) {
    moveDialHand(start, end);
    return;
  }
  app.dialKey = `${start}-${end}-${app.place.id}`;
  svg.replaceChildren();
  const cx = 200,
    cy = 200,
    R = 150,
    w = 34;
  const ang = (t) => ((t - start) / (end - start)) * 2 * Math.PI - Math.PI / 2;
  const pt = (a, r) => [cx + r * Math.cos(a), cy + r * Math.sin(a)];
  const defs = svgEl("defs");
  const glow = svgEl("filter", { id: "dial-glow", x: "-50%", y: "-50%", width: "200%", height: "200%" });
  glow.append(svgEl("feGaussianBlur", { stdDeviation: 6 }));
  defs.append(glow);
  svg.append(defs);
  svg.append(svgEl("circle", { cx, cy, r: R + w / 2 + 14, fill: "rgba(255,255,255,0.03)", stroke: "rgba(255,255,255,0.08)" }));
  // Regime arcs.
  const i0 = Math.max(0, Math.floor((start - sp.start) / sp.step)),
    i1 = Math.min(sp.samples.ms.length - 1, Math.ceil((end - sp.start) / sp.step));
  let runStart = i0;
  const arcs = [];
  for (let i = i0 + 1; i <= i1 + 1; i++) {
    if (i > i1 || sp.samples.regime[i] !== sp.samples.regime[runStart]) {
      arcs.push([sp.samples.ms[runStart], sp.samples.ms[Math.min(i, i1)], sp.samples.regime[runStart]]);
      runStart = i;
    }
  }
  const g = svgEl("g");
  for (const [a, b, code] of arcs) {
    const t0 = Math.max(a, start),
      t1 = Math.min(b, end);
    if (t1 <= t0) continue;
    const a0 = ang(t0),
      a1 = ang(t1);
    const large = a1 - a0 > Math.PI ? 1 : 0;
    const [x0, y0] = pt(a0, R),
      [x1, y1] = pt(a1, R);
    const path = svgEl("path", {
      d: `M${x0},${y0} A${R},${R} 0 ${large} 1 ${x1},${y1}`,
      fill: "none",
      stroke: REGIMES[code].colour,
      "stroke-width": w,
      class: "dial-seg",
    });
    path.append(svgEl("title", {}, `${REGIMES[code].name}: ${shortDate(t0)} to ${shortDate(t1)}`));
    path.addEventListener("click", (e) => {
      const r = svg.getBoundingClientRect();
      const x = ((e.clientX - r.left) / r.width) * 400 - cx,
        y = ((e.clientY - r.top) / r.height) * 400 - cy;
      let a = Math.atan2(y, x) + Math.PI / 2;
      if (a < 0) a += 2 * Math.PI;
      stop();
      setTime(start + (a / (2 * Math.PI)) * (end - start), { heavy: true });
    });
    g.append(path);
  }
  svg.append(g);
  // Day ticks.
  const days = (end - start) / DAY_MS;
  for (let k = 0; k <= Math.floor(days); k++) {
    const a = ang(start + k * DAY_MS);
    const [x0, y0] = pt(a, R - w / 2 - 4),
      [x1, y1] = pt(a, R - w / 2 - (k % 5 === 0 ? 14 : 8));
    svg.append(svgEl("line", { x1: x0, y1: y0, x2: x1, y2: y1, stroke: "rgba(255,255,255,0.45)", "stroke-width": k % 5 === 0 ? 2 : 1 }));
    if (k % 5 === 0 && k > 0) {
      const [tx, ty] = pt(a, R - w / 2 - 28);
      svg.append(svgEl("text", { x: tx, y: ty + 4, "text-anchor": "middle", fill: "rgba(255,255,255,0.65)", "font-size": 12, "font-weight": 600 }, `day ${k}`));
    }
  }
  // Event markers on the outside.
  const placed = [];
  for (const e of ev) {
    if (e.ms < start || e.ms > end) continue;
    if (!["sunset", "dusk", "dawn", "full-earth", "new-earth", "eclipse-begins", "sunrise"].includes(e.type)) continue;
    const a = ang(e.ms);
    if (placed.some((p) => Math.abs(p - a) < 0.17)) continue;
    placed.push(a);
    const [x, y] = pt(a, R + w / 2 + 16);
    const fo = svgEl("g", { transform: `translate(${x - 13},${y - 13})` });
    fo.innerHTML = icon(e.type).replace("<svg ", '<svg width="26" height="26" ');
    fo.append(svgEl("title", {}, `${EVENT_INFO[e.type].name}: ${longDate(e.ms)} ${clock(e.ms)} UTC`));
    svg.append(fo);
  }
  // Hand and centre.
  svg.append(svgEl("g", { id: "dial-hand" }));
  const centre = svgEl("g", { id: "dial-centre" });
  svg.append(centre);
  app.dialAnchored = anchored;
  moveDialHand(start, end);
  renderCycleFacts(start, end, arcs, anchored);
}
function moveDialHand(start, end) {
  const hand = $("dial-hand"),
    centre = $("dial-centre");
  if (!hand) return;
  const a = ((app.ms - start) / (end - start)) * 2 * Math.PI - Math.PI / 2;
  const cx = 200,
    cy = 200;
  hand.replaceChildren(
    svgEl("line", { x1: cx + 70 * Math.cos(a), y1: cy + 70 * Math.sin(a), x2: cx + 186 * Math.cos(a), y2: cy + 186 * Math.sin(a), stroke: "#fff", "stroke-width": 3, "stroke-linecap": "round" }),
    svgEl("circle", { cx: cx + 150 * Math.cos(a), cy: cy + 150 * Math.sin(a), r: 9, fill: "#fff", stroke: "#0b1026", "stroke-width": 3 }),
  );
  const day = (app.ms - start) / DAY_MS;
  const code = regimeOf(app.state);
  centre.replaceChildren(
    svgEl("text", { x: cx, y: cy - 22, "text-anchor": "middle", fill: "rgba(255,255,255,0.7)", "font-size": 13, "font-weight": 700, "letter-spacing": "0.14em" }, app.dialAnchored ? "SINCE SUNRISE" : "THIS CYCLE"),
    svgEl("text", { x: cx, y: cy + 22, "text-anchor": "middle", fill: "#fff", "font-size": 50, "font-weight": 800 }, app.dialAnchored ? `Day ${Math.floor(day) + 1}` : `${fmtNumber(day, 1)} d`),
    svgEl("text", { x: cx, y: cy + 48, "text-anchor": "middle", fill: REGIMES[code].colour, "font-size": 15, "font-weight": 700 }, `${REGIMES[code].name} · of ${fmtNumber((end - start) / DAY_MS, 1)}`),
  );
}
function renderCycleFacts(start, end, arcs, anchored) {
  const total = (code) => arcs.filter((a) => a[2] === code).reduce((s, a) => s + Math.max(0, Math.min(a[1], end) - Math.max(a[0], start)), 0);
  const box = $("cycle-facts");
  box.replaceChildren();
  const facts = [
    [0, "Daylight", total(0), "The Sun crosses the sky in about two weeks."],
    [1, "Golden twilight", total(1), "Bright enough to work outside without lamps after the Sun has set or before it rises."],
    [2, "Earthlit night", total(2), "Earth outshines the last of the twilight."],
    [3, "Night", total(3) + total(4), "Light below the end of Earth's civil twilight."],
  ];
  for (const [code, name, ms, text] of facts) {
    if (ms <= 0 && code !== 0) continue;
    box.append(
      el("div", { class: "fact", style: `--c:${REGIMES[code].colour}` }, el("b", {}, ms > 0 ? duration(ms) : "none"), el("span", {}, name), el("p", {}, text)),
    );
  }
  $("cycle-title").textContent = anchored ? `A day here lasts ${fmtNumber((end - start) / DAY_MS, 1)} Earth days` : "Light through one lunar cycle";
  $("cycle-lede").textContent = `${app.place.name}, from one sunrise to the next. The Moon turns once for every trip around Earth, so each day and night lasts about two weeks. Choose any part of the ring to jump there.`;
}

// ---------- Month ----------
function renderMonth() {
  const m = app.month;
  if (!m) return;
  $("month-label").textContent = fmtDate(m.start, { month: "long", year: "numeric" });
  $("month-title").textContent = `${fmtDate(m.start, { month: "long" })} at ${app.place.name}`;
  const grid = $("month-grid");
  grid.replaceChildren();
  for (const w of ["Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun"]) grid.append(el("div", { class: "weekday", role: "columnheader" }, w));
  const first = (new Date(m.start).getUTCDay() + 6) % 7;
  for (let i = 0; i < first; i++) grid.append(el("div", { class: "day empty", "aria-hidden": "true" }));
  const today = new Date().toISOString().slice(0, 10),
    selected = new Date(app.ms).toISOString().slice(0, 10);
  for (const day of m.daily) {
    const iso = new Date(day.ms).toISOString().slice(0, 10);
    const stops = day.regimes.map((code, i) => {
      const L = day.strip[Math.min(day.strip.length - 1, i)];
      const k = code === 0 ? 0.8 + 0.2 * clamp(Math.log10(L) / 5, 0, 1) : code === 1 ? 0.6 + 0.35 * clamp((Math.log10(L) - 0.47) / 3.3, 0, 1) : 0.9;
      const [r, g, b] = REGIME_RGB[code];
      return `rgb(${(r * k) | 0},${(g * k) | 0},${(b * k) | 0}) ${((i + 0.5) / day.regimes.length) * 100}%`;
    });
    const cls = ["day"];
    if (iso === today) cls.push("today");
    if (iso === selected) cls.push("selected");
    const tile = el("button", { type: "button", class: cls.join(" "), "aria-label": `${longDate(day.ms)}: ${dayDescription(day)}` });
    tile.append(el("span", { class: "band-fill", style: `background:linear-gradient(90deg,${stops.join(",")})` }));
    tile.append(el("span", { class: "num" }, String(new Date(day.ms).getUTCDate())));
    tile.append(el("span", { class: "mini" }, `${lux(day.max)} lux max`));
    const icons = el("span", { class: "icons" });
    const seen = new Set();
    for (const e of day.events) {
      if (!["sunrise", "sunset", "dusk", "dawn", "full-earth", "new-earth", "eclipse-begins"].includes(e.type) || seen.has(e.type)) continue;
      seen.add(e.type);
      const i = el("span", { title: `${EVENT_INFO[e.type].name} ${clock(e.ms)} UTC` });
      i.innerHTML = icon(e.type);
      icons.append(i);
    }
    tile.append(icons);
    tile.addEventListener("click", () => {
      stop();
      setTime(day.ms + 12 * HOUR, { heavy: true });
      window.scrollTo({ top: 0, behavior: "smooth" });
    });
    grid.append(tile);
  }
}
function dayDescription(day) {
  const names = [...new Set(day.regimes)].map((c) => REGIMES[c].short.toLowerCase());
  const ev = day.events.map((e) => EVENT_INFO[e.type]?.name).filter(Boolean);
  return `${names.join(", ")}; ${lux(day.min)} to ${lux(day.max)} lux${ev.length ? "; " + ev.join(", ") : ""}`;
}
function legend() {
  const box = $("regime-legend");
  for (const r of REGIMES) box.append(el("span", {}, el("i", { style: `background:${r.colour}` }), r.name));
  for (const t of ["sunrise", "sunset", "full-earth", "new-earth"]) {
    const s = el("span");
    s.innerHTML = icon(t).replace("<svg ", '<svg width="20" height="20" ');
    s.append(EVENT_INFO[t].name);
    box.append(s);
  }
}

// ---------- Places ----------
let globe;
const TABS = ["Study coasts", "Landing sites and heritage", "Seas", "Island summits"];
let activeTab = "Landing sites and heritage";
function renderTabs() {
  const box = $("place-tabs");
  box.replaceChildren();
  for (const t of TABS) {
    const b = el("button", { type: "button", role: "tab", "aria-selected": String(t === activeTab) }, t === "Landing sites and heritage" ? "Landing sites" : t);
    b.addEventListener("click", () => {
      activeTab = t;
      renderTabs();
      renderPlaceList();
    });
    box.append(b);
  }
}
function placeRow(p, onPick) {
  const li = el("li");
  const b = el("button", { type: "button", "aria-current": String(app.place && app.place.id === p.id) });
  b.append(
    el("span", { class: "pdot", style: `background:${GROUP_COLOURS[p.group]}` }),
    el("span", { class: "pname" }, p.name),
    el("span", { class: "pmeta" }, p.year ? String(p.year) : p.status === "submerged" ? "under the sea" : coord(p.latitude, p.longitude)),
  );
  b.addEventListener("click", () => onPick(p));
  li.append(b);
  return li;
}
function renderPlaceList() {
  const list = $("place-list");
  list.replaceChildren();
  for (const p of allPlaces().filter((p) => p.group === activeTab)) list.append(placeRow(p, (q) => choosePlace(q, { spin: true })));
}
function renderPlaceCard() {
  const p = app.place,
    card = $("place-card");
  card.replaceChildren();
  card.append(el("h3", {}, p.name));
  const tag = [p.group === "Your point" ? "Your point" : p.group.replace(/s$/, ""), coord(p.latitude, p.longitude)];
  if (p.kind) tag.unshift(`${p.kind[0].toUpperCase()}${p.kind.slice(1)}${p.year ? `, ${p.year}` : ""}`);
  card.append(el("p", { class: "tagline" }, tag.join(" · ")));
  if (p.marks) card.append(el("p", { class: "marks" }, p.marks));
  else if (p.description) card.append(el("p", { class: "marks" }, p.description));
  if (p.group === "Seas") card.append(el("p", { class: "marks" }, `${fmtNumber(p.area_km2 / 1000, 0)} thousand km², on average ${fmtNumber(p.mean_depth_m)} m deep and ${fmtNumber(p.max_depth_m)} m at its deepest.${p.also?.length ? ` Takes in ${p.also.slice(0, 4).join(", ")}.` : ""}`));
  if (p.group === "Island summits") card.append(el("p", { class: "marks" }, `An island of ${fmtNumber(p.area_km2)} km², rising ${fmtNumber(p.summit_m)} m above the sea${p.also?.length ? `, with ${p.also.join(", ")}` : ""}.`));
  const chips = el("div", { class: "chips" });
  chips.append(el("span", { class: "chip" }, setting(p)));
  const reach = earthReach();
  if (reach) chips.append(el("span", { class: "chip" }, el("i", { class: "dot earth" }), reach === "always" ? "Earth always up" : reach === "never" ? "Earth never rises" : "Earth rises and sets"));
  chips.append(el("span", { class: "chip" }, el("i", { class: "dot", style: `background:${REGIMES[regimeOf(app.state)].colour}` }), REGIMES[regimeOf(app.state)].name, " now"));
  card.append(chips);
}
function choosePlace(p, opts = {}) {
  if (p.custom && !p.name) p = customPlace(p.longitude, p.latitude);
  app.place = p;
  app.sectors = sectorsFor(p.longitude, p.latitude);
  app.span = null;
  app.month = null;
  app.monthKey = null;
  app.viewYear = null;
  app.yearLocked = false;
  app.dialKey = null;
  spanToken++;
  globe.selected = p;
  if (opts.spin) globe.spinTo(p.longitude, p.latitude, (done) => renderGlobe(done));
  renderPlaceList();
  if (!opts.keepView) {
    const s = stateAt(app.ms, p.longitude, p.latitude, astronomy, transfer, true, "UTC");
    app.state = s;
    aimCamera();
  }
  scheduleUpdate({ heavy: true });
  requestSpan().then(() => renderPlaceCard());
  renderPlaceCard();
}
let globeQueued = false;
function renderGlobe() {
  if (globeQueued) return;
  globeQueued = true;
  requestAnimationFrame(() => {
    globeQueued = false;
    globe.render(transfer);
  });
}
function updateGlobeLight() {
  const g = ephemeris(julianDate(app.ms, settings.scale), astronomy);
  // Earth's light at the point under it, for shading the near side.
  const sub = [Math.atan2(g.earth[1], g.earth[0]) / rad, Math.asin(g.earth[2]) / rad];
  const s = stateAt(app.ms, sub[0], sub[1], astronomy, transfer, true, settings.scale, 4);
  globe.setLight(g.sun, g.earth, s.solar, s.earth);
  renderGlobe();
}
function bindPlaces() {
  globe = new Globe($("globe"), moonTexture, {
    onPick: (p) => {
      if (p.custom) p = customPlace(p.longitude, p.latitude);
      choosePlace(p);
    },
    onHover: (p, e) => {
      const tip = $("globe-tip");
      if (!p) {
        tip.hidden = true;
        return;
      }
      const r = $("globe").getBoundingClientRect(),
        n = globe.canvas.width;
      const q = globe.project(p.longitude, p.latitude, n);
      tip.textContent = p.name;
      tip.style.left = `${(q.x / n) * r.width}px`;
      tip.style.top = `${(q.y / n) * r.height}px`;
      tip.hidden = false;
    },
  });
  globe.onDrag = renderGlobe;
  globe.drawOverlayOnly = renderGlobe;
  globe.places = allPlaces().map((p) => Object.assign(p, { colour: GROUP_COLOURS[p.group] }));
  new ResizeObserver(renderGlobe).observe($("globe"));
  renderTabs();
  // Dialog.
  const dialog = $("place-dialog");
  const fill = () => {
    const q = $("place-search").value.trim().toLowerCase();
    const list = $("dialog-list");
    list.replaceChildren();
    for (const group of TABS) {
      const items = allPlaces().filter((p) => p.group === group && (!q || p.name.toLowerCase().includes(q) || (p.also ?? []).some((a) => a.toLowerCase().includes(q)) || (p.water_body ?? "").toLowerCase().includes(q)));
      if (!items.length) continue;
      list.append(el("li", { class: "group" }, group));
      for (const p of items)
        list.append(
          placeRow(p, (q2) => {
            dialog.close();
            choosePlace(q2, { spin: true });
          }),
        );
    }
  };
  $("place-button").addEventListener("click", () => {
    $("place-search").value = "";
    fill();
    dialog.showModal();
    $("place-search").focus();
  });
  $("place-search").addEventListener("input", fill);
}

// ---------- Twilight comparison ----------
function renderCompare() {
  const angle = Number($("sun-angle").value);
  $("sun-angle-value").textContent =
    angle > 0 ? `${fmtNumber(angle, 1)}° above the horizon` : angle === 0 ? "on the horizon" : `${fmtNumber(-angle, 1)}° below the horizon`;
  const luxAt = (h, a) => {
    const s = h.suns,
      L = h.horizontal_lux;
    for (let i = 0; i < s.length - 1; i++)
      if (a >= s[i] && a <= s[i + 1]) {
        const t = (a - s[i]) / (s[i + 1] - s[i]);
        return Math.exp(Math.log(Math.max(L[i], 1e-30)) * (1 - t) + Math.log(Math.max(L[i + 1], 1e-30)) * t);
      }
    return 0;
  };
  const moonLux = luxAt(assets.sky.moon, angle),
    earthLux = luxAt(assets.sky.earth, angle);
  const tone = new ToneCurve(Math.max(moonLux, earthLux));
  renderDome($("dome-moon"), moonAtlas, angle, tone);
  renderDome($("dome-earth"), earthAtlas, angle, tone);
  const when = (degPerHour) => {
    const h = Math.abs(angle) / degPerHour;
    if (angle === 0) return "at sunset";
    const text = h >= 1 ? duration(h * HOUR) : `${Math.round(h * 60)} minutes`;
    return angle < 0 ? `${text} after sunset` : `${text} before sunset`;
  };
  $("dome-moon-text").textContent = `${lux(moonLux)} lux · ${when(SYNODIC_DEG_PER_HOUR)}`;
  $("dome-earth-text").textContent = `${lux(earthLux)} lux · ${when(15)}`;
  drawFadeChart(angle);
}
function drawFadeChart(angle) {
  const svg = $("fade-chart");
  svg.replaceChildren();
  const W = 520,
    H = 300,
    L = 56,
    Rm = 14,
    T = 16,
    B = 40;
  const hours = 100;
  const x = (h) => L + (h / hours) * (W - L - Rm);
  const y = (v) => T + ((5 - Math.log10(Math.max(v, 1e-4))) / 9) * (H - T - B);
  for (let e = -4; e <= 5; e++) {
    svg.append(svgEl("line", { x1: L, x2: W - Rm, y1: y(10 ** e), y2: y(10 ** e), stroke: "rgba(255,255,255,0.08)" }));
    svg.append(svgEl("text", { x: L - 8, y: y(10 ** e) + 4, "text-anchor": "end", fill: "rgba(255,255,255,0.55)", "font-size": 11 }, e >= 3 ? `${10 ** (e - 3)}k` : e < 0 ? `${10 ** e}` : `${10 ** e}`));
  }
  for (let h = 0; h <= hours; h += 12) {
    svg.append(svgEl("text", { x: x(h), y: H - B + 18, "text-anchor": "middle", fill: "rgba(255,255,255,0.55)", "font-size": 11 }, `${h} h`));
  }
  svg.append(svgEl("text", { x: (L + W) / 2, y: H - 4, "text-anchor": "middle", fill: "rgba(255,255,255,0.7)", "font-size": 12, "font-weight": 600 }, "Hours after sunset"));
  // Threshold.
  svg.append(svgEl("line", { x1: L, x2: W - Rm, y1: y(PRACTICAL_DUSK_LUX), y2: y(PRACTICAL_DUSK_LUX), stroke: "#b07cff", "stroke-dasharray": "5 5", "stroke-width": 1.5 }));
  svg.append(svgEl("text", { x: W - Rm - 4, y: y(PRACTICAL_DUSK_LUX) - 6, "text-anchor": "end", fill: "#c9a8ff", "font-size": 11, "font-weight": 600 }, "End of civil twilight on Earth, 2.98 lux"));
  const moonCurve = [],
    earthCurve = [];
  const moonLuxAt = (deg) =>
    interpolate(transfer.direct_elevation_deg, transfer.solar.direct_lux, deg) + interpolate(transfer.diffuse_elevation_deg, transfer.solar.diffuse_lux, deg);
  const eh = assets.sky.earth;
  const earthLuxAt = (deg) => {
    for (let i = 0; i < eh.suns.length - 1; i++)
      if (deg >= eh.suns[i] && deg <= eh.suns[i + 1]) {
        const t = (deg - eh.suns[i]) / (eh.suns[i + 1] - eh.suns[i]);
        return Math.exp(Math.log(Math.max(eh.horizontal_lux[i], 1e-30)) * (1 - t) + Math.log(Math.max(eh.horizontal_lux[i + 1], 1e-30)) * t);
      }
    return 0;
  };
  let moonCross = null,
    earthCross = null;
  for (let h = 0; h <= hours; h += 0.25) {
    const m = moonLuxAt(-h * SYNODIC_DEG_PER_HOUR);
    moonCurve.push(`${moonCurve.length ? "L" : "M"}${x(h).toFixed(1)},${y(m).toFixed(1)}`);
    if (moonCross === null && m < PRACTICAL_DUSK_LUX) moonCross = h;
  }
  for (let h = 0; h <= 2; h += 0.01) {
    const e = earthLuxAt(-h * 15);
    if (e < 1e-4) break;
    earthCurve.push(`${earthCurve.length ? "L" : "M"}${x(h).toFixed(1)},${y(e).toFixed(1)}`);
    if (earthCross === null && e < PRACTICAL_DUSK_LUX) earthCross = h;
  }
  svg.append(svgEl("path", { d: moonCurve.join(" "), fill: "none", stroke: "url(#moon-grad)", "stroke-width": 3.5, "stroke-linecap": "round" }));
  const defs = svgEl("defs");
  const grad = svgEl("linearGradient", { id: "moon-grad", x1: 0, x2: 1, y1: 0, y2: 0 });
  grad.append(svgEl("stop", { offset: "0%", "stop-color": "#ffc53d" }), svgEl("stop", { offset: "60%", "stop-color": "#ff8a5c" }), svgEl("stop", { offset: "100%", "stop-color": "#ff5d8f" }));
  defs.append(grad);
  svg.prepend(defs);
  svg.append(svgEl("path", { d: earthCurve.join(" "), fill: "none", stroke: "#38c6f4", "stroke-width": 3.5, "stroke-linecap": "round" }));
  if (moonCross !== null) {
    svg.append(svgEl("circle", { cx: x(moonCross), cy: y(PRACTICAL_DUSK_LUX), r: 5, fill: "#ff8a5c" }));
    svg.append(svgEl("text", { x: x(moonCross) - 8, y: y(PRACTICAL_DUSK_LUX) + 20, "text-anchor": "end", fill: "#ffb08c", "font-size": 12.5, "font-weight": 700 }, `Open Moon: ${Math.round(moonCross)} hours`));
  }
  if (earthCross !== null) {
    svg.append(svgEl("circle", { cx: x(earthCross), cy: y(PRACTICAL_DUSK_LUX), r: 5, fill: "#38c6f4" }));
    svg.append(svgEl("text", { x: x(earthCross) + 10, y: y(PRACTICAL_DUSK_LUX) + 36, fill: "#7fdcff", "font-size": 12.5, "font-weight": 700 }, `Earth: ${Math.round(earthCross * 60)} minutes`));
  }
  // The slider's moment on both curves.
  if (angle <= 0) {
    const hm = -angle / SYNODIC_DEG_PER_HOUR;
    if (hm <= hours) svg.append(svgEl("line", { x1: x(hm), x2: x(hm), y1: T, y2: H - B, stroke: "rgba(255,255,255,0.5)", "stroke-dasharray": "2 4" }));
  }
}

// ---------- Year ----------
function renderYear() {
  const yr = app.year;
  if (!yr) return;
  const c = $("year-canvas");
  c.width = yr.days;
  c.height = 24;
  const ctx = c.getContext("2d");
  const img = ctx.createImageData(yr.days, 24);
  for (let d = 0; d < yr.days; d++)
    for (let h = 0; h < 24; h++) {
      const i = d * 24 + h,
        code = yr.regime[i],
        L = yr.total[i];
      const k =
        code === 0 ? 0.62 + 0.38 * clamp(Math.log10(L) / 5, 0, 1) : code === 1 ? 0.45 + 0.5 * clamp((Math.log10(L) - 0.47) / 3.3, 0, 1) : code === 2 ? 0.55 + 0.4 * clamp(L / 4, 0, 1) : code === 3 ? 0.65 : 0.55;
      const [r, g, b] = REGIME_RGB[code];
      const o = (h * yr.days + d) * 4;
      img.data[o] = r * k;
      img.data[o + 1] = g * k;
      img.data[o + 2] = b * k;
      img.data[o + 3] = 255;
    }
  ctx.putImageData(img, 0, 0);
  const months = $("year-months");
  months.replaceChildren();
  for (let m = 0; m < 12; m++) {
    const t = Date.UTC(app.viewYear, m, 15);
    months.append(el("span", { style: `left:${((t - yr.start) / DAY_MS / yr.days) * 100}%` }, fmtDate(t, { month: "short" })));
  }
  $("year-title").textContent = `A year of light at ${app.place.name}`;
}
function bindYear() {
  const c = $("year-canvas"),
    tip = $("year-tip");
  const at = (e) => {
    const r = c.getBoundingClientRect();
    const d = Math.floor(((e.clientX - r.left) / r.width) * app.year.days),
      h = Math.floor(((e.clientY - r.top) / r.height) * 24);
    return [clamp(d, 0, app.year.days - 1), clamp(h, 0, 23), r];
  };
  c.addEventListener("pointermove", (e) => {
    if (!app.year) return;
    const [d, h, r] = at(e);
    const i = d * 24 + h;
    const t = app.year.start + d * DAY_MS + h * HOUR;
    tip.textContent = `${fmtDate(t, { day: "numeric", month: "short" })} ${String(h).padStart(2, "0")}:00 · ${REGIMES[app.year.regime[i]].short} · ${lux(app.year.total[i])} lux`;
    tip.style.left = `${clamp(e.clientX - r.left, 90, r.width - 90)}px`;
    tip.style.top = `${e.clientY - r.top}px`;
    tip.hidden = false;
  });
  c.addEventListener("pointerleave", () => (tip.hidden = true));
  c.addEventListener("click", (e) => {
    if (!app.year) return;
    const [d, h] = at(e);
    stop();
    app.yearLocked = true;
    setTime(app.year.start + d * DAY_MS + (h + 0.5) * HOUR, { heavy: true });
    window.scrollTo({ top: 0, behavior: "smooth" });
  });
  $("year-prev").addEventListener("click", () => {
    app.viewYear = Math.max(2000, app.viewYear - 1);
    app.yearLocked = true;
    requestYear();
  });
  $("year-next").addEventListener("click", () => {
    app.viewYear = Math.min(2500, app.viewYear + 1);
    app.yearLocked = true;
    requestYear();
  });
}

// ---------- Sharing and calendar export ----------
function queryUpdate() {
  const p = app.place;
  const q = new URLSearchParams();
  if (p.custom) {
    q.set("lat", p.latitude);
    q.set("lon", p.longitude);
  } else q.set("place", p.id);
  q.set("at", new Date(Math.round(app.ms / 60000) * 60000).toISOString().replace(":00.000Z", "Z"));
  history.replaceState(null, "", `${location.pathname}?${q}${location.hash}`);
  $("almanac-link").href = `almanac.html?lat=${p.latitude}&lon=${p.longitude}&at=${new Date(app.ms).toISOString()}`;
  $("almanac-link-2").href = $("almanac-link").href;
}
function toast(text) {
  const t = $("toast");
  t.textContent = text;
  t.hidden = false;
  clearTimeout(toast.timer);
  toast.timer = setTimeout(() => (t.hidden = true), 2600);
}
async function share() {
  const url = location.href;
  try {
    if (navigator.share && matchMedia("(pointer: coarse)").matches) {
      await navigator.share({ title: document.title, url });
      return;
    }
    await navigator.clipboard.writeText(url);
    toast("Link copied");
  } catch {
    toast("Copy the address bar to share this view");
  }
}
async function exportCalendar() {
  const p = app.place;
  const b = $("ics");
  b.disabled = true;
  const label = b.lastChild.textContent;
  b.lastChild.textContent = " Preparing your calendar…";
  try {
    const start = Date.now() > minMs && Date.now() < maxMs - 366 * DAY_MS ? Math.max(app.ms, Date.now()) : app.ms;
    const result = await ask("span", {
      settings: { ...settings, longitude: p.longitude, latitude: p.latitude },
      start,
      end: Math.min(maxMs, start + 365 * DAY_MS),
      step: 2 * HOUR,
      order: 6,
    });
    const stamp = (ms) => new Date(ms).toISOString().replace(/[-:]/g, "").replace(/\.\d+Z$/, "Z");
    const lines = ["BEGIN:VCALENDAR", "VERSION:2.0", "PRODID:-//Terluna//Open Moon Skies//EN", "CALSCALE:GREGORIAN", `X-WR-CALNAME:Open Moon · ${p.name}`];
    for (const e of result.events) {
      if (!["sunrise", "sunset", "dusk", "dawn", "full-earth", "new-earth", "eclipse-begins"].includes(e.type)) continue;
      lines.push(
        "BEGIN:VEVENT",
        `UID:${e.type}-${Math.round(e.ms / 1000)}-${p.latitude}-${p.longitude}@terluna`,
        `DTSTAMP:${stamp(Date.now())}`,
        `DTSTART:${stamp(e.ms)}`,
        `DTEND:${stamp(e.ms + 15 * 60000)}`,
        `SUMMARY:${EVENT_INFO[e.type].name} · ${p.name} (Open Moon)`,
        `DESCRIPTION:${EVENT_INFO[e.type].name} at ${p.name}\\, ${coord(p.latitude, p.longitude)}\\, on the Open Moon. Clear-sky calendar from Terluna's lighting model.`,
        "END:VEVENT",
      );
    }
    lines.push("END:VCALENDAR");
    const blob = new Blob([lines.join("\r\n")], { type: "text/calendar" });
    const a = el("a", { href: URL.createObjectURL(blob), download: `open-moon-${p.name.toLowerCase().replace(/[^a-z0-9]+/g, "-")}.ics` });
    a.click();
    setTimeout(() => URL.revokeObjectURL(a.href), 2000);
    toast("Calendar file ready");
  } catch (e) {
    toast(`Calendar export failed: ${e.message}`);
  } finally {
    b.disabled = false;
    b.lastChild.textContent = label;
  }
}

function credits() {
  const list = $("credits-list");
  const items = [
    ["Light, twilight and Earthlight: Terluna illumination calendar products (solved 1.2-atm column, titania shield, spectral transport)", "almanac.html"],
    ["Sky colour and brightness: Terluna solved spherical sky atlas, Open Moon and Earth control", null],
    ["Earth's spectrum and phase brightness: Glenar et al. 2019; Robinson et al. 2025", "https://doi.org/10.1016/j.icarus.2018.12.025"],
    ...assets.earth_credits.map((c) => [`Earth imagery: ${c.credit}`, c.source]),
    ["Moon terrain: LOLA topography over the GRAIL geoid, through Terluna's geography atlas at 28% water", null],
    ["Feature names: IAU Gazetteer of Planetary Nomenclature (USGS); heritage sites from Terluna's conservation register", null],
    [`Stars: ${stars.source.catalogue}`, stars.source.url],
    ["Ephemeris checks: NASA/JPL Horizons", "https://ssd.jpl.nasa.gov/horizons/"],
    ["Full-Moon illuminance on Earth: Kyba, Mohar & Posch 2017, How bright is moonlight?", "https://doi.org/10.1093/astrogeo/atx025"],
  ];
  for (const [text, href] of items) {
    const li = el("li");
    if (href) li.append(el("a", { href }, text));
    else li.append(text);
    list.append(li);
  }
  if (build) $("build-line").textContent = `Build ${build.schema} · ${Object.keys(build.files ?? {}).length} source files hashed`;
}

function showError(text) {
  $("error").hidden = false;
  $("error").textContent = text;
}

// ---------- Start ----------
async function start() {
  try {
    const query = new URLSearchParams(location.search);
    const at = Date.parse(query.get("at") ?? "");
    if (Number.isFinite(at)) app.ms = clamp(at, minMs, maxMs);
    let header;
    [astronomy, transfer, header, places, stars, build] = await Promise.all([
      json("data/astronomy.json"),
      json("data/transfer.json"),
      json("assets/assets.json"),
      json("assets/places.json"),
      json("assets/stars.json"),
      json("build.json").catch(() => null),
    ]);
    assets = header;
    if (astronomy.schema !== "terluna.illumination.calendar-astronomy/1" || transfer.schema !== "terluna.illumination.calendar-transfer/1")
      throw Error("These model files have an unsupported version.");
    const [moonBin, earthBin, et, mt, hm] = await Promise.all([
      binary("assets/sky-moon.bin"),
      binary("assets/sky-earth.bin"),
      imageData("assets/earth.jpg"),
      imageData("assets/moon.jpg"),
      imageData("assets/moon-height.png"),
    ]);
    moonAtlas = new SkyAtlas(assets.sky.moon, moonBin);
    earthAtlas = new SkyAtlas(assets.sky.earth, earthBin);
    earthTexture = et;
    moonTexture = mt;
    heightMap = hm;
    worker = new Worker(new URL("./explorer-worker.mjs", import.meta.url), { type: "module" });
    worker.postMessage({ type: "init", astronomy, transfer });
    worker.onmessage = ({ data }) => {
      const p = pending.get(data.id);
      if (!p) return;
      pending.delete(data.id);
      data.error ? p.reject(Error(data.error)) : p.resolve(data.result);
    };
    worker.onerror = () => showError("The calculation worker could not start. Serve the page over HTTP from one origin.");
    skyView = new SkyView($("sky"), { stars: prepareStars(stars), earthTexture });
    bindSky();
    bindTimeline();
    bindTransport();
    bindPlaces();
    bindYear();
    legend();
    credits();
    $("share").addEventListener("click", share);
    $("ics").addEventListener("click", exportCalendar);
    $("sun-angle").addEventListener("input", renderCompare);
    $("month-prev").addEventListener("click", () => shiftMonth(-1));
    $("month-next").addEventListener("click", () => shiftMonth(1));
    renderCompare();
    let place = null;
    if (query.has("place")) place = allPlaces().find((p) => p.id === query.get("place"));
    if (!place && query.has("lat") && query.has("lon")) {
      const lat = Number(query.get("lat")),
        lon = Number(query.get("lon"));
      if (Number.isFinite(lat) && Number.isFinite(lon) && Math.abs(lat) <= 90 && Math.abs(lon) <= 180)
        place = allPlaces().find((p) => Math.abs(p.latitude - lat) < 1e-4 && Math.abs(p.longitude - lon) < 1e-4) ?? customPlace(lon, lat);
    }
    place ??= allPlaces().find((p) => p.id === "coast-1") ?? allPlaces()[0];
    globe.lon0 = place.longitude;
    globe.lat0 = clamp(place.latitude, -60, 60);
    choosePlace(place);
  } catch (error) {
    console.error(error);
    showError(`The sky could not load. ${error.message}`);
    $("headline").textContent = "The sky could not load.";
  }
}
function shiftMonth(k) {
  const [y, m] = app.viewMonth;
  const d = new Date(Date.UTC(y, m + k, 1));
  if (d.getUTCFullYear() < 2000 || d.getUTCFullYear() > 2500) return;
  app.viewMonth = [d.getUTCFullYear(), d.getUTCMonth()];
  requestMonth();
}
start();
