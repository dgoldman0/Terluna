// Open Moon Skies: the sky at one place on the Open Moon, now and through the
// coming month. Light and positions come from the shared calendar evaluator.
import { stateAt, DAY_MS } from "./model/evaluator.mjs";
import { skyFrame } from "./sky-frame.mjs";
import { SkyAtlas, SkyView, sceneFrom, earthImage, ToneCurve, prepareStars } from "./sky-render.mjs";
import { CONDITIONS, MOMENTS, conditionOf, sentence, earthLine, brightness, fromNow } from "./words.mjs";

const $ = (id) => document.getElementById(id);
const HOUR = 3600000;
const MONTH = 30 * DAY_MS;
const PLAY_HOURS_PER_SECOND = 8;
const minMs = Date.UTC(2000, 0, 1),
  maxMs = Date.UTC(2501, 0, 1) - 60000;
const clamp = (x, a, b) => Math.max(a, Math.min(b, x));
const wrap = (d) => ((((d + 180) % 360) + 360) % 360) - 180;

// The foreground: open sea at coasts and seas, level land inland.
const SEA = new Uint8Array(72).fill(1),
  LAND = new Uint8Array(72);
// `ms` is the moment shown, `base` the start of the month bar, and `live` marks
// the real present moment.
const app = { base: Date.now(), ms: Date.now(), live: true, place: null, span: null, playing: false, userLooked: false };
let astronomy, transfer, assets, places, earthTexture, skyView, worker;

// ---------- Loading ----------
const json = (url) =>
  fetch(url).then((r) => {
    if (!r.ok) throw Error(`${url} is missing (${r.status})`);
    return r.json();
  });
// The Earth map and smaller copies of it, so small disks sample it smoothly.
async function textureLevels(url) {
  const img = new Image();
  img.src = url;
  await img.decode();
  const levels = [];
  for (let w = img.naturalWidth; w >= 64; w = Math.round(w / 2)) {
    const c = document.createElement("canvas");
    c.width = w;
    c.height = Math.round((w * img.naturalHeight) / img.naturalWidth);
    const ctx = c.getContext("2d", { willReadFrequently: true });
    ctx.imageSmoothingQuality = "high";
    ctx.drawImage(img, 0, 0, c.width, c.height);
    levels.push(ctx.getImageData(0, 0, c.width, c.height));
  }
  return levels;
}

// ---------- Dates in the visitor's own time ----------
const dayFormat = new Intl.DateTimeFormat(undefined, { weekday: "short", day: "numeric", month: "short" });
const timeFormat = new Intl.DateTimeFormat(undefined, { hour: "numeric", minute: "2-digit" });
const when = (ms) => `${dayFormat.format(ms)}, ${timeFormat.format(ms)}`;

// ---------- Icons ----------
const ICON = {
  sunrise: `<svg viewBox="0 0 32 32"><path d="M4 23h24" stroke="#ffd166" stroke-width="2.4" stroke-linecap="round"/><path d="M9 23a7 7 0 0 1 14 0" fill="#ffd166"/><path d="M16 3.5v6M12.6 6.8 16 3.5l3.4 3.3" fill="none" stroke="#ffd166" stroke-width="2.4" stroke-linecap="round" stroke-linejoin="round"/></svg>`,
  sunset: `<svg viewBox="0 0 32 32"><path d="M4 23h24" stroke="#ff9f5a" stroke-width="2.4" stroke-linecap="round"/><path d="M9 23a7 7 0 0 1 14 0" fill="#ff9f5a"/><path d="M16 3.5v6M12.6 6.2 16 9.5l3.4-3.3" fill="none" stroke="#ff9f5a" stroke-width="2.4" stroke-linecap="round" stroke-linejoin="round"/></svg>`,
  dusk: `<svg viewBox="0 0 32 32"><path d="M4 24h24M8 19h16M12 14h8" stroke="#b9a5ff" stroke-width="2.4" stroke-linecap="round"/><circle cx="23" cy="7" r="1.6" fill="#b9a5ff"/><circle cx="8" cy="9" r="1.2" fill="#b9a5ff"/></svg>`,
  dawn: `<svg viewBox="0 0 32 32"><path d="M4 24h24M8 19h16M12 14h8" stroke="#ffb27a" stroke-width="2.4" stroke-linecap="round"/><path d="M16 4v5" stroke="#ffb27a" stroke-width="2.4" stroke-linecap="round"/></svg>`,
  earthrise: `<svg viewBox="0 0 32 32"><path d="M4 23h24" stroke="#6fd3ff" stroke-width="2.4" stroke-linecap="round"/><path d="M9.5 23a6.5 6.5 0 0 1 13 0" fill="#6fd3ff"/><path d="M16 4v5.5M13 7 16 4l3 3" fill="none" stroke="#6fd3ff" stroke-width="2.4" stroke-linecap="round" stroke-linejoin="round"/></svg>`,
  earthset: `<svg viewBox="0 0 32 32"><path d="M4 23h24" stroke="#6fd3ff" stroke-width="2.4" stroke-linecap="round"/><path d="M9.5 23a6.5 6.5 0 0 1 13 0" fill="#6fd3ff"/><path d="M16 4v5.5M13 6.5 16 9.5l3-3" fill="none" stroke="#6fd3ff" stroke-width="2.4" stroke-linecap="round" stroke-linejoin="round"/></svg>`,
  "full-earth": `<svg viewBox="0 0 32 32"><circle cx="16" cy="16" r="10" fill="#6fd3ff"/><path d="M10 13c3-2 5 2 8 0s4-1 5 1M9 19c2 1 4-1 6 1" stroke="#0b2a44" stroke-width="1.6" fill="none" opacity=".45"/></svg>`,
  "new-earth": `<svg viewBox="0 0 32 32"><circle cx="16" cy="16" r="9.6" fill="none" stroke="#9fb7d6" stroke-width="2"/><path d="M16 6.4a9.6 9.6 0 0 1 0 19.2" fill="none" stroke="#9fb7d6" stroke-width="1.2" opacity=".55"/></svg>`,
  "eclipse-begins": `<svg viewBox="0 0 32 32"><circle cx="16" cy="16" r="11" fill="#ff7a59" opacity=".45"/><circle cx="16" cy="16" r="8.6" fill="#141a33" stroke="#ff7a59" stroke-width="2.2"/></svg>`,
};
const SHOWN = ["sunrise", "sunset", "dusk", "dawn", "earthrise", "earthset", "full-earth", "new-earth", "eclipse-begins"];

// ---------- Calculations ----------
let requestId = 0;
const pending = new Map();
function ask(payload) {
  const id = ++requestId;
  worker.postMessage({ id, type: "span", ...payload });
  return new Promise((resolve, reject) => pending.set(id, { resolve, reject }));
}
// Events come from a window around the month bar and the moment shown, and are
// recomputed whenever the moment moves near the window's edges.
let spanLoading = false;
function ensureSpan() {
  const p = app.place;
  if (!p || spanLoading) return;
  const from = Math.min(app.base, app.ms),
    to = Math.max(app.base + MONTH, app.ms);
  const needStart = Math.max(minMs, from - 15 * DAY_MS),
    needEnd = Math.min(maxMs, to + 35 * DAY_MS);
  const sp = app.span;
  if (sp && sp.placeId === p.id && sp.start <= needStart && sp.end >= needEnd) return;
  spanLoading = true;
  ask({
    settings: { longitude: p.longitude, latitude: p.latitude, earth: true, scale: "UTC" },
    // Whole-hour samples give the same events however the page was reached.
    start: Math.ceil(Math.max(minMs, from - 25 * DAY_MS) / HOUR) * HOUR,
    end: Math.min(maxMs, to + 45 * DAY_MS),
    step: HOUR,
    order: 8,
  })
    .then((result) => {
      if (p !== app.place) return;
      result.placeId = p.id;
      app.span = result;
      const el = result.samples.earthEl;
      app.reach = Math.min(...el) > 0.5 ? "always" : Math.max(...el) < -1 ? "never" : "sometimes";
      drawTrack();
      refresh(true);
    })
    .catch(fail)
    .finally(() => {
      spanLoading = false;
      ensureSpan();
    });
}
const events = () => app.span?.events ?? [];
const last = (type) => events().filter((e) => e.type === type && e.ms <= app.ms).at(-1)?.ms;
const next = (type) => events().find((e) => e.type === type && e.ms > app.ms)?.ms;

// ---------- The sky ----------
// Look toward what matters now: the Sun or its glow by day and in twilight,
// Earth at night. The horizon sits low in the picture when the target allows.
function targetView(s) {
  const r = skyView.canvas.getBoundingClientRect();
  const W = r.width,
    H = r.height,
    portrait = W / H < 1;
  const hfov = portrait ? clamp(46 * (W / H) / 0.5, 38, 60) : 72;
  const px = W / hfov;
  const c = conditionOf(s);
  const earthUp = s.earth_elevation_deg > 1;
  let el, az;
  if (c === 0 && s.sun_elevation_deg < 32) [el, az] = [s.sun_elevation_deg, s.sun_azimuth_deg];
  else if ((c === 0 || c >= 2) && earthUp) [el, az] = [s.earth_elevation_deg, s.earth_azimuth_deg];
  else if (c === 0) [el, az] = [Math.min(s.sun_elevation_deg, 40), s.sun_azimuth_deg];
  else [el, az] = [6, s.sun_azimuth_deg];
  // Keep the horizon and the target above the time controls and clear of the text.
  const hy = portrait ? 0.64 : 0.6,
    ty = portrait ? 0.52 : 0.26,
    tx = portrait ? 0.5 : 0.7;
  let pitch = Math.max(((hy - 0.5) * H) / px, el - ((0.5 - ty) * H) / px);
  pitch = Math.min(pitch, 90 - H / 2 / px);
  return { hfov, pitch, yaw: az - ((tx - 0.5) * W) / px };
}
function aim(smooth) {
  if (app.userLooked) return;
  const t = targetView(app.state),
    cam = skyView.camera;
  if (!smooth) {
    Object.assign(cam, t);
    return;
  }
  const k = 0.12;
  cam.yaw += wrap(t.yaw - cam.yaw) * k;
  cam.pitch += (t.pitch - cam.pitch) * k;
  cam.hfov += (t.hfov - cam.hfov) * k;
}
let drawing = false,
  lastText = 0;
function refresh(force = false) {
  if (drawing && !force) return;
  drawing = true;
  requestAnimationFrame(() => {
    drawing = false;
    draw(force);
  });
}
function draw(force) {
  const p = app.place;
  const s = stateAt(app.ms, p.longitude, p.latitude, astronomy, transfer, true, "UTC");
  app.state = s;
  const frame = skyFrame(app.ms, "UTC", p.longitude, p.latitude, astronomy);
  app.scene = sceneFrom(s, frame, skyView.atlas, p.setting === "land" ? LAND : SEA);
  aim(app.playing);
  skyView.scale = app.playing ? 6 : 2.5;
  skyView.fast = app.playing;
  skyView.render(app.scene);
  const t = performance.now();
  if (force || !app.playing || t - lastText > 250) {
    lastText = t;
    words();
    adaptInk();
  }
  drawCompass();
  if (app.playing) {
    $("scrub").value = String(Math.round(((app.ms - app.base) / MONTH) * 1000));
    $("when").textContent = when(app.ms);
    if (t - lastUrlWrite > 1000) {
      ensureSpan();
      writeUrl();
    }
  }
}

// Dark ink over bright skies, light ink otherwise: once for the text block and
// once for the footer line along the bottom edge.
function brightnessOf(img, x0, x1, y0, y1) {
  const w = img.width,
    d = img.data;
  let sum = 0,
    n = 0;
  for (let j = Math.floor(img.height * y0); j < img.height * y1; j += 2)
    for (let i = Math.floor(w * x0); i < w * x1; i += 2) {
      const k = (j * w + i) * 4;
      sum += 0.2126 * d[k] + 0.7152 * d[k + 1] + 0.0722 * d[k + 2];
      n++;
    }
  return n ? sum / n / 255 : 0;
}
function adaptInk() {
  const img = skyView.lowImage;
  if (!img) return;
  const toggle = (el, cls, mean) => {
    const on = el.classList.contains(cls);
    if (!on && mean > 0.66) el.classList.add(cls);
    else if (on && mean < 0.58) el.classList.remove(cls);
  };
  toggle(document.body, "bright", brightnessOf(img, 0, 0.5, 0.12, 0.6));
  toggle($("foot"), "on-bright", brightnessOf(img, 0, 1, 0.95, 1));
}

// ---------- Words ----------
function words() {
  const s = app.state,
    c = conditionOf(s),
    now = app.ms;
  $("condition").textContent = CONDITIONS[c].name;
  const p = app.place;
  const later = stateAt(now + HOUR, p.longitude, p.latitude, astronomy, transfer, false, "UTC");
  $("sentence").textContent = sentence({
    state: s,
    condition: c,
    now,
    rising: later.sun_elevation_deg > s.sun_elevation_deg,
    lastSunset: last("sunset"),
    lastSunrise: last("sunrise"),
    nextSunrise: next("sunrise"),
    nextSunset: next("sunset"),
    reach: app.reach,
  });
  $("earth-text").textContent = earthLine(s, app.reach, next("earthrise"), now);
  $("brightness").textContent = brightness(s.total, c, assets.earth_control);
  paintEarth();
  moments();
  $("back").hidden = app.live;
  $("when").textContent = app.live ? `Now · ${when(app.ms)}` : when(app.ms);
  document.title = `${CONDITIONS[c].name} · ${app.place.name} · Open Moon Skies`;
}
const earthTone = new ToneCurve(1);
function paintEarth() {
  const c = $("earth"),
    n = c.width;
  const hide = app.reach === "never";
  c.style.display = hide ? "none" : "";
  if (hide) return;
  earthTone.sigma = 0.55;
  earthTone.contrast = 1;
  earthTone.saturation = 1;
  earthTone.night = 0;
  const img = earthImage(app.scene, earthTexture, n, earthTone, { exposure: 4, noTint: true });
  const ctx = c.getContext("2d");
  ctx.clearRect(0, 0, n, n);
  ctx.fillStyle = "rgba(20,28,52,0.85)";
  ctx.beginPath();
  ctx.arc(n / 2, n / 2, n / 2 - 1, 0, Math.PI * 2);
  ctx.fill();
  const tmp = document.createElement("canvas");
  tmp.width = tmp.height = n;
  tmp.getContext("2d").putImageData(new ImageData(img.data, n, n), 0, 0);
  ctx.drawImage(tmp, 0, 0);
}
let momentsKey = "";
function moments() {
  const list = events().filter((e) => SHOWN.includes(e.type) && e.ms > app.ms).slice(0, 6);
  const key = `${list.map((e) => e.ms).join()}|${Math.floor(app.ms / HOUR)}|${app.span?.end ?? 0}`;
  if (key === momentsKey) return;
  momentsKey = key;
  const ol = $("moments");
  ol.replaceChildren();
  if (!list.length && app.span && app.span.end >= maxMs - HOUR) {
    const li = document.createElement("li");
    li.className = "end-note";
    li.textContent = "The calendar ends with the year 2500.";
    ol.append(li);
  }
  for (const e of list) {
    const li = document.createElement("li"),
      b = document.createElement("button");
    b.type = "button";
    b.innerHTML = ICON[e.type];
    const title = document.createElement("strong");
    title.textContent = MOMENTS[e.type];
    const date = document.createElement("span");
    date.className = "date";
    const day = document.createElement("span"),
      time = document.createElement("span");
    day.textContent = dayFormat.format(e.ms);
    time.textContent = timeFormat.format(e.ms);
    date.append(day, " ", time);
    const rel = document.createElement("span");
    rel.className = "in";
    rel.textContent = fromNow(e.ms - app.ms);
    b.append(title, date, rel);
    b.setAttribute("aria-label", `${MOMENTS[e.type]}, ${when(e.ms)}, ${fromNow(e.ms - app.ms)}. Show this moment.`);
    b.addEventListener("click", () => {
      stop();
      // Show each crossing from inside the brighter side: sunrise and sunset two
      // hours into the day, first light and fading twilight an hour into the twilight.
      const shift = { sunrise: 2 * HOUR, sunset: -2 * HOUR, dawn: HOUR, dusk: -HOUR }[e.type] ?? 0;
      setTime(e.ms + shift, { reaim: true });
    });
    li.append(b);
    ol.append(li);
  }
}

// ---------- The month bar ----------
function drawTrack() {
  const sp = app.span;
  if (!sp) {
    $("track-colours").style.background = "";
    return;
  }
  const stops = [];
  for (let k = 0; k <= 60; k++) {
    const i = Math.round((app.base + (k / 60) * MONTH - sp.start) / sp.step);
    const colour = i >= 0 && i < sp.samples.ms.length ? CONDITIONS[sp.samples.regime[i]].colour : "rgba(255,255,255,0.3)";
    stops.push(`${colour} ${((k / 60) * 100).toFixed(1)}%`);
  }
  $("track-colours").style.background = `linear-gradient(90deg, ${stops.join(",")})`;
}
// Show a moment. A moment outside the month bar moves the bar to start two days before it.
function setTime(ms, { reaim = false, live = false } = {}) {
  // Whole minutes, as the address records them.
  app.ms = clamp(live ? ms : Math.round(ms / 60000) * 60000, minMs, maxMs);
  app.live = live;
  if (app.ms < app.base || app.ms > app.base + MONTH) {
    app.base = clamp(app.ms - 2 * DAY_MS, minMs, maxMs - MONTH);
    drawTrack();
  }
  $("scrub").value = String(Math.round(((app.ms - app.base) / MONTH) * 1000));
  if (reaim) {
    app.userLooked = false;
    const p = app.place;
    app.state = stateAt(app.ms, p.longitude, p.latitude, astronomy, transfer, true, "UTC");
    aim(false);
  }
  refresh(true);
  ensureSpan();
  writeUrl();
}

// The address records the place, the moment unless it is the present, and the
// view once the sky has been dragged; opening the address restores them.
let lastUrlWrite = 0;
function writeUrl() {
  if (!app.place) return;
  const parts = [`place=${app.place.id}`];
  const iso = new Date(Math.round(app.ms / 60000) * 60000).toISOString().replace(":00.000Z", "Z");
  if (!app.live) parts.push(`at=${iso}`);
  if (app.userLooked) {
    const c = skyView.camera;
    parts.push(`view=${Math.round((((c.yaw % 360) + 360) % 360))},${Math.round(c.pitch)}`);
  }
  const url = `${location.pathname}?${parts.join("&")}`;
  if (url !== location.pathname + location.search) history.replaceState(null, "", url);
  lastUrlWrite = performance.now();
  $("almanac-link").href = `almanac.html?lat=${app.place.latitude}&lon=${app.place.longitude}&at=${iso}`;
}

// ---------- Directions ----------
// Compass letters along a strip just above the controls, each under the part
// of the sky it names; narrow views add the points between them.
const COMPASS = ["N", "NNE", "NE", "ENE", "E", "ESE", "SE", "SSE", "S", "SSW", "SW", "WSW", "W", "WNW", "NW", "NNW"];
let compassChips = null;
function drawCompass() {
  const box = $("compass");
  if (!compassChips)
    compassChips = COMPASS.map((text) => {
      const chip = document.createElement("span");
      chip.textContent = text;
      box.append(chip);
      return chip;
    });
  const c = skyView.canvas,
    W = c.width,
    H = c.height,
    d = skyView.dpr,
    cam = skyView.camera;
  const fine = cam.hfov < 60;
  const y = document.querySelector(".dock").getBoundingClientRect().top - 26;
  const text = document.querySelector(".now").getBoundingClientRect();
  COMPASS.forEach((_, k) => {
    const az = (k * 22.5 * Math.PI) / 180;
    const p = cam.project([Math.sin(az), Math.cos(az), 0], W, H);
    const x = p ? p[0] / d : -1;
    const clear = !(x > text.left - 24 && x < text.right + 24 && y - 12 < text.bottom && y + 12 > text.top);
    const show = (fine || k % 2 === 0) && p && x > 20 && x < W / d - 20 && clear;
    compassChips[k].style.display = show ? "" : "none";
    if (show) {
      compassChips[k].style.left = `${x}px`;
      compassChips[k].style.top = `${y}px`;
    }
  });
}

// ---------- Choosing a date ----------
const pad = (n) => String(n).padStart(2, "0");
// A datetime-local value in the visitor's time zone.
function localValue(ms) {
  const d = new Date(ms);
  return `${String(d.getFullYear()).padStart(4, "0")}-${pad(d.getMonth() + 1)}-${pad(d.getDate())}T${pad(d.getHours())}:${pad(d.getMinutes())}`;
}
function openPicker() {
  const input = $("when-input");
  input.min = localValue(minMs);
  input.max = localValue(maxMs);
  input.value = localValue(app.ms);
  $("picker").hidden = false;
  $("when-button").setAttribute("aria-expanded", "true");
  input.focus();
  try {
    input.showPicker?.();
  } catch {
    // The field stays open for typing where the browser has no picker to show.
  }
}
function closePicker() {
  $("picker").hidden = true;
  $("when-button").setAttribute("aria-expanded", "false");
}
function backToNow() {
  stop();
  app.base = Date.now();
  setTime(Date.now(), { reaim: true, live: true });
}

// ---------- Play ----------
let lastTick = 0;
function tick(t) {
  if (!app.playing) return;
  const dt = lastTick ? Math.min(t - lastTick, 80) : 16;
  lastTick = t;
  let ms = app.ms + (dt / 1000) * PLAY_HOURS_PER_SECOND * HOUR;
  if (ms > app.base + MONTH) ms = app.base;
  app.ms = ms;
  draw(false);
  requestAnimationFrame(tick);
}
function play() {
  app.playing = true;
  app.live = false;
  app.userLooked = false;
  lastTick = 0;
  $("play").classList.add("playing");
  $("play").setAttribute("aria-label", "Pause");
  requestAnimationFrame(tick);
}
function stop() {
  if (!app.playing) return;
  app.playing = false;
  app.ms = Math.round(app.ms / 60000) * 60000;
  $("play").classList.remove("playing");
  $("play").setAttribute("aria-label", "Play the coming month");
  refresh(true);
  writeUrl();
}

// ---------- Places ----------
function choosePlace(p, { view = null } = {}) {
  app.place = p;
  app.span = null;
  app.reach = null;
  momentsKey = "";
  drawTrack();
  $("place-name").textContent = p.name;
  $("place-latin").textContent = `${p.latin} · ${p.line}`;
  if (!view) {
    setTime(app.ms, { reaim: true, live: app.live });
    return;
  }
  app.userLooked = false;
  app.state = stateAt(app.ms, p.longitude, p.latitude, astronomy, transfer, true, "UTC");
  aim(false);
  skyView.camera.yaw = view.yaw;
  skyView.camera.pitch = view.pitch;
  app.userLooked = true;
  setTime(app.ms, { live: app.live });
}
function placeList() {
  const ul = $("place-list");
  ul.replaceChildren();
  for (const p of places) {
    const li = document.createElement("li"),
      b = document.createElement("button");
    b.value = p.id;
    b.setAttribute("aria-current", String(p === app.place));
    const strong = document.createElement("strong"),
      span = document.createElement("span");
    strong.textContent = p.name;
    span.textContent = p.line;
    b.append(strong, span);
    b.addEventListener("click", (e) => {
      e.preventDefault();
      $("places").close();
      if (p !== app.place) choosePlace(p);
    });
    li.append(b);
    ul.append(li);
  }
}

// ---------- Looking around ----------
function bindSky() {
  const c = $("sky");
  let drag = null;
  c.addEventListener("pointerdown", (e) => {
    c.setPointerCapture(e.pointerId);
    drag = { x: e.clientX, y: e.clientY, yaw: skyView.camera.yaw, pitch: skyView.camera.pitch };
    c.classList.add("dragging");
  });
  c.addEventListener("pointermove", (e) => {
    if (!drag) return;
    const r = c.getBoundingClientRect(),
      cam = skyView.camera,
      k = cam.hfov / r.width;
    app.userLooked = true;
    cam.yaw = drag.yaw - (e.clientX - drag.x) * k;
    cam.pitch = clamp(drag.pitch + (e.clientY - drag.y) * k, -15, 90 - (r.height / 2) * k);
    refresh();
  });
  const end = () => {
    if (drag && app.userLooked) writeUrl();
    drag = null;
    c.classList.remove("dragging");
  };
  c.addEventListener("pointerup", end);
  c.addEventListener("pointercancel", end);
  new ResizeObserver(() => {
    if (!app.state) return;
    aim(false);
    refresh(true);
  }).observe(c);
}

function fail(error) {
  console.error(error);
  const m = $("message");
  m.hidden = false;
  m.textContent = `The sky could not load: ${error.message}`;
  $("condition").textContent = "";
}

async function start() {
  if (location.protocol === "file:") return;
  try {
    const query = new URLSearchParams(location.search);
    const at = Date.parse(query.get("at") ?? "");
    if (Number.isFinite(at)) {
      app.ms = clamp(at, minMs, maxMs);
      app.live = false;
      app.base = clamp(app.ms - 2 * DAY_MS, minMs, maxMs - MONTH);
    }
    const look = (query.get("view") ?? "").split(",").map(Number);
    const view = look.length === 2 && look.every(Number.isFinite) ? { yaw: look[0], pitch: clamp(look[1], -15, 89) } : null;
    let stars, sky, placeData;
    [astronomy, transfer, assets, stars, placeData, sky, earthTexture] = await Promise.all([
      json("data/astronomy.json"),
      json("data/transfer.json"),
      json("assets/assets.json"),
      json("assets/stars.json"),
      json("assets/places.json"),
      fetch("assets/sky-moon.bin").then((r) => {
        if (!r.ok) throw Error(`assets/sky-moon.bin is missing (${r.status})`);
        return r.arrayBuffer();
      }),
      textureLevels("assets/earth.jpg"),
    ]);
    places = placeData.places;
    skyView = new SkyView($("sky"), { stars: prepareStars(stars), earthTexture });
    skyView.atlas = new SkyAtlas(assets.sky, sky);
    worker = new Worker(new URL("./explorer-worker.mjs", import.meta.url), { type: "module" });
    worker.postMessage({ type: "init", astronomy, transfer });
    worker.onmessage = ({ data }) => {
      const p = pending.get(data.id);
      if (!p) return;
      pending.delete(data.id);
      data.error ? p.reject(Error(data.error)) : p.resolve(data.result);
    };
    worker.onerror = () => fail(Error("the calculation worker could not start"));
    bindSky();
    $("play").addEventListener("click", () => (app.playing ? stop() : play()));
    $("scrub").addEventListener("input", () => {
      stop();
      setTime(app.base + (Number($("scrub").value) / 1000) * MONTH);
    });
    $("back").addEventListener("click", backToNow);
    $("when-button").addEventListener("click", () => ($("picker").hidden ? openPicker() : closePicker()));
    $("when-input").addEventListener("change", () => {
      const t = new Date($("when-input").value).getTime();
      if (!Number.isFinite(t)) return;
      stop();
      setTime(clamp(t, minMs, maxMs), { reaim: true });
    });
    $("picker-now").addEventListener("click", () => {
      closePicker();
      backToNow();
    });
    $("picker-done").addEventListener("click", closePicker);
    // Captured first: the browser's own date field keeps Escape for itself.
    document.addEventListener(
      "keydown",
      (e) => {
        if (e.key === "Escape" && !$("picker").hidden) {
          closePicker();
          $("when-button").focus();
        }
      },
      true,
    );
    document.addEventListener("pointerdown", (e) => {
      if (!$("picker").hidden && !e.target.closest(".when")) closePicker();
    });
    // The present moment moves on while the page is open.
    setInterval(() => {
      if (app.live && !app.playing && document.visibilityState === "visible") setTime(Date.now(), { live: true });
    }, 60000);
    $("place").addEventListener("click", () => {
      placeList();
      $("places").showModal();
    });
    document.body.classList.remove("loading");
    choosePlace(places.find((p) => p.id === query.get("place")) ?? places[0], { view });
  } catch (error) {
    fail(error);
  }
}
start();
