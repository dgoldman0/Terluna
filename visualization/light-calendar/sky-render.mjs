// Display renderer for the explorer: the sky from the solved spherical sky
// atlas, the Sun and Earth from the calendar state, stars from the Yale Bright
// Star Catalogue, and a level sea or land foreground. Absolute luminances
// (cd m^-2) go through one display tone curve; the tone curve, the eye's
// low-light desaturation, the sea's glitter and the ground's colour are display
// approximations.
import { apply, localVector, lonLat, starPrecession, mul } from "./sky-frame.mjs";

const rad = Math.PI / 180;
const clamp = (x, a, b) => Math.max(a, Math.min(b, x));
const STAR_PSF_SR = (1.5 * rad / 60) ** 2; // a point source spread over ~1.5 arcmin
export const RAYLEIGH_TAU_550 = 0.757; // solved 1.2-atm column (illumination/sky)
const GROUND_ALBEDO = 0.1; // the atlas and calendar ground albedo
const SEA_SLOPE_VARIANCE = 0.02; // Cox–Munk mean square slope for a light breeze (display)

function bracket(a, x) {
  if (x <= a[0]) return [0, 0];
  if (x >= a[a.length - 1]) return [a.length - 2, 1];
  let lo = 0,
    hi = a.length - 1;
  while (hi - lo > 1) {
    const m = (lo + hi) >> 1;
    if (a[m] <= x) lo = m;
    else hi = m;
  }
  return [lo, (x - a[lo]) / (a[lo + 1] - a[lo])];
}
const lerp = (a, b, t) => a + (b - a) * t;

// ---------- Sky atlas ----------
export class SkyAtlas {
  constructor(header, buffer) {
    this.h = header;
    this.data = new Uint16Array(buffer);
    this.ns = header.suns.length;
    this.ne = header.elevations.length;
    this.na = header.azimuths.length;
    this.daz = header.azimuths[1] - header.azimuths[0];
    const [lo, hi] = header.log_range;
    this.logScale = (hi - lo) / 65535;
    this.logLo = lo;
    // Sin-elevation index table for quick lookups.
    this.sinEl = header.elevations.map((e) => Math.sin(e * rad));
    this.cache = new Map();
  }
  // The [elevation][azimuth] slice for one Sun elevation: log10 Y, x, y.
  slice(sunEl) {
    const key = Math.round(sunEl * 200) / 200;
    if (this.cache.has(key)) return this.cache.get(key);
    const [i, w] = bracket(this.h.suns, Math.max(sunEl, this.h.suns[0]));
    const n = this.ne * this.na,
      out = new Float32Array(n * 3);
    for (let k = 0; k < n; k++) {
      const a = (i * n + k) * 3,
        b = ((i + 1) * n + k) * 3;
      out[k * 3] = lerp(this.data[a], this.data[b], w) * this.logScale + this.logLo;
      out[k * 3 + 1] = lerp(this.data[a + 1], this.data[b + 1], w) / 65535;
      out[k * 3 + 2] = lerp(this.data[a + 2], this.data[b + 2], w) / 65535;
    }
    if (this.cache.size > 64) this.cache.clear();
    this.cache.set(key, out);
    return out;
  }
  diffuseLux(sunEl) {
    const [i, w] = bracket(this.h.suns, Math.max(sunEl, this.h.suns[0]));
    return Math.exp(
      lerp(Math.log(this.h.diffuse_lux[i]), Math.log(this.h.diffuse_lux[i + 1]), w),
    );
  }
  beamXY(sunEl) {
    const [i, w] = bracket(this.h.suns, clamp(sunEl, 0, 90));
    const a = this.h.direct_xy[i],
      b = this.h.direct_xy[i + 1];
    return [lerp(a[0], b[0], w), lerp(a[1], b[1], w)];
  }
  // XYZ (cd m^-2) at view elevation/azimuth-from-source on a slice, times scale.
  sample(slice, el, relAz, scale, out) {
    const s = Math.sin(Math.max(0, el) * rad);
    const [j, v] = bracket(this.sinEl, s);
    let az = Math.abs(relAz) % 360;
    if (az > 180) az = 360 - az;
    const fa = az / this.daz,
      ia = Math.min(this.na - 2, Math.floor(fa)),
      u = fa - ia;
    const na = this.na;
    const p00 = (j * na + ia) * 3,
      p01 = p00 + 3,
      p10 = p00 + na * 3,
      p11 = p10 + 3;
    const w00 = (1 - v) * (1 - u),
      w01 = (1 - v) * u,
      w10 = v * (1 - u),
      w11 = v * u;
    const logY =
      slice[p00] * w00 + slice[p01] * w01 + slice[p10] * w10 + slice[p11] * w11;
    const x =
      slice[p00 + 1] * w00 + slice[p01 + 1] * w01 + slice[p10 + 1] * w10 + slice[p11 + 1] * w11;
    const y =
      slice[p00 + 2] * w00 + slice[p01 + 2] * w01 + slice[p10 + 2] * w10 + slice[p11 + 2] * w11;
    const Y = 10 ** logY * scale;
    out[0] += (x / y) * Y;
    out[1] += Y;
    out[2] += ((1 - x - y) / y) * Y;
  }
}

// ---------- Tone curve (display approximation) ----------
const M = [
  [3.2406, -1.5372, -0.4986],
  [-0.9689, 1.8758, 0.0415],
  [0.0557, -0.204, 1.057],
];
const encode = (c) =>
  c <= 0.0031308 ? 12.92 * c : 1.055 * Math.pow(c, 1 / 2.4) - 0.055;
export const ENCODE = new Uint8ClampedArray(4097);
for (let i = 0; i <= 4096; i++) ENCODE[i] = Math.round(255 * encode(i / 4096));

export class ToneCurve {
  // adaptLux: horizontal illuminance the eye is adapted to; boost: exposure factor.
  constructor(adaptLux, boost = 1) {
    const key = Math.max(adaptLux, 1e-7) / Math.PI;
    this.sigma = Math.max(1.12 * Math.pow(key, 0.78) / boost, 0.004);
    // Colour fades as the light falls toward rod vision.
    this.saturation = clamp((Math.log10(key) + 1.7) / 2.4, 0.18, 1);
    this.night = 1 - this.saturation;
  }
  // XYZ -> sRGB bytes into out[o..o+2].
  write(X, Y, Z, out, o) {
    if (!(Y > 0)) {
      out[o] = out[o + 1] = out[o + 2] = 0;
      return;
    }
    const v = Y / (Y + this.sigma);
    let r = (M[0][0] * X + M[0][1] * Y + M[0][2] * Z) / Y,
      g = (M[1][0] * X + M[1][1] * Y + M[1][2] * Z) / Y,
      b = (M[2][0] * X + M[2][1] * Y + M[2][2] * Z) / Y;
    const s = this.saturation;
    r = 1 + (r - 1) * s;
    g = 1 + (g - 1) * s;
    b = 1 + (b - 1) * s;
    // Rod vision reads as a cool grey-blue on screen.
    const n = this.night * 0.55;
    r *= 1 - 0.28 * n;
    g *= 1 - 0.08 * n;
    b *= 1 + 0.32 * n;
    r = Math.max(r, 0) * v;
    g = Math.max(g, 0) * v;
    b = Math.max(b, 0) * v;
    const m = Math.max(r, g, b);
    if (m > 1) {
      // Soft shoulder: blend toward white instead of hue-shifting clips.
      const t = clamp((m - 1) / m, 0, 1);
      r = r / m + (1 - r / m) * t;
      g = g / m + (1 - g / m) * t;
      b = b / m + (1 - b / m) * t;
    }
    out[o] = ENCODE[(clamp(r, 0, 1) * 4096) | 0];
    out[o + 1] = ENCODE[(clamp(g, 0, 1) * 4096) | 0];
    out[o + 2] = ENCODE[(clamp(b, 0, 1) * 4096) | 0];
  }
  value(Y) {
    return Y / (Y + this.sigma);
  }
}

export function xyToXYZ(x, y, Y) {
  return [(x / y) * Y, Y, ((1 - x - y) / y) * Y];
}

// ---------- Scene description from a calendar state ----------
export function sceneFrom(state, frame, atlas, sectors) {
  const sunEl = state.sun_elevation_deg,
    earthEl = state.earth_elevation_deg;
  const sunAtlas = atlas.diffuseLux(sunEl);
  const earthAtlas = atlas.diffuseLux(earthEl);
  const sunDir = localVector(sunEl, state.sun_azimuth_deg);
  const earthDir = localVector(earthEl, state.earth_azimuth_deg);
  // Direct beams as normal illuminance (lux) and their colour through the air.
  const sinS = Math.sin(Math.max(sunEl, 0.05) * rad),
    sinE = Math.sin(Math.max(earthEl, 0.05) * rad);
  return {
    state,
    frame,
    atlas,
    sectors,
    sunEl,
    sunAz: state.sun_azimuth_deg,
    earthEl,
    earthAz: state.earth_azimuth_deg,
    sunDir,
    earthDir,
    sunSlice: atlas.slice(sunEl),
    earthSlice: atlas.slice(Math.max(earthEl, -30)),
    sunScale: sunAtlas > 0 ? state.solar_diffuse / sunAtlas : 0,
    earthScale: earthAtlas > 0 ? state.earth_diffuse / earthAtlas : 0,
    sunNormalLux: state.solar_direct / sinS,
    earthNormalLux: state.earth_direct / sinE,
    sunBeamXY: atlas.beamXY(sunEl),
    earthBeamXY: atlas.beamXY(earthEl),
    zenithBeamXY: atlas.beamXY(90),
    tone: new ToneCurve(Math.max(state.total, 1e-6), 1.35),
  };
}

const tmp = new Float64Array(3);
// Sky XYZ for a view direction above the horizon.
export function skyXYZ(scene, el, az, out) {
  out[0] = out[1] = out[2] = 0;
  if (scene.sunScale > 0)
    scene.atlas.sample(scene.sunSlice, el, az - scene.sunAz, scene.sunScale, out);
  if (scene.earthScale > 0)
    scene.atlas.sample(scene.earthSlice, el, az - scene.earthAz, scene.earthScale, out);
  return out;
}

function fresnel(cosI) {
  const c = clamp(cosI, 0, 1);
  return 0.02 + 0.98 * Math.pow(1 - c, 5);
}

// Cox–Munk sun/Earth glint luminance for a downward view direction v.
function glint(v, source, normalLux) {
  if (source[2] <= 0 || normalLux <= 0) return 0;
  let hx = source[0] - v[0],
    hy = source[1] - v[1],
    hz = source[2] - v[2];
  const hn = Math.hypot(hx, hy, hz);
  hx /= hn;
  hy /= hn;
  hz /= hn;
  if (hz <= 0) return 0;
  const tan2 = (1 - hz * hz) / (hz * hz);
  const p = Math.exp(-tan2 / SEA_SLOPE_VARIANCE) / (Math.PI * SEA_SLOPE_VARIANCE);
  const cosV = Math.max(-v[2], 0.02);
  const cosW = source[0] * hx + source[1] * hy + source[2] * hz;
  return (fresnel(cosW) * normalLux * p) / (4 * cosV * hz ** 4);
}

const SOIL = [0.36, 0.3, 0.22]; // developed soil, linear sRGB reflectance shape (display)
const SOIL_Y = 0.2126 * SOIL[0] + 0.7152 * SOIL[1] + 0.0722 * SOIL[2];
const toXYZ = (r, g, b, out) => {
  out[0] = 0.4124 * r + 0.3576 * g + 0.1805 * b;
  out[1] = 0.2126 * r + 0.7152 * g + 0.0722 * b;
  out[2] = 0.0193 * r + 0.1192 * g + 0.9505 * b;
  return out;
};
const horizonTmp = new Float64Array(3),
  skyTmp = new Float64Array(3);

// Foreground XYZ for a view direction below the horizon.
export function groundXYZ(scene, el, az, v, out) {
  const sector = scene.sectors ? scene.sectors[((Math.round(az / 5) % 72) + 72) % 72] : 0;
  const horizon = skyXYZ(scene, 0.3, az, horizonTmp);
  const fade = Math.exp(-Math.abs(el) / 0.45);
  const E = scene.state.total;
  if (sector) {
    // Sea: wave facets tilted by about 7° (display) reflect higher, dimmer sky
    // with less than mirror reflectance; then the water body and glints.
    skyXYZ(scene, Math.min(90, -el + 12), az, out);
    const R = fresnel(Math.sin((-el + 7) * rad));
    out[0] *= R;
    out[1] *= R;
    out[2] *= R;
    const body = (0.015 * E) / Math.PI;
    out[0] += body * 0.8;
    out[1] += body;
    out[2] += body * 1.5;
    const gs = glint(v, scene.sunDir, scene.sunNormalLux);
    if (gs > 0) {
      const c = xyToXYZ(scene.sunBeamXY[0], scene.sunBeamXY[1], gs);
      out[0] += c[0];
      out[1] += c[1];
      out[2] += c[2];
    }
    const ge = glint(v, scene.earthDir, scene.earthNormalLux);
    if (ge > 0) {
      out[0] += ge * 0.92;
      out[1] += ge;
      out[2] += ge * 1.18;
    }
  } else {
    // Level land with the model's ground albedo under the sky's light colour.
    const sky = skyXYZ(scene, 60, az, skyTmp);
    const Ys = Math.max(sky[1], 1e-30);
    const X = sky[0] / Ys,
      Z = sky[2] / Ys;
    const r = Math.max(0, M[0][0] * X + M[0][1] + M[0][2] * Z) * (SOIL[0] / SOIL_Y),
      g = Math.max(0, M[1][0] * X + M[1][1] + M[1][2] * Z) * (SOIL[1] / SOIL_Y),
      b = Math.max(0, M[2][0] * X + M[2][1] + M[2][2] * Z) * (SOIL[2] / SOIL_Y);
    toXYZ(r, g, b, out);
    const k = (GROUND_ALBEDO * E) / Math.PI / Math.max(out[1], 1e-30);
    out[0] *= k;
    out[1] *= k;
    out[2] *= k;
  }
  // The distant ground fades into the horizon sky.
  out[0] = out[0] * (1 - fade) + horizon[0] * fade;
  out[1] = out[1] * (1 - fade) + horizon[1] * fade;
  out[2] = out[2] * (1 - fade) + horizon[2] * fade;
  return out;
}

// ---------- Camera ----------
// A panorama view: azimuth runs across the picture and elevation up it at one
// scale, so the horizon is straight and a wide sweep of sky fits at once. The
// Sun and Earth are drawn round at their true angular size.
export class Camera {
  constructor(yaw = 180, pitch = 25, hfov = 150) {
    this.yaw = yaw; // azimuth at the centre of the picture
    this.pitch = pitch; // elevation at the centre of the picture
    this.hfov = hfov; // degrees of azimuth across the picture
  }
  basis() {
    return null;
  }
  scale(W) {
    return W / this.hfov;
  }
  project(d, W, H) {
    const el = Math.asin(clamp(d[2], -1, 1)) / rad,
      az = Math.atan2(d[0], d[1]) / rad;
    const rel = ((((az - this.yaw) % 360) + 540) % 360) - 180;
    if (Math.abs(rel) > this.hfov / 2 + 20) return null;
    const s = this.scale(W);
    return [W / 2 + rel * s, H / 2 - (el - this.pitch) * s];
  }
  ray(px, py, W, H, b, out) {
    const s = this.scale(W);
    const el = clamp(this.pitch + (H / 2 - py) / s, -89.9, 89.9) * rad,
      az = (this.yaw + (px - W / 2) / s) * rad;
    out[0] = Math.cos(el) * Math.sin(az);
    out[1] = Math.cos(el) * Math.cos(az);
    out[2] = Math.sin(el);
    return out;
  }
  pixelsPerDegree(W) {
    return this.scale(W);
  }
  horizonY(H, W) {
    return H / 2 + this.pitch * this.scale(W);
  }
}

// ---------- Earth disk ----------
// Renders Earth's lit face, oriented with its north on the observer's sky,
// into an RGBA buffer of size n. Earth's surface luminance is scaled so the
// disk's integral equals the calendar's direct Earthlight; the cloud and
// surface pattern is a NASA Blue Marble composite.
export function earthImage(scene, texture, n, tone, options = {}) {
  // tone may be replaced by a local exposure below.
  const { state, frame } = scene;
  const e = scene.earthDir;
  // Screen basis at Earth: up toward the zenith, right = e x up.
  let up = [-e[0] * e[2], -e[1] * e[2], 1 - e[2] * e[2]];
  const ul = Math.hypot(...up);
  up = ul > 1e-9 ? up.map((x) => x / ul) : [0, 1, 0];
  const right = [
    e[1] * up[2] - e[2] * up[1],
    e[2] * up[0] - e[0] * up[2],
    e[0] * up[1] - e[1] * up[0],
  ];
  const s = scene.sunDir;
  const L = frame.localToEarth;
  const pixels = new Float32Array(n * n * 4);
  let sum = 0;
  const tw = texture.width,
    th = texture.height,
    td = texture.data;
  for (let py = 0; py < n; py++)
    for (let px = 0; px < n; px++) {
      const x = ((px + 0.5) / n) * 2 - 1,
        y = 1 - ((py + 0.5) / n) * 2,
        r2 = x * x + y * y;
      if (r2 > 1) continue;
      const z = Math.sqrt(1 - r2);
      const nl = [
        x * right[0] + y * up[0] - z * e[0],
        x * right[1] + y * up[1] - z * e[1],
        x * right[2] + y * up[2] - z * e[2],
      ];
      const illum = nl[0] * s[0] + nl[1] * s[1] + nl[2] * s[2];
      const k = (py * n + px) * 4;
      // Anti-aliased limb.
      const edge = clamp((1 - Math.sqrt(r2)) * n * 0.5, 0, 1);
      pixels[k + 3] = edge;
      if (illum <= 0) continue;
      const ef = apply(L, nl);
      const [lon, lat] = lonLat(ef);
      const tx = Math.min(tw - 1, Math.floor(((lon + 180) / 360) * tw)),
        ty = Math.min(th - 1, Math.floor(((90 - lat) / 180) * th));
      const t = (ty * tw + tx) * 4;
      const cr = (td[t] / 255) ** 2.2,
        cg = (td[t + 1] / 255) ** 2.2,
        cb = (td[t + 2] / 255) ** 2.2;
      // Soft terminator: the day side brightens over a few degrees.
      const lit = illum < 0.06 ? illum * (illum / 0.06) : illum;
      pixels[k] = cr * lit;
      pixels[k + 1] = cg * lit;
      pixels[k + 2] = cb * lit;
      sum += (0.2126 * cr + 0.7152 * cg + 0.0722 * cb) * lit;
    }
  // Scale to the calendar's direct Earthlight: sum(L dOmega) = E_normal.
  const radius = state.earth_radius_deg * rad;
  const dOmega = ((2 * radius) / n) ** 2;
  const K = options.exposure
    ? options.exposure
    : sum > 0
      ? scene.earthNormalLux / (sum * dOmega)
      : 0;
  // Air between the observer and Earth reddens its light (relative to the zenith beam).
  const tint = options.noTint ? [1, 1, 1] : beamTint(scene.earthBeamXY, scene.zenithBeamXY);
  let behind = options.behind ?? [0, 0, 0];
  let visibility = 1;
  if (options.local) {
    // Local adaptation (display): the eye resolves detail on the bright disk,
    // so Earth gets its own exposure and is added onto the sky by its contrast.
    let lit = 0;
    for (let i = 0; i < n * n; i++) if (pixels[i * 4] + pixels[i * 4 + 1] + pixels[i * 4 + 2] > 0) lit++;
    const meanLit = lit ? (K * sum) / lit : 0;
    const local = new ToneCurve(1);
    local.sigma = Math.max(0.6 * meanLit, 1e-9);
    local.saturation = 1;
    local.night = 0;
    const contrast = meanLit / Math.max(behind[1], 1e-12);
    visibility = clamp(contrast / (contrast + 0.3), 0, 1);
    tone = local;
    behind = [0, 0, 0];
  }
  const out = new Uint8ClampedArray(n * n * 4);
  const tmpOut = new Uint8ClampedArray(3);
  for (let i = 0; i < n * n; i++) {
    const a = pixels[i * 4 + 3];
    if (a <= 0) continue;
    const r = pixels[i * 4] * K * tint[0],
      g = pixels[i * 4 + 1] * K * tint[1],
      b = pixels[i * 4 + 2] * K * tint[2];
    // linear sRGB -> XYZ
    const X = 0.4124 * r + 0.3576 * g + 0.1805 * b + behind[0],
      Y = 0.2126 * r + 0.7152 * g + 0.0722 * b + behind[1],
      Z = 0.0193 * r + 0.1192 * g + 0.9505 * b + behind[2];
    tone.write(X, Y, Z, tmpOut, 0);
    out[i * 4] = tmpOut[0];
    out[i * 4 + 1] = tmpOut[1];
    out[i * 4 + 2] = tmpOut[2];
    out[i * 4 + 3] = Math.round(a * 255);
  }
  return { data: out, size: n, luminanceScale: K, up, right, visibility };
}

function xyToRGB(xy) {
  const [X, Y, Z] = xyToXYZ(xy[0], xy[1], 1);
  return [
    M[0][0] * X + M[0][1] * Y + M[0][2] * Z,
    M[1][0] * X + M[1][1] * Y + M[1][2] * Z,
    M[2][0] * X + M[2][1] * Y + M[2][2] * Z,
  ];
}
function beamTint(xy, ref) {
  if (!xy || !(xy[1] > 0)) return [1, 1, 1];
  const a = xyToRGB(xy),
    b = xyToRGB(ref);
  const t = a.map((v, i) => Math.max(v, 0) / Math.max(b[i], 1e-6));
  const y = 0.2126 * t[0] + 0.7152 * t[1] + 0.0722 * t[2];
  return t.map((v) => v / y);
}

// ---------- Stars ----------
export function prepareStars(catalogue) {
  const stars = [];
  for (const [ra, dec, mag, bv] of catalogue.stars)
    stars.push({
      v: [Math.cos(dec) * Math.cos(ra), Math.cos(dec) * Math.sin(ra), Math.sin(dec)],
      mag,
      colour: starColour(bv),
      lux: 2.54e-6 * 10 ** (-0.4 * mag),
    });
  return stars;
}
function starColour(bv) {
  // Ballesteros (2012) temperature, then a gentle blackbody-like tint.
  const t = 4600 * (1 / (0.92 * bv + 1.7) + 1 / (0.92 * bv + 0.62));
  const k = clamp((t - 3000) / 9000, 0, 1);
  return [1, 0.78 + 0.2 * k, 0.58 + 0.42 * k];
}
export function airmass(el) {
  const s = Math.sin(Math.max(el, 0) * rad);
  return 1 / (s + 0.133 * Math.exp(-s / 0.1));
}

// ---------- The sky window ----------
export class SkyView {
  constructor(canvas, assets) {
    this.canvas = canvas;
    this.ctx = canvas.getContext("2d");
    this.assets = assets;
    this.camera = new Camera();
    this.low = document.createElement("canvas");
    this.lowCtx = this.low.getContext("2d");
    this.earthCanvas = document.createElement("canvas");
    this.scale = 3;
    this.labels = true;
  }
  resize() {
    const r = this.canvas.getBoundingClientRect();
    const dpr = Math.min(window.devicePixelRatio || 1, 2);
    const W = Math.max(1, Math.round(r.width * dpr)),
      H = Math.max(1, Math.round(r.height * dpr));
    if (this.canvas.width !== W || this.canvas.height !== H) {
      this.canvas.width = W;
      this.canvas.height = H;
    }
    this.dpr = dpr;
    const step = Math.max(2, Math.round(this.scale * dpr));
    const w = Math.ceil(W / step),
      h = Math.ceil(H / step);
    if (this.low.width !== w || this.low.height !== h) {
      this.low.width = w;
      this.low.height = h;
      this.lowImage = this.lowCtx.createImageData(w, h);
    }
  }
  render(scene) {
    this.scene = scene;
    this.resize();
    const { canvas, ctx, low, camera } = this;
    const W = canvas.width,
      H = canvas.height;
    const b = camera.basis();
    const img = this.lowImage,
      data = img.data,
      w = low.width,
      h = low.height;
    const xyz = new Float64Array(3),
      tone = scene.tone;
    const v = [0, 0, 0];
    const sx = W / w,
      sy = H / h;
    for (let j = 0; j < h; j++) {
      for (let i = 0; i < w; i++) {
        camera.ray((i + 0.5) * sx, (j + 0.5) * sy, W, H, b, v);
        const el = Math.asin(v[2]) / rad;
        const az = ((Math.atan2(v[0], v[1]) / rad) + 360) % 360;
        if (el >= 0) skyXYZ(scene, el, az, xyz);
        else groundXYZ(scene, el, az, v, xyz);
        tone.write(xyz[0], xyz[1], xyz[2], data, (j * w + i) * 4);
        data[(j * w + i) * 4 + 3] = 255;
      }
    }
    this.lowCtx.putImageData(img, 0, 0);
    ctx.imageSmoothingEnabled = true;
    ctx.imageSmoothingQuality = "high";
    ctx.drawImage(low, 0, 0, W, H);
    this.drawStars(scene, b, W, H);
    this.clipAboveHorizon(b, W, H, () => {
      this.drawSun(scene, b, W, H);
      this.drawEarth(scene, b, W, H);
    });
    this.positions = {
      sun: camera.project(scene.sunDir, W, H, b),
      earth: camera.project(scene.earthDir, W, H, b),
    };
  }
  clipAboveHorizon(b, W, H, draw) {
    const { ctx } = this;
    ctx.save();
    ctx.beginPath();
    ctx.rect(0, 0, W, this.camera.horizonY(H, W));
    ctx.clip();
    draw();
    ctx.restore();
  }
  drawStars(scene, b, W, H) {
    const { ctx, camera } = this;
    const stars = this.assets.stars;
    if (!stars) return;
    const m = mul(scene.frame.eclToLocal, starPrecession(scene.frame.jdTT));
    const tone = scene.tone;
    const ppd = camera.pixelsPerDegree(W);
    const xyz = new Float64Array(3);
    ctx.save();
    ctx.globalCompositeOperation = "lighter";
    for (const star of stars) {
      const d = apply(m, star.v);
      if (d[2] <= 0.002) continue;
      const p = camera.project(d, W, H, b);
      if (!p || p[0] < -4 || p[1] < -4 || p[0] > W + 4 || p[1] > H + 4) continue;
      const el = Math.asin(d[2]) / rad;
      const T = Math.exp(-RAYLEIGH_TAU_550 * airmass(el));
      const Ls = (star.lux * T) / STAR_PSF_SR;
      skyXYZ(scene, el, ((Math.atan2(d[0], d[1]) / rad) + 360) % 360, xyz);
      const contrast = tone.value(xyz[1] + Ls) - tone.value(xyz[1]);
      if (contrast < 0.025) continue;
      const a = clamp(contrast * 1.4, 0, 1);
      const radius = clamp(0.6 + 1.6 * a, 0.6, 2.4) * this.dpr * (ppd > 60 ? 1.4 : 1);
      const c = star.colour;
      ctx.fillStyle = `rgba(${(255 * c[0]) | 0},${(255 * c[1]) | 0},${(255 * c[2]) | 0},${a})`;
      ctx.beginPath();
      ctx.arc(p[0], p[1], radius, 0, Math.PI * 2);
      ctx.fill();
    }
    ctx.restore();
  }
  drawSun(scene, b, W, H) {
    if (scene.sunEl + scene.state.sun_radius_deg <= 0) return;
    const p = this.camera.project(scene.sunDir, W, H, b);
    if (!p) return;
    const { ctx } = this;
    const ppd = this.camera.pixelsPerDegree(W);
    const r = Math.max(2 * this.dpr, scene.state.sun_radius_deg * ppd);
    const [x, y] = scene.sunBeamXY;
    const rgb = xyToRGB([x, y]).map((v) => Math.max(v, 0));
    const m = Math.max(...rgb);
    const c = rgb.map((v) => Math.round(255 * Math.pow(v / m, 0.45)));
    // Glare (display): a soft bloom scaled with the beam against the sky.
    const glare = clamp(Math.log10(1 + scene.sunNormalLux / Math.max(scene.state.total, 1)) / 1.2, 0.15, 1);
    const g = ctx.createRadialGradient(p[0], p[1], r * 0.5, p[0], p[1], r * 2 + ppd * 7 * glare);
    g.addColorStop(0, `rgba(${c[0]},${c[1]},${c[2]},${0.75 * glare})`);
    g.addColorStop(0.25, `rgba(${c[0]},${c[1]},${c[2]},${0.22 * glare})`);
    g.addColorStop(1, `rgba(${c[0]},${c[1]},${c[2]},0)`);
    ctx.save();
    ctx.globalCompositeOperation = "lighter";
    ctx.fillStyle = g;
    ctx.beginPath();
    ctx.arc(p[0], p[1], r * 2 + ppd * 7 * glare, 0, Math.PI * 2);
    ctx.fill();
    ctx.restore();
    ctx.fillStyle = `rgb(${255},${Math.min(255, c[1] + 40)},${Math.min(255, c[2] + 60)})`;
    ctx.beginPath();
    ctx.arc(p[0], p[1], r, 0, Math.PI * 2);
    ctx.fill();
  }
  drawEarth(scene, b, W, H) {
    if (scene.earthEl + scene.state.earth_radius_deg <= 0) return;
    const p = this.camera.project(scene.earthDir, W, H, b);
    if (!p || !this.assets.earthTexture) return;
    const ppd = this.camera.pixelsPerDegree(W);
    const D = scene.state.earth_radius_deg * 2 * ppd;
    if (p[0] < -D || p[1] < -D || p[0] > W + D || p[1] > H + D) return;
    const n = clamp(Math.round(D), 8, 640);
    const behind = skyXYZ(scene, Math.max(scene.earthEl, 0), scene.earthAz, new Float64Array(3));
    const img = earthImage(scene, this.assets.earthTexture, n, scene.tone, { behind, local: true });
    const ec = this.earthCanvas;
    ec.width = n;
    ec.height = n;
    const ectx = ec.getContext("2d");
    ectx.putImageData(new ImageData(img.data, n, n), 0, 0);
    // Rotate so Earth's local "up" (toward zenith) points along its screen direction.
    const q = this.camera.project(
      scene.earthDir.map((v, i) => v + img.up[i] * 0.01),
      W,
      H,
      b,
    );
    const angle = q ? Math.atan2(q[0] - p[0], -(q[1] - p[1])) : 0;
    const { ctx } = this;
    ctx.save();
    ctx.translate(p[0], p[1]);
    ctx.rotate(angle);
    ctx.imageSmoothingEnabled = true;
    // The lit face is added onto the sky as light.
    ctx.globalCompositeOperation = "lighter";
    ctx.globalAlpha = img.visibility;
    ctx.drawImage(ec, -D / 2, -D / 2, D, D);
    ctx.restore();
    this.earthScreen = { x: p[0], y: p[1], d: D, angle };
  }
}

// ---------- Fisheye dome for the twilight comparison ----------
export function renderDome(canvas, atlas, sunEl, tone) {
  const ctx = canvas.getContext("2d");
  const n = canvas.width;
  const img = ctx.createImageData(n, n);
  const slice = atlas.slice(sunEl);
  const xyz = new Float64Array(3);
  const R = n / 2;
  for (let j = 0; j < n; j++)
    for (let i = 0; i < n; i++) {
      const x = (i + 0.5 - R) / R,
        y = (j + 0.5 - R) / R,
        r = Math.hypot(x, y);
      const k = (j * n + i) * 4;
      if (r > 1) continue;
      const el = 90 * (1 - r);
      // Sun's azimuth at the bottom of the dome.
      const az = (Math.atan2(x, y) / rad + 360) % 360;
      xyz[0] = xyz[1] = xyz[2] = 0;
      atlas.sample(slice, el, az, 1, xyz);
      tone.write(xyz[0], xyz[1], xyz[2], img.data, k);
      img.data[k + 3] = Math.round(255 * clamp((1 - r) * R, 0, 1));
    }
  ctx.putImageData(img, 0, 0);
}
