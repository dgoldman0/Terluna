// An orthographic Moon globe for choosing places. The surface is the geography
// atlas (seas by depth, land by height with relief); the light on it is the
// calendar's Sun and Earth directions with the solved clear-sky light at each
// point's Sun elevation.
import { interpolate } from "./model/evaluator.mjs";

const rad = Math.PI / 180;
const clamp = (x, a, b) => Math.max(a, Math.min(b, x));

export class Globe {
  constructor(canvas, texture, options = {}) {
    this.canvas = canvas;
    this.ctx = canvas.getContext("2d");
    this.texture = texture;
    this.lon0 = options.lon0 ?? -10;
    this.lat0 = options.lat0 ?? 10;
    this.places = [];
    this.selected = null;
    this.hover = null;
    this.onPick = options.onPick ?? (() => {});
    this.onHover = options.onHover ?? (() => {});
    this.buffer = document.createElement("canvas");
    this.bctx = this.buffer.getContext("2d");
    this.bind();
  }
  setLight(sunSel, earthSel, solarLux, earthLux) {
    this.sun = sunSel;
    this.earth = earthSel;
    this.solarLux = solarLux;
    this.earthLux = earthLux;
  }
  resize() {
    const r = this.canvas.getBoundingClientRect();
    const dpr = Math.min(window.devicePixelRatio || 1, 2);
    const n = Math.max(64, Math.round(r.width * dpr));
    if (this.canvas.width !== n) {
      this.canvas.width = n;
      this.canvas.height = n;
    }
    this.dpr = dpr;
    // The surface is drawn at a modest resolution and scaled; overlays are sharp.
    const m = Math.min(n, 560);
    if (this.buffer.width !== m) {
      this.buffer.width = m;
      this.buffer.height = m;
      this.image = this.bctx.createImageData(m, m);
    }
  }
  basis() {
    const l = this.lon0 * rad,
      p = this.lat0 * rad;
    const c = [Math.cos(p) * Math.cos(l), Math.cos(p) * Math.sin(l), Math.sin(p)];
    const e = [-Math.sin(l), Math.cos(l), 0];
    const n = [-Math.sin(p) * Math.cos(l), -Math.sin(p) * Math.sin(l), Math.cos(p)];
    return { c, e, n };
  }
  project(lon, lat, size) {
    const b = this.basis(),
      v = [
        Math.cos(lat * rad) * Math.cos(lon * rad),
        Math.cos(lat * rad) * Math.sin(lon * rad),
        Math.sin(lat * rad),
      ];
    const z = v[0] * b.c[0] + v[1] * b.c[1] + v[2] * b.c[2];
    const x = v[0] * b.e[0] + v[1] * b.e[1] + v[2] * b.e[2];
    const y = v[0] * b.n[0] + v[1] * b.n[1] + v[2] * b.n[2];
    const R = size * 0.47;
    return { x: size / 2 + x * R, y: size / 2 - y * R, front: z > 0 };
  }
  unproject(px, py, size) {
    const R = size * 0.47;
    const x = (px - size / 2) / R,
      y = (size / 2 - py) / R,
      r2 = x * x + y * y;
    if (r2 > 1) return null;
    const z = Math.sqrt(1 - r2),
      b = this.basis();
    const v = [0, 1, 2].map((i) => x * b.e[i] + y * b.n[i] + z * b.c[i]);
    return [Math.atan2(v[1], v[0]) / rad, Math.asin(clamp(v[2], -1, 1)) / rad];
  }
  render(transfer) {
    this.resize();
    const m = this.buffer.width,
      data = this.image.data,
      tex = this.texture,
      tw = tex.width,
      th = tex.height,
      td = tex.data;
    const b = this.basis(),
      R = m * 0.47,
      sun = this.sun,
      earth = this.earth;
    // Light levels by Sun elevation from the calendar's clear-sky table.
    const el = transfer.direct_elevation_deg,
      fe = transfer.diffuse_elevation_deg;
    const lightAt = (e) =>
      interpolate(el, transfer.solar.direct_lux, e) + interpolate(fe, transfer.solar.diffuse_lux, e);
    const lut = new Float32Array(721);
    for (let i = 0; i <= 720; i++) lut[i] = lightAt(i / 4 - 90);
    this.duskDepression = 0;
    for (let i = 360; i >= 0; i--)
      if (lut[i] < 2.98) {
        this.duskDepression = 90 - i / 4;
        break;
      }
    const lnMax = Math.log10(90000),
      lnMin = Math.log10(0.003);
    for (let j = 0; j < m; j++)
      for (let i = 0; i < m; i++) {
        const k = (j * m + i) * 4;
        const x = (i + 0.5 - m / 2) / R,
          y = (m / 2 - j - 0.5) / R,
          r2 = x * x + y * y;
        if (r2 > 1) {
          data[k + 3] = 0;
          continue;
        }
        const z = Math.sqrt(1 - r2);
        const v0 = x * b.e[0] + y * b.n[0] + z * b.c[0],
          v1 = x * b.e[1] + y * b.n[1] + z * b.c[1],
          v2 = x * b.e[2] + y * b.n[2] + z * b.c[2];
        const lon = Math.atan2(v1, v0) / rad,
          lat = Math.asin(clamp(v2, -1, 1)) / rad;
        const tx = Math.min(tw - 1, Math.floor(((lon + 180) / 360) * tw)),
          ty = Math.min(th - 1, Math.floor(((90 - lat) / 180) * th));
        const t = (ty * tw + tx) * 4;
        let r = td[t],
          g = td[t + 1],
          bl = td[t + 2];
        if (sun) {
          const sunEl = Math.asin(clamp(v0 * sun[0] + v1 * sun[1] + v2 * sun[2], -1, 1)) / rad;
          const f = (sunEl + 90) * 4,
            i0 = Math.min(719, Math.max(0, Math.floor(f)));
          const solar = lut[i0] + (lut[i0 + 1] - lut[i0]) * (f - i0);
          let earthLight = 0;
          if (earth && this.earthLux > 0) {
            const ee = v0 * earth[0] + v1 * earth[1] + v2 * earth[2];
            if (ee > -0.02) earthLight = this.earthLux * clamp(ee * 1.2 + 0.1, 0, 1);
          }
          const total = solar + earthLight;
          const level = clamp((Math.log10(Math.max(total, 1e-6)) - lnMin) / (lnMax - lnMin), 0, 1);
          const bright = 0.16 + 0.84 * level ** 1.6;
          // Tints: warm in the twilight band, blue where Earthlight leads.
          let tr = 1,
            tg = 1,
            tb = 1;
          if (sunEl < 0 && solar >= earthLight) {
            const w = clamp(1 - Math.abs(sunEl + 14) / 30, 0, 1) * 0.55;
            tr = 1 + 0.25 * w;
            tg = 1 - 0.05 * w;
            tb = 1 - 0.35 * w;
          } else if (earthLight > solar) {
            tr = 0.72;
            tg = 0.92;
            tb = 1.25;
          } else if (sunEl < 0) {
            tr = 0.7;
            tg = 0.78;
            tb = 1.15;
          }
          r = r * bright * tr;
          g = g * bright * tg;
          bl = bl * bright * tb;
        }
        // A soft limb.
        const limb = 0.72 + 0.28 * Math.sqrt(z);
        data[k] = r * limb;
        data[k + 1] = g * limb;
        data[k + 2] = bl * limb;
        data[k + 3] = Math.round(255 * clamp((1 - Math.sqrt(r2)) * R * 1.5, 0, 1));
      }
    this.bctx.putImageData(this.image, 0, 0);
    const { ctx, canvas } = this;
    const n = canvas.width;
    ctx.clearRect(0, 0, n, n);
    // Halo of the air.
    const halo = ctx.createRadialGradient(n / 2, n / 2, n * 0.45, n / 2, n / 2, n * 0.5);
    halo.addColorStop(0, "rgba(140,190,255,0.35)");
    halo.addColorStop(1, "rgba(140,190,255,0)");
    ctx.fillStyle = halo;
    ctx.fillRect(0, 0, n, n);
    ctx.imageSmoothingEnabled = true;
    ctx.drawImage(this.buffer, 0, 0, n, n);
    this.drawOverlay();
  }
  drawOverlay() {
    const { ctx, canvas } = this;
    const n = canvas.width,
      d = this.dpr;
    // Graticule.
    ctx.save();
    ctx.strokeStyle = "rgba(255,255,255,0.13)";
    ctx.lineWidth = 1 * d;
    for (let lat = -60; lat <= 60; lat += 30) this.line(Array.from({ length: 73 }, (_, i) => [i * 5 - 180, lat]), n);
    for (let lon = -180; lon < 180; lon += 30) this.line(Array.from({ length: 37 }, (_, i) => [lon, i * 5 - 90]), n);
    ctx.restore();
    // The terminator, and the line where sunlight alone falls to the end of
    // practical twilight (2.98 lux).
    if (this.sun) {
      const circle = (radiusDeg) => {
        const s = this.sun;
        const a0 = Math.abs(s[2]) < 0.9 ? [0, 0, 1] : [1, 0, 0];
        let a = [s[1] * a0[2] - s[2] * a0[1], s[2] * a0[0] - s[0] * a0[2], s[0] * a0[1] - s[1] * a0[0]];
        const al = Math.hypot(...a);
        a = a.map((v) => v / al);
        const b = [s[1] * a[2] - s[2] * a[1], s[2] * a[0] - s[0] * a[2], s[0] * a[1] - s[1] * a[0]];
        const r = radiusDeg * rad;
        return Array.from({ length: 145 }, (_, k) => {
          const t = (k / 144) * 2 * Math.PI;
          const v = [0, 1, 2].map(
            (i) => Math.cos(r) * s[i] + Math.sin(r) * (Math.cos(t) * a[i] + Math.sin(t) * b[i]),
          );
          return [Math.atan2(v[1], v[0]) / rad, Math.asin(clamp(v[2], -1, 1)) / rad];
        });
      };
      ctx.save();
      ctx.lineWidth = 2 * d;
      ctx.strokeStyle = "rgba(255,214,120,0.95)";
      this.line(circle(90), n);
      if (this.duskDepression) {
        ctx.setLineDash([6 * d, 5 * d]);
        ctx.strokeStyle = "rgba(196,170,255,0.9)";
        this.line(circle(90 + this.duskDepression), n);
      }
      ctx.restore();
    }
    // Places.
    for (const place of this.places) {
      const p = this.project(place.longitude, place.latitude, n);
      if (!p.front) continue;
      const active = place === this.hover;
      ctx.fillStyle = place.colour ?? "rgba(255,255,255,0.85)";
      ctx.strokeStyle = "rgba(10,14,30,0.8)";
      ctx.lineWidth = 1.5 * d;
      ctx.beginPath();
      ctx.arc(p.x, p.y, (active ? 5 : 3.2) * d, 0, Math.PI * 2);
      ctx.fill();
      ctx.stroke();
    }
    // Under the Sun and under the Earth.
    const mark = (v, colour, label) => {
      if (!v) return;
      const lon = Math.atan2(v[1], v[0]) / rad,
        lat = Math.asin(v[2]) / rad;
      const p = this.project(lon, lat, n);
      if (!p.front) return;
      ctx.save();
      ctx.fillStyle = colour;
      ctx.shadowColor = colour;
      ctx.shadowBlur = 12 * d;
      ctx.beginPath();
      ctx.arc(p.x, p.y, 6 * d, 0, Math.PI * 2);
      ctx.fill();
      ctx.restore();
      ctx.font = `${600} ${11 * d}px system-ui, sans-serif`;
      ctx.fillStyle = "rgba(255,255,255,0.92)";
      ctx.strokeStyle = "rgba(5,8,20,0.75)";
      ctx.lineWidth = 3 * d;
      ctx.strokeText(label, p.x + 10 * d, p.y + 4 * d);
      ctx.fillText(label, p.x + 10 * d, p.y + 4 * d);
    };
    mark(this.sun, "#ffc53d", "Sun overhead");
    mark(this.earth, "#38c6f4", "Earth overhead");
    if (this.selected) {
      const p = this.project(this.selected.longitude, this.selected.latitude, n);
      if (p.front) {
        ctx.save();
        ctx.strokeStyle = "#fff";
        ctx.lineWidth = 2.5 * d;
        ctx.shadowColor = "rgba(0,0,0,0.6)";
        ctx.shadowBlur = 6 * d;
        ctx.beginPath();
        ctx.arc(p.x, p.y, 10 * d, 0, Math.PI * 2);
        ctx.stroke();
        ctx.fillStyle = "#ff5d8f";
        ctx.beginPath();
        ctx.arc(p.x, p.y, 4.5 * d, 0, Math.PI * 2);
        ctx.fill();
        ctx.restore();
      }
    }
  }
  line(points, n) {
    const { ctx } = this;
    ctx.beginPath();
    let drawing = false;
    for (const [lon, lat] of points) {
      const p = this.project(lon, lat, n);
      if (!p.front) {
        drawing = false;
        continue;
      }
      if (drawing) ctx.lineTo(p.x, p.y);
      else ctx.moveTo(p.x, p.y);
      drawing = true;
    }
    ctx.stroke();
  }
  placeAt(px, py) {
    const n = this.canvas.width;
    let best = null,
      bestD = 12 * this.dpr;
    for (const place of this.places) {
      const p = this.project(place.longitude, place.latitude, n);
      if (!p.front) continue;
      const dd = Math.hypot(p.x - px, p.y - py);
      if (dd < bestD) {
        best = place;
        bestD = dd;
      }
    }
    return best;
  }
  spinTo(lon, lat, done) {
    const from = [this.lon0, this.lat0];
    let dl = ((lon - from[0] + 540) % 360) - 180;
    const to = [from[0] + dl, clamp(lat, -70, 70)];
    const start = performance.now();
    const step = (now) => {
      const t = clamp((now - start) / 700, 0, 1),
        e = t < 0.5 ? 2 * t * t : 1 - (-2 * t + 2) ** 2 / 2;
      this.lon0 = from[0] + (to[0] - from[0]) * e;
      this.lat0 = from[1] + (to[1] - from[1]) * e;
      done?.(t >= 1);
      if (t < 1) requestAnimationFrame(step);
    };
    requestAnimationFrame(step);
  }
  bind() {
    const c = this.canvas;
    let drag = null;
    const pos = (e) => {
      const r = c.getBoundingClientRect();
      return [(e.clientX - r.left) * this.dpr, (e.clientY - r.top) * this.dpr];
    };
    c.addEventListener("pointerdown", (e) => {
      c.setPointerCapture(e.pointerId);
      drag = { x: e.clientX, y: e.clientY, lon: this.lon0, lat: this.lat0, moved: false };
    });
    c.addEventListener("pointermove", (e) => {
      if (drag) {
        const dx = e.clientX - drag.x,
          dy = e.clientY - drag.y;
        if (Math.hypot(dx, dy) > 4) drag.moved = true;
        const k = 180 / (c.getBoundingClientRect().width * 0.94);
        this.lon0 = drag.lon - dx * k;
        this.lat0 = clamp(drag.lat + dy * k, -85, 85);
        this.onDrag?.();
        return;
      }
      const [x, y] = pos(e);
      const place = this.placeAt(x, y);
      if (place !== this.hover) {
        this.hover = place;
        this.onHover(place, e);
        this.drawOverlayOnly?.();
      }
    });
    const end = (e) => {
      if (!drag) return;
      const moved = drag.moved;
      drag = null;
      if (moved) return;
      const [x, y] = pos(e);
      const place = this.placeAt(x, y);
      if (place) this.onPick(place);
      else {
        const ll = this.unproject(x, y, c.width);
        if (ll) this.onPick({ custom: true, longitude: ll[0], latitude: ll[1] });
      }
    };
    c.addEventListener("pointerup", end);
    c.addEventListener("pointercancel", () => (drag = null));
    c.addEventListener("pointerleave", () => {
      if (this.hover) {
        this.hover = null;
        this.onHover(null);
      }
    });
    c.addEventListener("keydown", (e) => {
      const keys = { ArrowLeft: [-10, 0], ArrowRight: [10, 0], ArrowUp: [0, 10], ArrowDown: [0, -10] };
      if (keys[e.key]) {
        e.preventDefault();
        this.lon0 += keys[e.key][0];
        this.lat0 = clamp(this.lat0 + keys[e.key][1], -85, 85);
        this.onDrag?.();
      }
    });
  }
}
