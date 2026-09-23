/* Conservative illustrative surface columns, millimetres over each grid cell.
 * Canopy interception, film, soil and ponded stores have separate ledgers.
 * Static D8 routing on a 4 m grid resolves downhill redistribution and boundary
 * export. Within-cell pond morphology and a connected free surface are unresolved.
 */
import { OM } from './om.js';
import L from '../world/landscape.js';
import C from './core.js';
class SurfaceWater {
  constructor(options = {}) {
    this.n = options.n ?? 65;
    this.cell = options.cell ?? 4;
    this.minX = options.minX ?? -128;
    this.minZ = options.minZ ?? -64;
    if (
      !Number.isInteger(this.n) ||
      this.n < 3 ||
      this.n > 257 ||
      !Number.isFinite(this.cell) ||
      this.cell <= 0 ||
      !Number.isFinite(this.minX) ||
      !Number.isFinite(this.minZ)
    )
      throw new RangeError('Invalid surface-water grid');
    const N = this.n * this.n;
    this.meta = { n: this.n, cell: this.cell, minX: this.minX, minZ: this.minZ };
    this.elevation = new Float64Array(N);
    this.cover = new Float32Array(N);
    this.capacity = new Float32Array(N);
    this.permeability = new Float32Array(N);
    this.roof = new Uint8Array(N);
    this.slope = new Float32Array(N);
    this.initialSoil = new Float64Array(N);
    for (let j = 0; j < this.n; j++)
      for (let i = 0; i < this.n; i++) {
        const k = j * this.n + i,
          x = this.minX + i * this.cell,
          z = this.minZ + j * this.cell,
          f = L.sample(x, z);
        this.elevation[k] = f.elevation;
        this.cover[k] = f.canopy;
        this.capacity[k] = 12 + 75 * f.soilDepth;
        this.permeability[k] = 0.0006 + 0.006 * f.weights[0] + 0.002 * f.weights[3];
        this.roof[k] = C.roofMask(x, z);
        this.slope[k] = f.slope;
        this.initialSoil[k] = f.elevation > 0 ? this.capacity[k] * f.moisture * 0.65 : 0;
      }
    this.flow = L.drainageGrid(this.elevation, this.n, this.cell);
    this.reset();
  }
  reset() {
    const N = this.n * this.n;
    for (const k of ['leaf', 'film', 'soil', 'pond']) this[k] = new Float64Array(N);
    this.soil.set(this.initialSoil);
    this.initialStorage = this.initialSoil.reduce((a, b) => a + b, 0);
    this.time = 0;
    this.elapsed = 0;
    this.kind = 'clear';
    this.totalInput = 0;
    this.evaporation = 0;
    this.deepDrainage = 0;
    this.boundaryExport = 0;
    this.roofExport = 0;
    this.steps = 0;
    this.version = (this.version || 0) + 1;
    this.buffer = new Float32Array(N * 4);
  }
  advance(dt, w) {
    if (!Number.isFinite(dt) || dt < 0)
      throw new RangeError('Water step must be finite and nonnegative');
    if (dt > 10) {
      for (let t = 0; t < dt; t += 10) this.advance(Math.min(10, dt - t), w);
      return;
    }
    if (dt === 0) return;
    const N = this.n * this.n,
      rain = (w.rain / 3600) * dt,
      evap =
        0.000065 *
        (1 - w.humidity) *
        (1 + w.wind * 0.2) *
        Math.exp((w.temperature - 290) / 25) *
        dt;
    this.totalInput += rain * N;
    for (let k = 0; k < N; k++) {
      const inputRain = this.roof[k] ? 0 : rain;
      if (this.roof[k]) this.roofExport += rain;
      if (this.elevation[k] <= 0) {
        this.boundaryExport += inputRain;
        continue;
      }
      const intercept = Math.min(
        inputRain * this.cover[k],
        Math.max(0, 0.35 * this.cover[k] - this.leaf[k]),
      );
      this.leaf[k] += intercept;
      this.film[k] += inputRain - intercept;
      const le = Math.min(this.leaf[k], evap * 1.5),
        drip = Math.min(this.leaf[k] - le, (this.leaf[k] - le) * 0.0005 * dt);
      this.leaf[k] -= le + drip;
      this.film[k] += drip;
      this.evaporation += le;
      const infiltrate = Math.min(
        this.film[k],
        this.permeability[k] * dt * Math.max(0.05, 1 - this.soil[k] / this.capacity[k]),
        this.capacity[k] - this.soil[k],
      );
      this.film[k] -= infiltrate;
      this.soil[k] += infiltrate;
      const pondInfiltration = Math.min(
        this.pond[k],
        Math.max(
          0,
          this.permeability[k] * dt * Math.max(0.05, 1 - this.soil[k] / this.capacity[k]) -
            infiltrate,
        ),
        this.capacity[k] - this.soil[k],
      );
      this.pond[k] -= pondInfiltration;
      this.soil[k] += pondInfiltration;
      const fe = Math.min(this.film[k], evap * (1 - 0.7 * this.cover[k]));
      this.film[k] -= fe;
      this.evaporation += fe;
      const se = Math.min(this.soil[k], evap * 0.25 * (1 - 0.6 * this.cover[k]));
      this.soil[k] -= se;
      this.evaporation += se;
      const deep = this.soil[k] * (1 - Math.exp(-0.000005 * dt));
      this.soil[k] -= deep;
      this.deepDrainage += deep;
      const excess = Math.max(0, this.film[k] - 0.24);
      this.film[k] -= excess;
      this.pond[k] += excess;
    }
    // Process downstream first so transferred water advances at most one cell per
    // substep. This avoids a timestep-dependent instantaneous catchment sweep.
    const order = this.flow.order;
    for (let q = order.length - 1; q >= 0; q--) {
      const k = order[q],
        receiver = this.flow.receiver[k];
      if (this.pond[k] <= 0) continue;
      const pe = Math.min(this.pond[k], evap);
      this.pond[k] -= pe;
      this.evaporation += pe;
      const depression = Math.max(0, this.flow.filled[k] - this.elevation[k]) * 1000;
      const mobile = Math.max(0, this.pond[k] - Math.min(depression, 50));
      const out = mobile * (1 - Math.exp(-(0.005 + 0.035 * Math.sqrt(this.slope[k])) * dt));
      this.pond[k] -= out;
      if (receiver < 0 || this.elevation[receiver] <= 0) this.boundaryExport += out;
      else this.pond[receiver] += out;
    }
    this.elapsed += dt;
    this.steps++;
    this.version++;
  }
  seek(time, kind = 'episode') {
    if (time < 0 || !Number.isFinite(time))
      throw new RangeError('Water time must be finite and nonnegative');
    if (time < this.time || kind !== this.kind) {
      this.reset();
      this.kind = kind;
    }
    for (let t = this.time; t < time - 1e-9; t += 10) {
      const dt = Math.min(10, time - t);
      this.advance(dt, C.weatherAt(t + dt * 0.5, kind));
    }
    this.time = time;
    return this;
  }
  get ledger() {
    let storage = 0;
    for (const a of [this.leaf, this.film, this.soil, this.pond]) for (const v of a) storage += v;
    const exported = this.boundaryExport + this.roofExport + this.deepDrainage,
      residual = this.initialStorage + this.totalInput - this.evaporation - exported - storage;
    return {
      units: 'mm-cell; multiply by cell area / 1000 for m3',
      cellArea_m2: this.cell * this.cell,
      initialStorage: this.initialStorage,
      input: this.totalInput,
      storage,
      evaporated: this.evaporation,
      deepDrainage: this.deepDrainage,
      boundaryExport: this.boundaryExport,
      roofExport: this.roofExport,
      residual,
      relativeResidual: Math.abs(residual) / Math.max(1, this.initialStorage + this.totalInput),
      steps: this.steps,
      elapsedSeconds: this.elapsed,
      forcingTimeSeconds: this.time,
    };
  }
  textureData() {
    for (let k = 0; k < this.n * this.n; k++) {
      this.buffer[k * 4] = C.clamp(this.film[k] / 0.24);
      this.buffer[k * 4 + 1] = C.clamp(this.soil[k] / this.capacity[k]);
      this.buffer[k * 4 + 2] = this.pond[k];
      this.buffer[k * 4 + 3] = C.clamp(this.leaf[k] / Math.max(0.01, 0.35 * this.cover[k]));
    }
    return this.buffer;
  }
  sample(x, z) {
    const channel = (i) => {
      const u = C.clamp((x - this.minX) / this.cell, 0, this.n - 1),
        v = C.clamp((z - this.minZ) / this.cell, 0, this.n - 1),
        a = Math.min(this.n - 2, Math.floor(u)),
        b = Math.min(this.n - 2, Math.floor(v)),
        tx = u - a,
        tz = v - b,
        at = (xx, zz) => this.buffer[(zz * this.n + xx) * 4 + i];
      return C.mix(
        C.mix(at(a, b), at(a + 1, b), tx),
        C.mix(at(a, b + 1), at(a + 1, b + 1), tx),
        tz,
      );
    };
    this.textureData();
    return {
      film: channel(0),
      soilSaturation: channel(1),
      ponded_mm: channel(2),
      leafWetness: channel(3),
    };
  }
}
OM.SurfaceWater = SurfaceWater;
export default SurfaceWater;
