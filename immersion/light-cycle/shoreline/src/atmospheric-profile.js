/* A1 shared thermodynamic/optical state. SI throughout.
 * Existing sounding nodes are preserved. The upper continuation is an explicit
 * optical boundary assumption, never a prediction of the upper lunar atmosphere.
 * A1_METHODS.md specifies the dry-air-equivalent Rayleigh closure and exclusions.
 */
(function (root) {
'use strict';
const W = root.OpenMoonWeatherColumn || (typeof require === 'function' ? require('./weather-column.js') : null);
const KB = 1.380649e-23;
const NREF = 101325 / (KB * 288.15);
const DU = 2.687e20;
function positive(x, label) { if (!Number.isFinite(x) || x <= 0) throw new RangeError(label + ' must be positive and finite'); return x; }
function canonical(value) {
  if (Array.isArray(value)) return '[' + value.map(canonical).join(',') + ']';
  if (value && typeof value === 'object') return '{' + Object.keys(value).sort().map(k => JSON.stringify(k) + ':' + canonical(value[k])).join(',') + '}';
  if (typeof value === 'number' && !Number.isFinite(value)) throw new TypeError('Nonfinite profile data');
  return JSON.stringify(value);
}
function freeze(x) { if (x && typeof x === 'object' && !Object.isFrozen(x)) { Object.values(x).forEach(freeze); Object.freeze(x); } return x; }
function derived(z, p, T, r, planet) {
  const e = p * r / (W.EPS + r), Tv = W.virtualT(T, r);
  return {z, p, T, r, e, Tv, rho: p / (W.K.Rd * Tv), g: planet.g0 * (planet.radius / (planet.radius + z)) ** 2,
    numberDry_m3: (p - e) / (KB * T), numberWater_m3: e / (KB * T), numberTotal_m3: p / (KB * T)};
}
function ozoneShape(z, spec) {
  if (z <= spec.bottom_m || z >= spec.top_m) return 0;
  const f = z < spec.peak_m ? (z - spec.bottom_m) / (spec.peak_m - spec.bottom_m) : (spec.top_m - z) / (spec.top_m - spec.peak_m);
  return 2 * f / (spec.top_m - spec.bottom_m);
}
function fromColumn(column, options = {}) {
  if (!column || column.schema !== 'open-moon-column/1') throw new TypeError('Expected a weather-column/1 state');
  const planet = {...column.planet}, last = column.rows.at(-1);
  const top = options.top_m ?? (column.world === 'earth' ? 120000 : 600000);
  const step = options.upperStep_m ?? (column.world === 'earth' ? 200 : 1000);
  const targetT = options.upperTemperature_K ?? last.T;
  const transition = options.temperatureTransition_m ?? (column.world === 'earth' ? 15000 : 80000);
  const waterScale = options.waterDecay_m ?? (column.world === 'earth' ? 8000 : 45000);
  for (const [label, value] of Object.entries({top, step, targetT, transition, waterScale})) positive(value, label);
  if (top <= last.z || targetT < 150 || targetT > 2000) throw new RangeError('Upper boundary or temperature outside declared domain');
  const ozone = {column_DU: column.world === 'moon_no_ozone' ? 0 : 300, bottom_m: 10000, peak_m: 25000, top_m: 40000, ...(options.ozone || {})};
  if (!Number.isFinite(ozone.column_DU) || ozone.column_DU < 0 || !(0 <= ozone.bottom_m && ozone.bottom_m < ozone.peak_m && ozone.peak_m < ozone.top_m && ozone.top_m < top)) throw new RangeError('Invalid selected ozone profile');
  const rows = column.rows.map(v => derived(v.z, v.p, v.T, v.r, planet));
  function upper(z) { const dz = z - last.z; return {T: last.T + (targetT - last.T) * -Math.expm1(-dz / transition), r: last.r * Math.exp(-dz / waterScale)}; }
  let z = last.z, logP = Math.log(last.p);
  // Two-point Gauss integration of hydrostatic d(log p)/dz, with spherical g.
  while (z < top) {
    const next = Math.min(top, z + step), half = (next - z) / 2, mid = (next + z) / 2;
    let integral = 0;
    for (const sign of [-1, 1]) { const h = mid + sign * half / Math.sqrt(3), s = upper(h); integral += planet.g0 * (planet.radius / (planet.radius + h)) ** 2 / (W.K.Rd * W.virtualT(s.T, s.r)); }
    logP -= half * integral; const s = upper(next); rows.push(derived(next, Math.exp(logP), s.T, s.r, planet)); z = next;
  }
  const state = {schema: 'open-moon-atmospheric-profile/1', world: column.world, regime: column.key, planet,
    sounding: {inputs: JSON.parse(JSON.stringify(column.inputs)), end_m: last.z, endPressure_Pa: last.p},
    upper: {kind: 'hydrostatic-temperature-relaxation', top_m: top, step_m: step, targetTemperature_K: targetT, transition_m: transition, waterDecay_m: waterScale},
    optics: {rayleigh: 'total molecular number density; dry-air-equivalent cross section', referenceNumberDensity_m3: NREF,
      referenceTemperature_K: 288.15, referencePressure_Pa: 101325, referenceRayleighAt1um_m1: 1.24062e-6,
      ozone, ozoneTemperature_K: 233, aerosolExtinction: 'zero in reference', waterVapourAbsorption: 'omitted', scatteringOrders: 'solver-specific'}, rows};
  freeze(state);
  return attach(state);
}
function attach(state) {
  const rows = state.rows;
  function sample(z) {
    if (!Number.isFinite(z) || z < 0 || z > state.upper.top_m) throw new RangeError('Altitude outside profile domain');
    let lo = 0, hi = rows.length - 1;
    while (hi - lo > 1) { const mid = (lo + hi) >> 1; if (rows[mid].z <= z) lo = mid; else hi = mid; }
    const a = rows[lo], b = rows[hi], t = (z - a.z) / (b.z - a.z);
    const r = derived(z, Math.exp(Math.log(a.p) + t * Math.log(b.p / a.p)), a.T + t * (b.T - a.T), a.r + t * (b.r - a.r), state.planet);
    r.airRelative = r.numberTotal_m3 / NREF;
    r.ozone_m3 = ozoneShape(z, state.optics.ozone) * state.optics.ozone.column_DU * DU;
    return r;
  }
  return Object.freeze({state, sample, serialize: () => canonical(state), fingerprint: async () => {
    const text = canonical(state);
    if (typeof require === 'function') return require('node:crypto').createHash('sha256').update(text).digest('hex');
    const bytes = await crypto.subtle.digest('SHA-256', new TextEncoder().encode(text));
    return Array.from(new Uint8Array(bytes), b => b.toString(16).padStart(2, '0')).join('');
  }});
}
function hydrate(state) {
  if (state?.schema !== 'open-moon-atmospheric-profile/1' || !Array.isArray(state.rows) || state.rows.length < 2) throw new TypeError('Invalid profile');
  let z = -1, p = Infinity;
  for (const row of state.rows) { if (!(row.z > z && row.p < p && row.T > 0 && row.r >= 0) || !Object.values(row).every(Number.isFinite)) throw new RangeError('Invalid profile rows'); z = row.z; p = row.p; }
  if (state.rows[0].z !== 0 || state.rows.at(-1).z !== state.upper.top_m) throw new RangeError('Profile boundary mismatch');
  return attach(freeze(JSON.parse(JSON.stringify(state))));
}
const API = {KB, NREF, DU, canonical, fromColumn, hydrate, ozoneShape};
root.OpenMoonAtmosphericProfile = API; if (typeof module !== 'undefined') module.exports = API;
})(globalThis);
