/* Altitude / local-solar-angle transport for the inherited optical atmosphere.
 * Molecular profile is held fixed for the reference comparison (H=48.4/8 km).
 * RGB-band proxy; this does not replace the spectral clear-sky atlas. The optional
 * ozone layer is explicitly selected, with 300 DU and effective RGB absorption.
 */
import { EARTH_RADIUS, MOON_RADIUS } from '../../shared/constants.js';
const BETA = Object.freeze([5.8e-6, 13.5e-6, 33.1e-6]);
function profile(world) {
  if (!['moon', 'earth', 'moon_no_ozone'].includes(world))
    throw new RangeError('Unknown optical world');
  return {
    R: world === 'earth' ? EARTH_RADIUS : MOON_RADIUS,
    H: world === 'earth' ? 8000 : 48400,
    density: world === 'earth' ? 1 : 1.2,
    top: world === 'earth' ? 120000 : 600000,
    ozone: world === 'moon_no_ozone' ? 0 : 1,
  };
}
function altitude(R, h, mu, s) {
  return Math.sqrt((R + h) ** 2 + s * s + 2 * (R + h) * mu * s) - R;
}
function exitDistance(R, h, mu, top) {
  const r = R + h;
  return -r * mu + Math.sqrt(r * r * mu * mu + (R + top) ** 2 - r * r);
}
function blocked(R, h, mu) {
  return mu < 0 && (R + h) * Math.sqrt(Math.max(0, 1 - mu * mu)) < R;
}
// Unit-integral, selected 10–40 km triangular ozone number profile.
function ozoneDensity(h) {
  return h < 10000 || h > 40000
    ? 0
    : (h < 25000 ? (h - 10000) / 15000 : (40000 - h) / 15000) / 15000;
}
function columns(world, h, mu, steps = 96) {
  const p = profile(world);
  if (
    !Number.isFinite(h) ||
    h < 0 ||
    h >= p.top ||
    !Number.isFinite(mu) ||
    mu < -1 ||
    mu > 1 ||
    !Number.isInteger(steps) ||
    steps < 8
  )
    throw new RangeError('Invalid optical ray');
  if (blocked(p.R, h, mu)) return { molecular: Infinity, ozone: Infinity, blocked: true };
  const end = exitDistance(p.R, h, mu, p.top);
  let molecular = 0,
    ozone = 0;
  // Nonuniform quadrature retains near-ground resolution on upward long rays.
  for (let i = 0; i < steps; i++) {
    const a = end * (i / steps) ** 2,
      b = end * ((i + 1) / steps) ** 2,
      z = altitude(p.R, h, mu, (a + b) / 2);
    molecular += Math.exp(-z / p.H) * (b - a) * p.density;
    ozone += ozoneDensity(z) * (b - a);
  }
  return { molecular, ozone: ozone * p.ozone, blocked: false };
}
const OZONE_TAU = Object.freeze([0.015, 0.035, 0.003]);
function transmission(world, h, mu, steps = 96) {
  const c = columns(world, h, mu, steps);
  if (c.blocked) return [0, 0, 0];
  return BETA.map((b, i) => Math.exp(-b * c.molecular - OZONE_TAU[i] * c.ozone));
}
function segmentColumn(world, h, mu, length, steps = 16) {
  const p = profile(world);
  let value = 0;
  for (let i = 0; i < steps; i++)
    value +=
      (p.density * Math.exp(-altitude(p.R, h, mu, (length * (i + 0.5)) / steps) / p.H) * length) /
      steps;
  return value;
}
function buildTable(world, maxHeight, width = 128, height = 64) {
  const data = new Float32Array(width * height * 4),
    max = Math.max(1000, maxHeight);
  for (let j = 0; j < height; j++)
    for (let i = 0; i < width; i++) {
      const z = (max * j) / (height - 1),
        v = (2 * i) / (width - 1) - 1,
        mu = Math.sign(v) * v * v,
        t = transmission(world, z, mu);
      data.set([...t, 1], (j * width + i) * 4);
    }
  return { data, width, height, maxHeight: max };
}
const API = {
  BETA,
  OZONE_TAU,
  profile,
  altitude,
  exitDistance,
  blocked,
  ozoneDensity,
  columns,
  transmission,
  segmentColumn,
  buildTable,
};
export default API;
