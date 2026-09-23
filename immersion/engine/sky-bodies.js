/* Sun, Earth and star geometry for the experience, from the baked site-sky product.
 *
 * The illumination domain owns the model (illumination/ephemeris.py). Its product
 * carries the site, the constants and golden samples; this module evaluates the same
 * closed-form mean-orbit formulas at any time, so playback can run past one month
 * without a jump, and tests/unit/sky-bodies.test.js checks it against the samples.
 *
 * Frames: the product uses local east-north-up; the scene uses x = west, y = up,
 * z = north, so a direction (e, n, u) becomes (-e, u, n).
 */
const TAU = 2 * Math.PI;
const rad = (deg) => (deg * Math.PI) / 180;

function unit(lat, lon) {
  return [Math.cos(lat) * Math.cos(lon), Math.cos(lat) * Math.sin(lon), Math.sin(lat)];
}
const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
const apply = (m, v) => m.map((row) => dot(row, v));
const multiply = (a, b) =>
  a.map((row) => [0, 1, 2].map((j) => row[0] * b[0][j] + row[1] * b[1][j] + row[2] * b[2][j]));

export function lambertPhase(alpha) {
  return (Math.sin(alpha) + (Math.PI - alpha) * Math.cos(alpha)) / Math.PI;
}

export class SkyBodies {
  constructor(product) {
    if (product?.schema !== 'terluna.illumination.site-sky/1')
      throw new Error('Unexpected site-sky schema: ' + product?.schema);
    const k = product.constants,
      s = product.site;
    this.product = product;
    this.site = s;
    this.period = k.synodic_month_s;
    this.omega = {
      synodic: TAU / k.synodic_month_s,
      sidereal: TAU / k.sidereal_month_s,
      anomalistic: TAU / k.anomalistic_month_s,
      draconic: TAU / k.draconic_month_s,
      earth: TAU / k.earth_sidereal_day_s,
    };
    this.k = k;
    const lat = rad(s.latitude_deg),
      lon = rad(s.longitude_deg);
    this.basis = [
      [-Math.sin(lon), Math.cos(lon), 0],
      [-Math.sin(lat) * Math.cos(lon), -Math.sin(lat) * Math.sin(lon), Math.cos(lat)],
      unit(lat, lon),
    ];
    const c = Math.cos(k.obliquity_rad),
      e = Math.sin(k.obliquity_rad);
    this.eqToEcl = [
      [1, 0, 0],
      [0, c, e],
      [0, -e, c],
    ];
  }
  /** State at t seconds after local noon, in the product's east-north-up frame. */
  state(t) {
    const s = this.site,
      k = this.k,
      w = this.omega,
      lon = rad(s.longitude_deg);
    const sun = unit(0, lon - w.synodic * t);
    const earth = s.libration
      ? unit(
          k.libration_latitude_rad * Math.sin(w.draconic * t + rad(s.libration_phase_latitude_deg)),
          k.libration_longitude_rad *
            Math.sin(w.anomalistic * t + rad(s.libration_phase_longitude_deg)),
        )
      : [1, 0, 0];
    const cosSE = Math.max(-1, Math.min(1, dot(sun, earth))),
      alpha = Math.acos(-cosSE);
    const theta = rad(s.sun_longitude_at_noon_deg) - lon + w.sidereal * t;
    const eclToBody = [
      [Math.cos(theta), Math.sin(theta), 0],
      [-Math.sin(theta), Math.cos(theta), 0],
      [0, 0, 1],
    ];
    return {
      sun: apply(this.basis, sun),
      earth: apply(this.basis, earth),
      earthFraction: (1 - cosSE) / 2,
      earthPhaseAngle: alpha,
      earthlightRatio:
        k.earth_geometric_albedo *
        (k.earth_radius_m / k.earth_moon_distance_m) ** 2 *
        lambertPhase(alpha),
      celestialToEnu: multiply(multiply(this.basis, eclToBody), this.eqToEcl),
      earthRotation: rad(s.earth_rotation_at_noon_deg) + w.earth * t,
    };
  }
  /** The same state with directions and the celestial matrix in the scene frame. */
  scene(t) {
    const st = this.state(t),
      toScene = ([e, n, u]) => [-e, u, n];
    const m = st.celestialToEnu;
    return {
      ...st,
      sun: toScene(st.sun),
      earth: toScene(st.earth),
      celestialToScene: [m[0].map((v) => -v), m[2], m[1]],
    };
  }
}
