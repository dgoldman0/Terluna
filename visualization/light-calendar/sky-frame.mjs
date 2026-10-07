// Orientation of the stars and of Earth's surface in an observer's local sky.
// The calendar evaluator gives the Sun and Earth directions; this module turns
// the ecliptic of date and Earth's rotating surface into the same local frame so
// the explorer can place stars and show which side of Earth faces the observer.
// The lunar rotation uses the evaluator's own mean arguments (Meeus ch. 47 and 53,
// without nutation or physical libration); Earth's rotation uses the IAU 1982
// mean sidereal time with UT1 taken as UTC.
import { julianDate, ttOffset } from "./model/evaluator.mjs";

const rad = Math.PI / 180;

const rz = (t) => [
  [Math.cos(t), -Math.sin(t), 0],
  [Math.sin(t), Math.cos(t), 0],
  [0, 0, 1],
];
const rx = (t) => [
  [1, 0, 0],
  [0, Math.cos(t), -Math.sin(t)],
  [0, Math.sin(t), Math.cos(t)],
];
export const mul = (a, b) =>
  a.map((row) => [0, 1, 2].map((j) => row.reduce((s, v, k) => s + v * b[k][j], 0)));
export const apply = (m, v) => m.map((row) => row[0] * v[0] + row[1] * v[1] + row[2] * v[2]);
export const transpose = (m) => [0, 1, 2].map((i) => m.map((row) => row[i]));

// Local frame rows: east, north, up at the observer (as in evaluator geometry()).
export function localFrame(longitude, latitude) {
  const lon = longitude * rad,
    lat = latitude * rad;
  return [
    [-Math.sin(lon), Math.cos(lon), 0],
    [-Math.sin(lat) * Math.cos(lon), -Math.sin(lat) * Math.sin(lon), Math.cos(lat)],
    [Math.cos(lat) * Math.cos(lon), Math.cos(lat) * Math.sin(lon), Math.sin(lat)],
  ];
}

// Ecliptic-of-date direction (from the Moon toward a target) to selenographic.
export function eclipticToMoon(jdTT, data) {
  const t = (jdTT - 2451545) / 36525;
  const poly = (a, b, c, d, e) => rad * (a + t * (b + t * (c + t * (d + t * e))));
  const F = poly(93.272095, 483202.0175233, -0.0036539, -1 / 3526000, 1 / 863310000);
  const node = poly(125.0445479, -1934.1362891, 0.0020754, 1 / 467441, -1 / 60616000);
  return mul(rz(Math.PI - F), mul(rx(data.equator_to_ecliptic), rz(-node)));
}

export function obliquity(jdTT) {
  const t = (jdTT - 2451545) / 36525;
  return (23.439291 - 0.0130042 * t) * rad;
}

export function siderealTime(jdUT) {
  const d = jdUT - 2451545,
    t = d / 36525;
  const deg = 280.46061837 + 360.98564736629 * d + 0.000387933 * t * t - (t * t * t) / 38710000;
  return (((deg % 360) + 360) % 360) * rad;
}

// Matrices into the observer's local frame (east, north, up).
export function skyFrame(ms, scale, longitude, latitude, data) {
  const jdTT = julianDate(ms, scale);
  const utcMs = scale === "TT" ? ms - ttOffset(ms) * 1000 : ms;
  const jdUT = julianDate(utcMs, "TT");
  const moonToLocal = localFrame(longitude, latitude);
  const eclToLocal = mul(moonToLocal, eclipticToMoon(jdTT, data));
  const eps = obliquity(jdTT);
  // Earth-fixed -> equatorial of date (sidereal rotation) -> ecliptic of date.
  const earthToEcl = mul(rx(-eps), rz(siderealTime(jdUT)));
  return {
    jdTT,
    eclToLocal,
    equatorialToLocal: mul(eclToLocal, rx(-eps)),
    earthToLocal: mul(eclToLocal, earthToEcl),
    localToEarth: transpose(mul(eclToLocal, earthToEcl)),
  };
}

export function localVector(elevationDeg, azimuthDeg) {
  const e = elevationDeg * rad,
    a = azimuthDeg * rad;
  return [Math.cos(e) * Math.sin(a), Math.cos(e) * Math.cos(a), Math.sin(e)];
}

export function lonLat(v) {
  return [Math.atan2(v[1], v[0]) / rad, Math.asin(Math.max(-1, Math.min(1, v[2]))) / rad];
}

// The point on Earth at the centre of its disk, as seen from the observer.
export function subObserverPoint(frame, earthElevationDeg, earthAzimuthDeg) {
  const toward = localVector(earthElevationDeg, earthAzimuthDeg);
  return lonLat(apply(frame.localToEarth, toward.map((x) => -x)));
}

// J2000 star direction to the local frame: J2000 equator -> J2000 ecliptic,
// general precession in longitude to the date, then the Moon's orientation.
export function starPrecession(jdTT) {
  const t = (jdTT - 2451545) / 36525;
  const eps0 = 23.4392911 * rad,
    p = ((5029.0966 * t + 1.11113 * t * t) / 3600) * rad;
  return mul(rz(p), rx(-eps0));
}
