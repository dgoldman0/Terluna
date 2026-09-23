/* The development cove's authored weather: clear and fog presets and a four-hour
 * episode. Forcing inputs for the engine's wetness and cloud display, not forecasts. */
import C from '../engine/core.js';
const { clamp, mix, smooth, cloudTau, emptyLedger, stepLedger } = C;
/* Authored four-hour episode. These are forcing inputs, not weather forecasts. */
export const EPISODE = [
  [0, 0.1, 0.015, 1100, 850, 1.5, 0.48, 294, 0, 35000],
  [1800, 0.25, 0.04, 1100, 900, 2.4, 0.59, 293, 0, 30000],
  [3600, 0.52, 0.08, 1050, 1000, 3.8, 0.72, 292, 0, 22000],
  [5100, 0.88, 0.2, 900, 1100, 5.4, 0.87, 290, 1.5, 8000],
  [6300, 0.97, 0.29, 800, 1200, 6.7, 0.96, 289, 7, 3800],
  [8100, 0.99, 0.34, 750, 1300, 6.2, 0.98, 288.5, 10, 2500],
  [9600, 0.92, 0.19, 900, 1150, 4.6, 0.94, 289, 3.5, 5000],
  [10800, 0.7, 0.12, 1100, 1000, 3.2, 0.87, 290, 0.3, 9500],
  [12600, 0.37, 0.05, 1350, 950, 2.6, 0.75, 292, 0, 22000],
  [14400, 0.12, 0.02, 1500, 900, 1.7, 0.64, 293, 0, 35000],
];
export function weatherAt(t, kind = 'episode') {
  let a, b;
  if (kind === 'clear') a = b = [0, 0, 0, 1100, 850, 1.5, 0.48, 294, 0, 180000];
  else if (kind === 'fog') a = b = [0, 0.62, 0.08, 800, 800, 1.2, 0.99, 288, 0, 110];
  else {
    t = clamp(t, 0, 14400);
    let j = 1;
    while (j < EPISODE.length - 1 && EPISODE[j][0] < t) j++;
    a = EPISODE[j - 1];
    b = EPISODE[j];
  }
  const q = a === b ? 0 : smooth(a[0], b[0], t),
    v = a.map((x, i) => mix(x, b[i], q));
  return {
    coverage: v[1],
    lwc: v[2],
    cloudBase: v[3],
    thickness: v[4],
    wind: v[5],
    humidity: v[6],
    temperature: v[7],
    rain: v[8],
    visibility: v[9],
    radius: 12,
    tau: cloudTau(v[2], v[4], 12),
  };
}
export function windDistance(t, kind = 'episode', step = 10) {
  let distance = 0;
  for (let a = 0; a < t; a += step) {
    const h = Math.min(step, t - a);
    distance += weatherAt(a + h * 0.5, kind).wind * h;
  }
  return distance;
}
export function ledgerAt(t, kind = 'episode', step = 10) {
  let l = emptyLedger();
  for (let s = 0; s < t; s += step) {
    const dt = Math.min(step, t - s);
    l = stepLedger(l, weatherAt(s + dt * 0.5, kind), dt);
  }
  return l;
}
