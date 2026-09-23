/* Eye-like exposure for the experience: a perceptual camera convention, not a model
 * of vision. Exposure follows the metered light as L^-0.85, calibrated so clear lunar
 * noon (94.2 klx) keeps the reference exposure, and adapts over real time: quickly
 * toward brighter scenes, more slowly toward darker ones, as eyes do.
 */
export const REFERENCE_EXPOSURE = 0.244647046;
export const REFERENCE_LUX = 94200;
export const ADAPTATION_POWER = 0.85;
export const MAX_EXPOSURE = 25000;
export const BRIGHTEN_SECONDS = 0.35;
export const DARKEN_SECONDS = 2.0;

/** The exposure an adapted eye would settle on for this much light (lux). */
export function targetExposure(lux) {
  const l = Number.isFinite(lux) ? Math.max(0.004, lux) : REFERENCE_LUX;
  return Math.min(MAX_EXPOSURE, REFERENCE_EXPOSURE * Math.pow(REFERENCE_LUX / l, ADAPTATION_POWER));
}

/** Move the current exposure toward the target over dt seconds, in log space. */
export function adapt(current, target, dt) {
  if (!(current > 0) || !Number.isFinite(current)) return target;
  const tau = target < current ? BRIGHTEN_SECONDS : DARKEN_SECONDS;
  const k = 1 - Math.exp(-Math.max(0, dt) / tau);
  return Math.exp(Math.log(current) + (Math.log(target) - Math.log(current)) * k);
}
