/* The development cove: a placeholder world for building and testing the engine.
 *
 * Its landscape, viewpoints, rain shelter, path, weather episode and planting rules
 * were authored for development and make no claim about any Open Moon environment.
 * Environments grounded in the geography, climate and biosphere research will be
 * defined the same way, as world objects the experience registers on OM.world.
 */
import landscape from './landscape.js';
import { SHELTER, WAYPOINTS, pathDistance, roofMask } from './cove-sites.js';
import { EPISODE, ledgerAt, weatherAt, windDistance } from './cove-weather.js';
import { plan } from './cove-flora.js';

export default {
  id: 'development-cove',
  status: 'placeholder',
  landscape,
  waypoints: WAYPOINTS,
  shelter: SHELTER,
  roofMask,
  pathDistance,
  weather: { EPISODE, at: weatherAt, windDistance, ledgerAt },
  flora: plan,
};
