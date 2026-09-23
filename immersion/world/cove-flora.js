/* The development cove's planting rules: which trees, understory, grass and stones
 * go where, from the landscape's habitat fields. Authored suitability rules for engine
 * development, not a proposed ecology; the biosphere research should supply real ones. */
import C from '../engine/core.js';
import L from './landscape.js';
import { pathDistance, roofMask } from './cove-sites.js';

/* Trees stand on an 8 m lattice, at most one per cell. The detailed plan below
 * builds the cells in `detailed` as full tree models; the engine's far field draws
 * every other cell with the same rule (treeAt), out toward the horizon. */
export const TREE_LATTICE = {
  spacing: 8,
  originX: 0,
  originZ: 7,
  detailed: { i0: -15, i1: 17, j0: 0, j1: 17 },
};
// Outside the landscape's field grid the habitat rule never gives a tree probability
// above 0.45 * (1 - 0.42 * 0.8) + 0.025 (woodland under the default exposure, no
// wet margin), so a cell whose acceptance draw exceeds it is rejected unsampled.
const OPEN_COUNTRY_MAX_PROBABILITY = 0.45 * (1 - 0.42 * 0.8) + 0.025;

// Mean crown area of the rule's trees (m²): family 1 inland, 76% mature at 7.7-12.2 m
// and 24% juvenile at 2.1-5.2 m, crown radius 0.34 of height.
const MEAN_CROWN_AREA = Math.PI * 0.34 ** 2 * (0.76 * 100.7 + 0.24 * 14.1);

/** Expected canopy cover (fraction of ground under crowns) of the tree rule, from a habitat
 * classification (landscape.classify) at surface height h and slope. */
export function canopyCover(f, h, slope) {
  if (h < 1.2 || slope > 0.52 || f.soilDepth < 0.16) return 0;
  const p = 0.45 * f.community[1] + 0.25 * f.community[2] + 0.025 * f.community[0];
  return Math.min(1, (p * MEAN_CROWN_AREA) / (TREE_LATTICE.spacing * TREE_LATTICE.spacing));
}

/** The tree in lattice cell (i, j), or null. Deterministic per cell: the rule reads the
 * habitat before planting, as the detailed plan does, whatever canopy is set. */
export function treeAt(i, j) {
  const r = C.rng(Math.floor(L.hash(i, j, 873) * 4294967295)),
    x = (i + (r() - 0.5) * 0.8) * 8,
    z = 7 + (j + (r() - 0.5) * 0.8) * 8,
    accept = r();
  if (accept > OPEN_COUNTRY_MAX_PROBABILITY && !L.gridCoverage(x, z)) return null;
  if (L.height(x, z) < 1.2) return null;
  const f = L.sample(x, z, { canopy: false });
  if (
    pathDistance(x, z) < 3.7 ||
    roofMask(x, z, 4) ||
    f.elevation < 1.2 ||
    f.slope > 0.52 ||
    f.soilDepth < 0.16
  )
    return null;
  const probability = 0.45 * f.community[1] + 0.25 * f.community[2] + 0.025 * f.community[0];
  if (accept > probability) return null;
  const family = f.moisture > 0.45 ? 2 : f.exposure > 0.72 && f.elevation < 5 ? 0 : 1;
  const mature = r() > 0.24,
    h = mature ? (family === 0 ? 5.5 : 7.7) + r() * 4.5 : 2.1 + r() * 3.1;
  return {
    id: `tree:${i}:${j}`,
    x,
    z,
    y: C.groundHeight(x, z),
    height: h,
    family,
    ageClass: mature ? 'mature' : 'juvenile',
    crownRadius: h * (family === 1 ? 0.34 : 0.39),
    canopyOpacity: family === 2 ? 0.85 : 0.8,
    seed: Math.floor(r() * 4294967295),
    habitat: {
      soilDepth: f.soilDepth,
      moisture: f.moisture,
      exposure: f.exposure,
      suitability: probability,
    },
  };
}

export function plan() {
  L.initialize();
  L.setCanopies([]);
  const trees = [],
    rocks = [],
    gravel = [],
    tufts = [],
    understory = [];
  const { i0, i1, j0, j1 } = TREE_LATTICE.detailed;
  for (let j = j0; j <= j1; j++)
    for (let i = i0; i <= i1; i++) {
      const tree = treeAt(i, j);
      if (tree) trees.push(tree);
    }
  // Two managed edge trees are part of the authored planting scenario. Their
  // improved rooting zones occur in the shared substrate/material field.
  for (const [id, x, z, h, seed] of [
    ['managed:west-edge', -10, -1, 7.2, 381712],
    ['managed:east-edge', 18, 0, 7.2, 292013],
  ]) {
    const f = L.sample(x, z);
    if (f.elevation > 0.7 && f.soilDepth > 0.3)
      trees.push({
        id,
        x,
        z,
        y: C.groundHeight(x, z),
        height: h,
        family: 0,
        ageClass: 'mature',
        crownRadius: h * 0.39,
        canopyOpacity: 0.85,
        seed,
        managed: true,
        habitat: {
          soilDepth: f.soilDepth,
          moisture: f.moisture,
          exposure: f.exposure,
          suitability: f.community[1],
        },
      });
  }
  // A stable identifier priority plus crown-scale spacing, independent of render LOD.
  trees.sort((a, b) => a.seed - b.seed);
  const accepted = [];
  for (const t of trees)
    if (
      !accepted.some(
        (a) => Math.hypot(a.x - t.x, a.z - t.z) < (a.crownRadius + t.crownRadius) * 0.55,
      )
    )
      accepted.push(t);
  L.setCanopies(accepted);
  for (let j = -12; j < 34; j++)
    for (let i = -33; i < 34; i++) {
      const r = C.rng(Math.floor(L.hash(i, j, 218) * 4294967295)),
        x = (i + r()) * 4,
        z = (j + r()) * 4,
        f = L.sample(x, z);
      if (pathDistance(x, z) < 1.8 || roofMask(x, z, 1) || f.elevation < -1.0 || f.elevation > 24)
        continue;
      if (r() < f.weights[2] * 0.55 + f.weights[1] * 0.12) {
        const s = 0.32 + r() ** 1.3 * (1.4 + f.substrate * 2.4);
        if (pathDistance(x, z) < 2.0 + s * 1.2) continue;
        rocks.push({
          id: `rock:${i}:${j}`,
          x,
          z,
          s,
          angle: r() * C.TAU,
          family: Math.floor(r() * 6),
          stretch: 0.9 + r() * 0.7,
        });
      }
    }
  for (let j = -45; j < 92; j++)
    for (let i = -130; i < 130; i++) {
      const r = C.rng(Math.floor(L.hash(i, j, 94) * 4294967295)),
        x = (i + r()) * 0.5,
        z = (j + r()) * 0.5,
        f = L.sample(x, z);
      if (
        f.elevation < -0.55 ||
        f.elevation > 3 ||
        roofMask(x, z) ||
        r() > f.weights[1] * 0.8 + f.drainage * 0.18
      )
        continue;
      gravel.push({
        id: `pebble:${i}:${j}`,
        x,
        z,
        s: 0.018 + r() ** 3 * 0.09,
        angle: r() * C.TAU,
        tone: r(),
      });
    }
  for (let j = -2; j < 97; j++)
    for (let i = -83; i < 84; i++) {
      const r = C.rng(Math.floor(L.hash(i, j, 512) * 4294967295)),
        x = (i + r()) * 1.5,
        z = (j + r()) * 1.5;
      if (pathDistance(x, z) < 1.55 || roofMask(x, z, 1)) continue;
      const f = L.sample(x, z);
      const community = f.community[2] > 0.2 ? 2 : f.canopy > 0.35 ? 1 : 0;
      const density =
        (f.community[0] * 0.65 + f.community[2] * 0.6 + f.community[3] * 0.18) *
        (0.45 + 0.55 * L.noise(x * 0.085, z * 0.085, 91));
      if (r() < density)
        tufts.push({
          id: `tuft:${i}:${j}`,
          x,
          z,
          h: (community === 2 ? 0.38 : 0.19) + r() * 0.3,
          angle: r() * C.TAU,
          scale: 0.6 + r() * 0.9,
          community,
          tone: r(),
        });
      if (
        f.canopy > 0.28 &&
        f.soilDepth > 0.2 &&
        f.elevation > 1.8 &&
        r() < (0.035 + 0.12 * f.moisture) * f.canopy
      )
        understory.push({
          id: `understory:${i}:${j}`,
          x,
          z,
          h: 0.22 + r() * 0.4,
          angle: r() * C.TAU,
          scale: 0.5 + r() * 0.75,
        });
    }
  return { trees: accepted, rocks, gravel, tufts, understory, seed: 873, version: 'community-2' };
}
