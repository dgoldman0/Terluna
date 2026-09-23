/* The development cove's viewpoints, rain shelter and path. Authored for engine
 * development; not a proposed Open Moon site. */
const clamp = (x, a = 0, b = 1) => Math.max(a, Math.min(b, x));
export const SHELTER = { x: 20, z: 33, width: 11, depth: 8, roofHeight: 4.3 };
export function roofMask(x, z, pad = 0) {
  return Math.abs(x - SHELTER.x) < SHELTER.width / 2 + pad &&
    Math.abs(z - SHELTER.z) < SHELTER.depth / 2 + pad
    ? 1
    : 0;
}
export const WAYPOINTS = [
  { id: 'shore', name: 'Shoreline', x: 0, z: 9, yaw: 0, pitch: -0.04 },
  { id: 'shelter', name: 'Rain shelter', x: 20, z: 33, yaw: 0.05, pitch: -0.06 },
  { id: 'path', name: 'Woodland path', x: -21, z: 39, yaw: -0.7, pitch: -0.06 },
  { id: 'overlook', name: 'High overlook', x: -56, z: 47, yaw: 0.04, pitch: -0.12 },
];
export function pathDistance(x, z) {
  let best = 1e9;
  const a = [
    [0, -18],
    [0, 9],
    [7, 20],
    [20, 33],
    [-6, 43],
    [-30, 43],
    [-56, 47],
  ];
  for (let i = 1; i < a.length; i++) {
    const [ax, az] = a[i - 1],
      [bx, bz] = a[i],
      dx = bx - ax,
      dz = bz - az,
      t = clamp(((x - ax) * dx + (z - az) * dz) / (dx * dx + dz * dz));
    best = Math.min(best, Math.hypot(x - ax - t * dx, z - az - t * dz));
  }
  return best;
}
