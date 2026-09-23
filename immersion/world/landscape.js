/* Shared landscape model, revision 2. Coordinates and heights are metres.
 * Authored engineered landscape; drainage is derived from this surface.
 * No climate, erosion-time, or ecological-equilibrium prediction is implied.
 * This module is independent of Three.js and also runs under Node.
 */
const clamp = (x, a = 0, b = 1) => Math.max(a, Math.min(b, x));
const mix = (a, b, t) => a + (b - a) * t;
const smooth = (a, b, x) => {
  const t = clamp((x - a) / (b - a));
  return t * t * (3 - 2 * t);
};
function hash(x, z, s = 0) {
  let n = Math.imul(x | 0, 374761393) + Math.imul(z | 0, 668265263) + Math.imul(s | 0, 1442695041);
  n = Math.imul(n ^ (n >>> 13), 1274126177);
  return ((n ^ (n >>> 16)) >>> 0) / 4294967295;
}
function noise(x, z, s = 0) {
  const i = Math.floor(x),
    j = Math.floor(z);
  let u = x - i,
    v = z - j;
  u = u * u * (3 - 2 * u);
  v = v * v * (3 - 2 * v);
  return mix(
    mix(hash(i, j, s), hash(i + 1, j, s), u),
    mix(hash(i, j + 1, s), hash(i + 1, j + 1, s), u),
    v,
  );
}
function fbm(x, z, s = 0) {
  return (
    0.58 * noise(x, z, s) +
    0.28 * noise(x * 2.03 + 17.4, z * 2.03 + 8.1, s) +
    0.14 * noise(x * 4.13 + 31.7, z * 4.13 + 19.6, s)
  );
}
function bell(x, z, cx, cz, sx, sz) {
  return Math.exp(-(((x - cx) / sx) ** 2 + ((z - cz) / sz) ** 2));
}
function segmentDistance(x, z, a, b) {
  const dx = b[0] - a[0],
    dz = b[1] - a[1],
    t = clamp(((x - a[0]) * dx + (z - a[1]) * dz) / (dx * dx + dz * dz));
  return { d: Math.hypot(x - a[0] - t * dx, z - a[1] - t * dz), t };
}
// Connected shallow valley centre lines, selected to define the basin layout.
// The contributing area, depression fill, receivers and outlets are computed.
const VALLEYS = [
  {
    points: [
      [-42, 130],
      [-27, 85],
      [-18, 57],
      [-11, 34],
      [-7, 13],
      [-3, -10],
      [4, -45],
    ],
    width: 3.5,
    cut: 0.4,
  },
  {
    points: [
      [57, 124],
      [43, 93],
      [28, 62],
      [-11, 34],
    ],
    width: 3.0,
    cut: 0.22,
  },
  {
    points: [
      [-74, 112],
      [-60, 83],
      [-27, 85],
    ],
    width: 2.5,
    cut: 0.28,
  },
];
function valley(x, z) {
  let v = 0;
  for (const line of VALLEYS) {
    let d = Infinity;
    for (let i = 1; i < line.points.length; i++)
      d = Math.min(d, segmentDistance(x, z, line.points[i - 1], line.points[i]).d);
    v += line.cut * Math.exp(-((d / line.width) ** 2));
  }
  return v;
}
function substrateAffinity(x, z) {
  const local =
    0.95 * bell(x, z, -17, -7, 13, 21) +
    0.94 * bell(x, z, 23, -11, 16, 24) +
    0.84 * bell(x, z, -55, 45, 25, 29) +
    0.7 * bell(x, z, 94, 16, 29, 40);
  return clamp(local + 0.33 * smooth(0.55, 0.79, fbm(x * 0.025 + 9, z * 0.025 - 7, 34)));
}
function regionalHeight(x, z) {
  const warpX = (fbm(x * 0.0007, z * 0.0007, 83) - 0.5) * 165,
    warpZ = (fbm(x * 0.0006 + 40, z * 0.0006, 12) - 0.5) * 140;
  const qx = x + warpX,
    qz = z + warpZ;
  const ranges =
    205 * bell(qx, qz, -950, -1830, 670, 910) +
    258 * bell(qx, qz, 1270, -2720, 890, 940) +
    92 * bell(qx, qz, 2430, -2080, 780, 970);
  // Connected anisotropic ridges and erosional folds carried by the same field.
  const folded =
    0.47 +
    0.32 * fbm(qx * 0.008, qz * 0.004, 7) +
    0.49 * (1 - Math.abs(2 * noise(qx * 0.017 + qz * 0.007, qz * 0.011, 91) - 1)) ** 1.4;
  return ranges * folded + 38 * bell(x, z, -182, -265, 112, 236) * smooth(35, 115, -z);
}
function rawHeight(x, z) {
  if (!Number.isFinite(x) || !Number.isFinite(z))
    throw new TypeError('Landscape coordinates must be finite metres.');
  const coastal = 68 * Math.tanh((z + 11.5) / 940);
  const headlands =
    2.75 * bell(x, z, -48, -8, 28, 46) +
    2.4 * bell(x, z, 68, 9, 43, 52) -
    0.4 * bell(x, z, 0, -8, 24, 35);
  const overlook = 14 * bell(x, z, -57, 50, 27, 25);
  const inland = (fbm(x * 0.019 + 3, z * 0.013 - 9, 2) - 0.5) * 3.4 * smooth(8, 72, z);
  const regional = regionalHeight(x, z);
  const rock = substrateAffinity(x, z),
    land = smooth(-1.5, 2.5, coastal + headlands + overlook);
  const meso = (fbm(x * 0.38, z * 0.27, 33) - 0.5) * 0.26 * rock * land;
  const fractures =
    (Math.abs(noise(x * 1.25 + z * 0.38, z * 0.42, 29) - 0.5) - 0.25) * 0.15 * rock * land;
  const soilMicro = (noise(x * 1.7, z * 1.7, 101) - 0.5) * 0.028 * land;
  const channels = Math.abs(x) < 200 && z > -70 && z < 180 ? valley(x, z) : 0;
  return (
    coastal +
    headlands +
    overlook +
    inland +
    regional +
    meso +
    fractures +
    soilMicro -
    channels -
    0.56172690182407 * bell(x, z, 0, 9, 200, 200)
  );
}
class MinHeap {
  constructor() {
    this.a = [];
  }
  push(v) {
    const a = this.a;
    let i = a.length;
    a.push(v);
    while (i) {
      const p = (i - 1) >> 1;
      if (a[p].h < v.h || (a[p].h === v.h && a[p].id < v.id)) break;
      a[i] = a[p];
      i = p;
    }
    a[i] = v;
  }
  pop() {
    const a = this.a,
      first = a[0],
      last = a.pop();
    if (a.length) {
      let i = 0;
      while (2 * i + 1 < a.length) {
        let c = 2 * i + 1;
        if (
          c + 1 < a.length &&
          (a[c + 1].h < a[c].h || (a[c + 1].h === a[c].h && a[c + 1].id < a[c].id))
        )
          c++;
        if (last.h < a[c].h || (last.h === a[c].h && last.id < a[c].id)) break;
        a[i] = a[c];
        i = c;
      }
      a[i] = last;
    }
    return first;
  }
  get length() {
    return this.a.length;
  }
}
const NEIGHBOURS = [
  [-1, -1, Math.SQRT2],
  [0, -1, 1],
  [1, -1, Math.SQRT2],
  [-1, 0, 1],
  [1, 0, 1],
  [-1, 1, Math.SQRT2],
  [0, 1, 1],
  [1, 1, Math.SQRT2],
];
/** Priority flood supplies depression spill heights and an acyclic flat routing.
 * Sea cells and perimeter cells are explicit open boundaries. Contributing area
 * is a geometric routing diagnostic; rainfall/runoff volumes are separate.
 */
function drainageGrid(elevation, n, cell) {
  if (n < 3 || elevation.length !== n * n || !(cell > 0))
    throw new RangeError('Invalid drainage grid');
  const N = n * n,
    filled = new Float64Array(elevation),
    parent = new Int32Array(N).fill(-1),
    rank = new Int32Array(N).fill(-1),
    seen = new Uint8Array(N),
    heap = new MinHeap();
  for (let j = 0; j < n; j++)
    for (let i = 0; i < n; i++) {
      const id = j * n + i;
      if (i === 0 || j === 0 || i === n - 1 || j === n - 1 || elevation[id] <= 0) {
        seen[id] = 1;
        heap.push({ id, h: elevation[id] });
      }
    }
  let order = 0;
  while (heap.length) {
    const { id, h } = heap.pop();
    rank[id] = order++;
    const i = id % n,
      j = (id / n) | 0;
    for (const [dx, dz] of NEIGHBOURS) {
      const x = i + dx,
        z = j + dz;
      if (x < 0 || z < 0 || x >= n || z >= n) continue;
      const k = z * n + x;
      if (seen[k]) continue;
      seen[k] = 1;
      parent[k] = id;
      filled[k] = Math.max(elevation[k], h);
      heap.push({ id: k, h: filled[k] });
    }
  }
  const receiver = new Int32Array(N).fill(-1),
    incoming = new Int32Array(N),
    accumulation = new Float64Array(N).fill(cell * cell);
  for (let id = 0; id < N; id++) {
    const i = id % n,
      j = (id / n) | 0;
    if (i === 0 || j === 0 || i === n - 1 || j === n - 1 || elevation[id] <= 0) continue;
    let best = -1,
      bestSlope = 0;
    for (const [dx, dz, dist] of NEIGHBOURS) {
      const k = (j + dz) * n + i + dx,
        s = (filled[id] - filled[k]) / (cell * dist);
      if (s > bestSlope + 1e-13) {
        bestSlope = s;
        best = k;
      }
    }
    if (best < 0) best = parent[id];
    receiver[id] = best;
    if (best >= 0) incoming[best]++;
  }
  const queue = new Int32Array(N);
  let head = 0,
    tail = 0;
  for (let i = 0; i < N; i++) if (!incoming[i]) queue[tail++] = i;
  while (head < tail) {
    const id = queue[head++],
      r = receiver[id];
    if (r >= 0) {
      accumulation[r] += accumulation[id];
      if (--incoming[r] === 0) queue[tail++] = r;
    }
  }
  if (tail !== N) throw Error('Drainage receivers contain a cycle.');
  let outletArea = 0;
  for (let i = 0; i < N; i++) if (receiver[i] < 0) outletArea += accumulation[i];
  return {
    filled,
    receiver,
    accumulation,
    rank,
    order: queue,
    outletArea,
    totalArea: N * cell * cell,
    processed: tail,
  };
}
const GRID = { n: 257, cell: 2, minX: -256, minZ: -128 };
let data = null,
  canopies = [];
function gridSample(array, x, z, meta = GRID) {
  const u = clamp((x - meta.minX) / meta.cell, 0, meta.n - 1),
    v = clamp((z - meta.minZ) / meta.cell, 0, meta.n - 1),
    i = Math.min(meta.n - 2, Math.floor(u)),
    j = Math.min(meta.n - 2, Math.floor(v)),
    a = j * meta.n + i,
    tx = u - i,
    tz = v - j;
  return mix(
    mix(array[a], array[a + 1], tx),
    mix(array[a + meta.n], array[a + meta.n + 1], tx),
    tz,
  );
}
function gridChannel(array, x, z, channel) {
  const { n, cell, minX, minZ } = GRID,
    u = clamp((x - minX) / cell, 0, n - 1),
    v = clamp((z - minZ) / cell, 0, n - 1),
    i = Math.min(n - 2, Math.floor(u)),
    j = Math.min(n - 2, Math.floor(v)),
    tx = u - i,
    tz = v - j;
  const at = (a, b) => array[(b * n + a) * 4 + channel];
  return mix(mix(at(i, j), at(i + 1, j), tx), mix(at(i, j + 1), at(i + 1, j + 1), tx), tz);
}
function gridCoverage(x, z) {
  return (
    smooth(0, 12, x - GRID.minX) *
    smooth(0, 12, GRID.minX + (GRID.n - 1) * GRID.cell - x) *
    smooth(0, 12, z - GRID.minZ) *
    smooth(0, 12, GRID.minZ + (GRID.n - 1) * GRID.cell - z)
  );
}
function height(x, z) {
  if (!data) initialize();
  return (
    rawHeight(x, z) -
    (data ? gridSample(data.incision, x, z) * gridCoverage(x, z) : 0) +
    (data?.referenceCorrection || 0) * bell(x, z, 0, 9, 18, 18)
  );
}
function gradient(x, z, epsilon = 0.5) {
  return {
    x: (height(x + epsilon, z) - height(x - epsilon, z)) / (2 * epsilon),
    z: (height(x, z + epsilon) - height(x, z - epsilon)) / (2 * epsilon),
  };
}
function initialize() {
  if (data) return data;
  const { n, cell, minX, minZ } = GRID,
    N = n * n,
    base = new Float64Array(N);
  for (let j = 0; j < n; j++)
    for (let i = 0; i < n; i++) base[j * n + i] = rawHeight(minX + i * cell, minZ + j * cell);
  const first = drainageGrid(base, n, cell),
    cut = new Float32Array(N);
  for (let i = 0; i < N; i++)
    cut[i] = 0.27 * smooth(150, 1800, first.accumulation[i]) * smooth(0.15, 2.0, base[i]);
  const incision = new Float32Array(N);
  for (let j = 1; j < n - 1; j++)
    for (let i = 1; i < n - 1; i++) {
      const k = j * n + i;
      incision[k] = (cut[k] * 4 + cut[k - 1] + cut[k + 1] + cut[k - n] + cut[k + n]) / 8;
    }
  data = { incision };
  data.referenceCorrection = 1.3458966316214787 - height(0, 9);
  const elevation = new Float64Array(N);
  for (let j = 0; j < n; j++)
    for (let i = 0; i < n; i++) elevation[j * n + i] = height(minX + i * cell, minZ + j * cell);
  const hydrology = drainageGrid(elevation, n, cell),
    moisture = new Float32Array(N),
    canopy = new Float32Array(N),
    material = new Float32Array(N * 4),
    habitat = new Float32Array(N * 4),
    soil = new Float32Array(N),
    exposure = new Float32Array(N),
    slope = new Float32Array(N),
    curvature = new Float32Array(N);
  Object.assign(data, {
    ...GRID,
    elevation,
    hydrology,
    moisture,
    canopy,
    material,
    habitat,
    soil,
    exposure,
    slope,
    curvature,
  });
  for (let j = 0; j < n; j++)
    for (let i = 0; i < n; i++) {
      const k = j * n + i,
        x = minX + i * cell,
        z = minZ + j * cell,
        h = elevation[k],
        il = Math.max(0, i - 1),
        ir = Math.min(n - 1, i + 1),
        jl = Math.max(0, j - 1),
        jr = Math.min(n - 1, j + 1);
      const dx = (elevation[j * n + ir] - elevation[j * n + il]) / ((ir - il) * cell),
        dz = (elevation[jr * n + i] - elevation[jl * n + i]) / ((jr - jl) * cell),
        s = Math.hypot(dx, dz);
      const cv =
        (elevation[j * n + ir] +
          elevation[j * n + il] +
          elevation[jr * n + i] +
          elevation[jl * n + i] -
          4 * h) /
        (cell * cell);
      const route = smooth(100, 2200, hydrology.accumulation[k]);
      const ex = clamp(
        0.72 + 0.45 * s - 0.3 * smooth(0, 3, elevation[j * n + Math.max(0, i - 6)] - h),
      );
      const wet = clamp(
        0.2 + 0.43 * route + 0.2 * smooth(0, 0.045, cv) + 0.16 * (1 - smooth(0, 4, h)) - 0.2 * s,
      );
      const f = classify(x, z, h, s, wet, route, 0, ex);
      moisture[k] = wet;
      soil[k] = f.soilDepth;
      exposure[k] = ex;
      slope[k] = s;
      curvature[k] = cv;
      material.set(f.weights, k * 4);
      habitat.set(f.community, k * 4);
    }
  return data;
}
function classify(x, z, h, s, wet, route, cover, ex) {
  const plantedSoil =
    clamp(0.92 * bell(x, z, -10, 3, 9, 17) + 0.88 * bell(x, z, 19, 7, 12, 19)) *
    smooth(0.5, 1.15, h);
  const affinity = substrateAffinity(x, z),
    rock =
      clamp(affinity * 0.85 + smooth(0.2, 0.65, s) * 0.8 - 0.2 * route) * (1 - 0.52 * plantedSoil);
  const beach = (1 - smooth(1.15, 3.0, h)) * (1 - smooth(0.12, 0.38, s)) * (1 - 0.88 * plantedSoil);
  const gravel = beach * clamp(0.13 + 0.6 * rock + 0.14 * route),
    sand = beach * (1 - rock) * (1 - 0.45 * route),
    stone = rock * (1 - 0.35 * beach),
    loam = (1 - beach) * (1 - rock);
  const sum = Math.max(1e-9, sand + gravel + stone + loam),
    weights = [sand / sum, gravel / sum, stone / sum, loam / sum];
  const soilDepth = clamp(
    0.12 + 0.75 * weights[3] + 0.3 * route + 0.38 * plantedSoil - 0.33 * rock,
    0.015,
    1.2,
  );
  const viable = smooth(0.55, 1.65, h) * (1 - smooth(0.4, 0.8, s));
  const woodland =
    viable *
    smooth(0.13, 0.48, soilDepth) *
    Math.max(smooth(1.8, 4.7, h), plantedSoil * 0.95) *
    (1 - 0.42 * ex);
  const margin = viable * smooth(0.4, 0.78, wet) * (1 - 0.6 * cover);
  const grass = viable * (1 - 0.68 * woodland) * (1 - 0.88 * cover) * (1 - 0.7 * rock);
  const sparse = viable * rock * (1 - 0.65 * cover);
  return { weights, soilDepth, plantedSoil, community: [grass, woodland, margin, sparse] };
}
function sample(x, z) {
  initialize();
  const h = height(x, z),
    g = gradient(x, z, 1),
    s = Math.hypot(g.x, g.z),
    coverage = gridCoverage(x, z);
  const wet = mix(0.24, gridSample(data.moisture, x, z), coverage),
    route = smooth(100, 2200, gridSample(data.hydrology.accumulation, x, z)) * coverage;
  const canopy = gridSample(data.canopy, x, z) * coverage,
    ex = mix(0.8, gridSample(data.exposure, x, z), coverage);
  const f = classify(x, z, h, s, wet, route, canopy, ex);
  f.weights = f.weights.map((v, k) => mix(v, gridChannel(data.material, x, z, k), coverage));
  f.soilDepth = mix(f.soilDepth, gridSample(data.soil, x, z), coverage);
  const lightFactors = [1 - 0.88 * canopy, 1, 1 - 0.6 * canopy, 1 - 0.65 * canopy];
  f.community = f.community.map((v, k) =>
    mix(v, gridChannel(data.habitat, x, z, k) * lightFactors[k], coverage),
  );
  return {
    elevation: h,
    gradient: g,
    slope: s,
    curvature: gridSample(data.curvature, x, z) * coverage,
    exposure: ex,
    moisture: wet,
    drainage: route,
    contributingArea: gridSample(data.hydrology.accumulation, x, z) * coverage,
    canopy,
    lightAvailability: 1 - 0.78 * canopy,
    substrate: substrateAffinity(x, z),
    ...f,
  };
}
function setCanopies(trees) {
  initialize();
  canopies = trees.map((t) => ({
    x: t.x,
    z: t.z,
    r: t.crownRadius,
    opacity: t.canopyOpacity || 0.8,
  }));
  data.canopy.fill(0);
  for (const t of canopies) {
    const r = t.r * 1.5,
      minI = Math.max(0, Math.floor((t.x - r - GRID.minX) / GRID.cell)),
      maxI = Math.min(GRID.n - 1, Math.ceil((t.x + r - GRID.minX) / GRID.cell)),
      minJ = Math.max(0, Math.floor((t.z - r - GRID.minZ) / GRID.cell)),
      maxJ = Math.min(GRID.n - 1, Math.ceil((t.z + r - GRID.minZ) / GRID.cell));
    for (let j = minJ; j <= maxJ; j++)
      for (let i = minI; i <= maxI; i++) {
        const x = GRID.minX + i * GRID.cell,
          z = GRID.minZ + j * GRID.cell,
          q = ((x - t.x) / t.r) ** 2 + ((z - t.z) / t.r) ** 2,
          k = j * GRID.n + i;
        data.canopy[k] = 1 - (1 - data.canopy[k]) * Math.exp(-t.opacity * 1.8 * Math.exp(-q * 1.8));
      }
  }
  return data.canopy;
}
/** Local shoreline root query retained for UI/audio compatibility. The terrain
 * remains fully two-dimensional and supports islands and multiple crossings.
 */
function shorelineAtX(x, z0 = -100, z1 = 60) {
  let a = z0,
    b = z1,
    fa = height(x, a),
    fb = height(x, b);
  for (let n = 0; n < 10 && fa > 0; n++) {
    a -= Math.max(100, b - a);
    fa = height(x, a);
  }
  for (let n = 0; n < 10 && fb < 0; n++) {
    b += Math.max(100, b - a);
    fb = height(x, b);
  }
  if (fa > 0 || fb < 0) throw new RangeError('No bracketed coastal crossing at this longitude.');
  for (let i = 0; i < 40; i++) {
    const m = (a + b) * 0.5;
    if (height(x, m) > 0) b = m;
    else a = m;
  }
  return (a + b) * 0.5;
}
function materialFieldRGBA() {
  initialize();
  const out = new Float32Array(GRID.n * GRID.n * 4);
  out.set(data.material);
  return out;
}
function environmentFieldRGBA() {
  initialize();
  const out = new Float32Array(GRID.n * GRID.n * 4);
  for (let k = 0; k < GRID.n * GRID.n; k++) {
    out[k * 4] = data.moisture[k];
    out[k * 4 + 1] = data.canopy[k];
    out[k * 4 + 2] = smooth(100, 2200, data.hydrology.accumulation[k]);
    out[k * 4 + 3] = data.soil[k];
  }
  return out;
}
const api = {
  version: 'landscape-2',
  GRID,
  VALLEYS,
  clamp,
  mix,
  smooth,
  hash,
  noise,
  fbm,
  rawHeight,
  height,
  gradient,
  initialize,
  sample,
  classify,
  setCanopies,
  shorelineAtX,
  gridSample,
  gridChannel,
  gridCoverage,
  drainageGrid,
  materialFieldRGBA,
  environmentFieldRGBA,
  get data() {
    return data;
  },
};
export default api;
