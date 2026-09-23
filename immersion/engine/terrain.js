/* Camera-following nested regular terrain grids. Fine edges collapse onto the
 * parent lattice; a 24%-wide transition interpolates the parent triangle.
 * Every pass renders these same CPU-displaced vertices, including scene depth,
 * reflections and shadows. Heights and collision share OpenMoonLandscape.
 */
import { OM } from './om.js';
import C from './core.js';
import L from '../world/landscape.js';
function triangle(a, b, c, d, u, v) {
  return u + v <= 1 ? a + (b - a) * u + (c - a) * v : d + (c - d) * (1 - u) + (b - d) * (1 - v);
}
class TerrainSystem {
  constructor(T, scene, material, world = 'moon', options = {}) {
    this.T = T;
    this.scene = scene;
    this.material = material;
    this.world = world;
    this.cells = options.cells || 128;
    this.levelCount = options.levels || 11;
    this.baseStep = options.step || 0.25;
    this.levels = [];
    this.updates = 0;
    this.lastUpdateMs = 0;
    this.maxUpdateMs = 0;
    this.heightEvaluations = 0;
    this.parentScratch = { h: 0, nx: 0, ny: 1, nz: 0 };
    if (
      !Number.isInteger(this.cells) ||
      this.cells < 16 ||
      this.cells % 4 ||
      !Number.isInteger(this.levelCount) ||
      this.levelCount < 1 ||
      this.levelCount > 16 ||
      !Number.isFinite(this.baseStep) ||
      this.baseStep <= 0
    )
      throw new RangeError(
        'Terrain needs a positive spacing, 1–16 levels and a cell count divisible by four.',
      );
    for (let level = 0; level < this.levelCount; level++) {
      const step = this.baseStep * 2 ** level,
        n = this.cells + 1,
        geo = new T.BufferGeometry();
      for (const [name, size] of [
        ['position', 3],
        ['normal', 3],
        ['omWeights', 4],
        ['omCanopy', 1],
      ])
        geo.setAttribute(
          name,
          new T.BufferAttribute(new Float32Array(n * n * size), size).setUsage(T.DynamicDrawUsage),
        );
      geo.setIndex(new T.BufferAttribute(new Uint32Array(this.cells * this.cells * 6), 1));
      const mesh = new T.Mesh(geo, material);
      mesh.receiveShadow = true;
      mesh.castShadow = level < 6;
      mesh.name = 'terrain-lod-' + level;
      mesh.userData.landscape = true;
      scene.add(mesh);
      this.levels.push({
        level,
        step,
        mesh,
        geo,
        cache: new Map(),
        key: '',
        originKey: '',
        topologies: new Map(),
        bounds: null,
      });
    }
    this.update(0, 9, true);
  }
  sample(level, x, z) {
    const entry = this.levels[Math.min(level, this.levelCount - 1)],
      s = this.baseStep * 2 ** level;
    const ix = Math.round(x / s),
      iz = Math.round(z / s),
      key = (ix + 131072) * 262145 + (iz + 131072);
    let val = entry.cache.get(key + level / 32);
    if (val) return val;
    const h = C.groundHeight(x, z, this.world),
      eps = Math.max(0.125, s * 0.42),
      g = C.surfaceGradient(x, z, eps),
      R = C.worldRadius(this.world),
      nx = -g.x + x / R,
      nz = -g.z + z / R,
      norm = Math.hypot(nx, 1, nz);
    const f = L.classify(
      x,
      z,
      h + C.curvatureSag(x, z, this.world),
      Math.hypot(g.x, g.z),
      0.24,
      0,
      0,
      0.8,
    );
    val = { h, nx: nx / norm, ny: 1 / norm, nz: nz / norm, weights: f.weights };
    entry.cache.set(key + level / 32, val);
    this.heightEvaluations++;
    return val;
  }
  parentSample(level, x, z) {
    const step = this.baseStep * 2 ** (level + 1),
      ix = Math.floor(x / step),
      iz = Math.floor(z / step),
      x0 = ix * step,
      z0 = iz * step,
      u = (x - x0) / step,
      v = (z - z0) / step;
    const a = this.sample(level + 1, x0, z0),
      b = this.sample(level + 1, x0 + step, z0),
      c = this.sample(level + 1, x0, z0 + step),
      d = this.sample(level + 1, x0 + step, z0 + step);
    const out = this.parentScratch;
    out.h = triangle(a.h, b.h, c.h, d.h, u, v);
    out.nx = triangle(a.nx, b.nx, c.nx, d.nx, u, v);
    out.ny = triangle(a.ny, b.ny, c.ny, d.ny, u, v);
    out.nz = triangle(a.nz, b.nz, c.nz, d.nz, u, v);
    return out;
  }
  update(x, z, force = false) {
    if (!Number.isFinite(x) || !Number.isFinite(z))
      throw new TypeError('Terrain observer coordinates must be finite.');
    const started = typeof performance !== 'undefined' ? performance.now() : Date.now(),
      N = this.cells,
      n = N + 1;
    let changed = false;
    for (const entry of this.levels) {
      const { level, step, geo } = entry,
        cx = Math.floor(x / (step * 2)) * step * 2,
        cz = Math.floor(z / (step * 2)) * step * 2,
        r = N * step * 0.5;
      const bounds = { minX: cx - r, maxX: cx + r, minZ: cz - r, maxZ: cz + r, cx, cz, r },
        child = level ? this.levels[level - 1].bounds : null;
      const key = [cx, cz, child?.minX, child?.minZ, this.world].join(':');
      entry.bounds = bounds;
      if (!force && key === entry.key) continue;
      changed = true;
      entry.key = key;
      const originKey = [cx, cz, this.world].join(':');
      const moved = force || originKey !== entry.originKey;
      entry.originKey = originKey;
      const p = geo.attributes.position.array,
        normal = geo.attributes.normal.array,
        weights = geo.attributes.omWeights.array;
      let minY = Infinity,
        maxY = -Infinity;
      if (moved) {
        for (let j = 0; j <= N; j++)
          for (let i = 0; i <= N; i++) {
            let wx = bounds.minX + i * step,
              wz = bounds.minZ + j * step;
            // Even parent lattice on each boundary. Degenerate triangles are omitted.
            if (level < this.levelCount - 1) {
              if (j === 0 || j === N) wx = Math.floor(wx / (2 * step)) * 2 * step;
              if (i === 0 || i === N) wz = Math.floor(wz / (2 * step)) * 2 * step;
            }
            const k = j * n + i,
              a = this.sample(level, wx, wz),
              alpha =
                level === this.levelCount - 1
                  ? 0
                  : C.smooth(0.76, 1, Math.max(Math.abs(wx - cx), Math.abs(wz - cz)) / r),
              b = alpha > 0 ? this.parentSample(level, wx, wz) : a;
            p[k * 3] = wx;
            p[k * 3 + 1] = C.mix(a.h, b.h, alpha);
            p[k * 3 + 2] = wz;
            minY = Math.min(minY, p[k * 3 + 1]);
            maxY = Math.max(maxY, p[k * 3 + 1]);
            if (alpha === 0) {
              normal[k * 3] = a.nx;
              normal[k * 3 + 1] = a.ny;
              normal[k * 3 + 2] = a.nz;
              weights.set(a.weights, k * 4);
              continue;
            }
            const nx = C.mix(a.nx, b.nx, alpha),
              ny = C.mix(a.ny, b.ny, alpha),
              nz = C.mix(a.nz, b.nz, alpha),
              len = Math.hypot(nx, ny, nz);
            normal[k * 3] = nx / len;
            normal[k * 3 + 1] = ny / len;
            normal[k * 3 + 2] = nz / len;
            weights.set(a.weights, k * 4);
          }
      }
      const topologyKey = child
        ? [(child.minX - bounds.minX) / step, (child.minZ - bounds.minZ) / step].join(':')
        : 'full';
      let topology = entry.topologies.get(topologyKey);
      if (!topology) {
        const indices = new Uint32Array(N * N * 6);
        let count = 0;
        const add = (a, b, c) => {
          const area =
            (p[b * 3 + 2] - p[a * 3 + 2]) * (p[c * 3] - p[a * 3]) -
            (p[b * 3] - p[a * 3]) * (p[c * 3 + 2] - p[a * 3 + 2]);
          if (area > 1e-9) {
            indices[count++] = a;
            indices[count++] = b;
            indices[count++] = c;
          }
        };
        for (let j = 0; j < N; j++)
          for (let i = 0; i < N; i++) {
            const wx = bounds.minX + i * step,
              wz = bounds.minZ + j * step;
            if (
              child &&
              wx >= child.minX &&
              wx + step <= child.maxX &&
              wz >= child.minZ &&
              wz + step <= child.maxZ
            )
              continue;
            const a = j * n + i;
            // The double-collapsed upper-right corner has an interior vertex on its
            // diagonal. Split that triangle explicitly to eliminate the T-junction.
            if (level < this.levelCount - 1 && i === N - 1 && j === N - 1) {
              add(a + 1, a, a + n + 1);
              add(a, a + n, a + n + 1);
            } else {
              add(a, a + n, a + 1);
              add(a + 1, a + n, a + n + 1);
            }
          }
        topology = indices.slice(0, count);
        entry.topologies.set(topologyKey, topology);
      }
      geo.index.array.set(topology);
      geo.setDrawRange(0, topology.length);
      geo.index.needsUpdate = true;
      if (moved) {
        for (const name of ['position', 'normal', 'omWeights'])
          geo.attributes[name].needsUpdate = true;
        if (!geo.boundingBox) geo.boundingBox = new this.T.Box3();
        if (!geo.boundingSphere) geo.boundingSphere = new this.T.Sphere();
        geo.boundingBox.min.set(bounds.minX, minY, bounds.minZ);
        geo.boundingBox.max.set(bounds.maxX, maxY, bounds.maxZ);
        geo.boundingBox.getBoundingSphere(geo.boundingSphere);
      }
      entry.activeTriangles = topology.length / 3;
      // Bound caches after long traversal; recently sampled lattice is retained.
      if (entry.cache.size > n * n * 5) {
        const keep = new Map();
        for (const [key, v] of entry.cache) {
          keep.set(key, v);
          if (keep.size > n * n * 3) keep.delete(keep.keys().next().value);
        }
        entry.cache = keep;
      }
    }
    if (changed) {
      this.updates++;
      this.lastUpdateMs =
        (typeof performance !== 'undefined' ? performance.now() : Date.now()) - started;
      this.maxUpdateMs = Math.max(this.maxUpdateMs, this.lastUpdateMs);
    }
    return changed;
  }
  setWorld(world) {
    if (world === this.world) return;
    this.world = world;
    for (const l of this.levels) {
      l.cache.clear();
      l.key = '';
    }
  }
  get stats() {
    return {
      levels: this.levelCount,
      nearSpacing_m: this.baseStep,
      outerHalfExtent_m: (this.cells * this.baseStep * 2 ** (this.levelCount - 1)) / 2,
      vertices: this.levels.reduce((s, l) => s + l.geo.attributes.position.count, 0),
      triangles: this.levels.reduce((s, l) => s + l.activeTriangles, 0),
      updates: this.updates,
      lastUpdateMs: this.lastUpdateMs,
      maxUpdateMs: this.maxUpdateMs,
      heightEvaluations: this.heightEvaluations,
      cachedSamples: this.levels.reduce((s, l) => s + l.cache.size, 0),
    };
  }
  dispose() {
    for (const l of this.levels) {
      this.scene.remove(l.mesh);
      l.geo.dispose();
      l.cache.clear();
    }
  }
}
OM.TerrainSystem = TerrainSystem;
OM.triangleHeight = triangle;
export { TerrainSystem, triangle };
