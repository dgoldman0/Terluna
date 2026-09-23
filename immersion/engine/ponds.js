/* Sub-grid pond reconstruction from conservative column storage.
 * Fine terrain spill basins pool their assigned volumes at one solved level.
 * Unresolved mobile water and above-spill excess remain explicitly accounted.
 * This reconstruction changes no SurfaceWater ledger or geological field.
 */
import { OM } from './om.js';
import C from './core.js';
import L from '../world/landscape.js';
function triangleVolume(level, sorted, area) {
  const [a, b, c] = sorted;
  if (level <= a) return 0;
  if (level >= c) return area * (level - (a + b + c) / 3);
  if (level < b) {
    const d = level - a;
    return (area * d * d * d) / (3 * (b - a) * (c - a));
  }
  const d = c - level;
  return Math.max(0, area * (level - (a + b + c) / 3 + (d * d * d) / (3 * (c - a) * (c - b))));
}
function clipTriangle(points, level) {
  let polygon = points.map((p) => [p[0], p[1], level - p[2]]),
    out = [];
  for (let i = 0; i < polygon.length; i++) {
    const a = polygon[i],
      b = polygon[(i + 1) % polygon.length],
      insideA = a[2] >= 0,
      insideB = b[2] >= 0;
    if (insideA) out.push(a);
    if (insideA !== insideB) {
      const t = a[2] / (a[2] - b[2]);
      out.push([a[0] + t * (b[0] - a[0]), a[1] + t * (b[1] - a[1]), 0]);
    }
  }
  return out;
}
class PondModel {
  constructor(water, { subdivisions = 4 } = {}) {
    if (!Number.isInteger(subdivisions) || subdivisions < 2 || subdivisions > 16)
      throw new RangeError('Pond subdivisions must be 2–16');
    this.water = water;
    this.subdivisions = subdivisions;
    this.step = water.cell / subdivisions;
    this.minX = water.minX - water.cell / 2;
    this.minZ = water.minZ - water.cell / 2;
    this.cells = water.n * subdivisions;
    this.n = this.cells + 1;
    this.version = -1;
    const n = this.n,
      step = this.step,
      heights = new Float64Array(n * n);
    for (let j = 0; j < n; j++)
      for (let i = 0; i < n; i++)
        heights[j * n + i] = L.height(this.minX + i * step, this.minZ + j * step);
    // Flood vertex elevations so spill levels and clipped surface triangles use
    // the same terrain. Diagonal routing matches the landscape's D8 convention.
    const flow = L.drainageGrid(heights, n, step),
      labels = new Int32Array(n * n).fill(-1),
      groups = [],
      queue = [];
    for (let k = 0; k < heights.length; k++) {
      if (labels[k] >= 0 || heights[k] <= 0 || flow.filled[k] - heights[k] < 0.002) continue;
      const id = groups.length,
        spill = flow.filled[k],
        group = {
          id,
          spill,
          triangles: [],
          columns: new Map(),
          min: Infinity,
          capacity: 0,
          volume: 0,
          level: spill,
          x: 0,
          z: 0,
          area: 0,
        };
      queue.length = 0;
      queue.push(k);
      labels[k] = id;
      for (let q = 0; q < queue.length; q++) {
        const index = queue[q],
          i = index % n,
          j = Math.floor(index / n);
        for (const [dx, dz] of [
          [-1, 0],
          [1, 0],
          [0, -1],
          [0, 1],
          [-1, -1],
          [1, -1],
          [-1, 1],
          [1, 1],
        ]) {
          const x = i + dx,
            z = j + dz;
          if (x < 0 || x >= n || z < 0 || z >= n) continue;
          const a = z * n + x;
          if (
            labels[a] < 0 &&
            heights[a] > 0 &&
            flow.filled[a] - heights[a] >= 0.002 &&
            Math.abs(flow.filled[a] - spill) < 1e-7
          ) {
            labels[a] = id;
            queue.push(a);
          }
        }
      }
      groups.push(group);
    }
    const area = (step * step) / 2;
    for (let j = 0; j < this.cells; j++)
      for (let i = 0; i < this.cells; i++) {
        const a = j * n + i,
          column = Math.floor(j / subdivisions) * water.n + Math.floor(i / subdivisions);
        for (const tri of [
          [a, a + n, a + 1],
          [a + 1, a + n, a + n + 1],
        ]) {
          const candidates = tri
            .filter((k) => labels[k] >= 0)
            .sort((a, b) => heights[a] - heights[b]);
          if (!candidates.length) continue;
          const group = groups[labels[candidates[0]]],
            points = tri.map((k) => [
              this.minX + (k % n) * step,
              this.minZ + Math.floor(k / n) * step,
              heights[k],
            ]),
            sorted = points.map((p) => p[2]).sort((a, b) => a - b);
          const capacity = triangleVolume(group.spill, sorted, area);
          if (capacity < 1e-10) continue;
          group.triangles.push({ points, sorted, area, column, capacity });
          group.capacity += capacity;
          group.min = Math.min(group.min, sorted[0]);
          group.columns.set(column, (group.columns.get(column) || 0) + capacity);
          group.x += (points.reduce((s, p) => s + p[0], 0) / 3) * capacity;
          group.z += (points.reduce((s, p) => s + p[1], 0) / 3) * capacity;
        }
      }
    this.groups = groups.filter((g) => g.capacity > 1e-8);
    this.columnCapacity = new Float64Array(water.n * water.n);
    for (const g of this.groups) {
      g.x /= g.capacity;
      g.z /= g.capacity;
      for (const [k, v] of g.columns) this.columnCapacity[k] += v;
    }
    this.levelBuffer = new Float32Array(n * n * 4);
    this.static = {
      basins: this.groups.length,
      spacing_m: step,
      triangles: this.groups.reduce((s, g) => s + g.triangles.length, 0),
    };
    this.update(true);
  }
  update(force = false) {
    if (!force && this.version === this.water.version) return false;
    this.version = this.water.version;
    const area = this.water.cell ** 2,
      water = this.water;
    let input = 0,
      mobile = 0,
      overflow = 0,
      represented = 0,
      active = 0,
      error = 0;
    for (let k = 0; k < water.pond.length; k++) {
      const v = (water.pond[k] * area) / 1000;
      input += v;
      if (this.columnCapacity[k] === 0) mobile += v;
    }
    for (const g of this.groups) {
      let assigned = 0;
      for (const [k, capacity] of g.columns)
        assigned += (((water.pond[k] * area) / 1000) * capacity) / this.columnCapacity[k];
      g.assigned = assigned;
      g.volume = Math.min(g.capacity, assigned);
      overflow += Math.max(0, assigned - g.capacity);
      if (g.volume < 1e-12) {
        g.level = g.min;
        g.area = 0;
        continue;
      }
      let lo = g.min,
        hi = g.spill;
      for (let i = 0; i < 38; i++) {
        const mid = (lo + hi) / 2;
        let v = 0;
        for (const t of g.triangles) v += triangleVolume(mid, t.sorted, t.area);
        if (v < g.volume) lo = mid;
        else hi = mid;
      }
      g.level = (lo + hi) / 2;
      let actual = 0;
      for (const t of g.triangles) actual += triangleVolume(g.level, t.sorted, t.area);
      error = Math.max(error, Math.abs(actual - g.volume));
      represented += actual;
      active++;
    }
    this.levelBuffer.fill(0);
    for (const g of this.groups) {
      if (g.volume < 1e-12) continue;
      for (const t of g.triangles)
        for (const p of t.points) {
          const i = Math.round((p[0] - this.minX) / this.step),
            j = Math.round((p[1] - this.minZ) / this.step),
            k = (j * this.n + i) * 4;
          this.levelBuffer[k] = g.level;
          this.levelBuffer[k + 3] = 1;
        }
    }
    this.summary = {
      ...this.static,
      activeBasins: active,
      columnPondVolume_m3: input,
      reconstructedVolume_m3: represented,
      mobileOutsideBasins_m3: mobile,
      aboveSpillVolume_m3: overflow,
      accountingResidual_m3: input - represented - mobile - overflow,
      maxBasinSolveError_m3: error,
    };
    return true;
  }
  meshData(world = 'moon') {
    const positions = [],
      depths = [],
      normals = [];
    let polygons = 0;
    for (const g of this.groups) {
      if (g.volume < 1e-12) continue;
      g.area = 0;
      for (const t of g.triangles) {
        const poly = clipTriangle(t.points, g.level);
        if (poly.length < 3) continue;
        for (let k = 1; k < poly.length - 1; k++) {
          const triangle = [poly[0], poly[k], poly[k + 1]];
          const [a, b, c] = triangle,
            area = Math.abs((b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])) / 2;
          if (area < 1e-12) continue;
          g.area += area;
          for (const p of triangle) {
            positions.push(p[0], g.level - C.curvatureSag(p[0], p[1], world), p[1]);
            depths.push(Math.max(0, p[2]));
            const R = C.worldRadius(world),
              l = Math.hypot(p[0] / R, 1, p[1] / R);
            normals.push(p[0] / R / l, 1 / l, p[1] / R / l);
          }
          polygons++;
        }
      }
    }
    return {
      positions: new Float32Array(positions),
      depths: new Float32Array(depths),
      normals: new Float32Array(normals),
      triangles: polygons,
    };
  }
  best() {
    return this.groups.reduce((best, g) => (g.volume > (best?.volume || 0) ? g : best), null);
  }
}
function createPonds(T, scene, atm, state) {
  const model = new PondModel(state.surfaceWater),
    geometry = new T.BufferGeometry();
  const uniforms = { ...atm.uniforms };
  const headTexture = new T.DataTexture(
    model.levelBuffer,
    model.n,
    model.n,
    T.RGBAFormat,
    T.FloatType,
  );
  headTexture.minFilter = headTexture.magFilter = T.NearestFilter;
  headTexture.generateMipmaps = false;
  headTexture.needsUpdate = true;
  state.fieldUniforms.uPondLevels.value = headTexture;
  state.fieldUniforms.uPondBounds.value.set(model.minX, model.minZ, model.step, model.n);
  const material = OM.makePondMaterial
    ? OM.makePondMaterial(T, atm, state)
    : new T.ShaderMaterial({
        uniforms,
        transparent: true,
        depthWrite: false,
        side: T.DoubleSide,
        toneMapped: true,
        vertexShader: `attribute float pondDepth;varying float vDepth;varying vec3 vWorld;void main(){vDepth=pondDepth;vWorld=position;gl_Position=projectionMatrix*modelViewMatrix*vec4(position,1.);}`,
        fragmentShader:
          OM.SKY_UNIFORMS +
          OM.CLEAR_LOOKUP +
          `
   varying float vDepth;varying vec3 vWorld;
   void main(){
    vec3 n=normalize(vec3(vWorld.x/uR,1.,vWorld.z/uR)),v=normalize(cameraPosition-vWorld);
    float nv=max(.03,dot(n,v)),F=.02037+.97963*pow(1.-nv,5.);
    float absorption=1.-exp(-max(0.,vDepth)*.28/nv),alpha=F+(1.-F)*absorption;
    vec3 reflected=omClear(reflect(-v,n)),scatter=uDiffuse/OM_PI*vec3(.011,.039,.042);
    vec3 colour=(reflected*F+scatter*(1.-F)*absorption)/max(.0001,alpha);
    vec3 h=normalize(v+uSun);float nh=max(0.,dot(n,h));colour+=uDirect*pow(nh,4000.)*.12;
    vec3 tr=exp(-uLocalExtinction*length(cameraPosition-vWorld));colour=mix(omAirColour(normalize(vWorld-cameraPosition)),colour,tr);
    gl_FragColor=vec4(colour,alpha*smoothstep(0.,.0015,vDepth));
    #ifdef TONE_MAPPING
     gl_FragColor.rgb*=uWhiteBalance;
    #endif
    #include <tonemapping_fragment>
    #include <colorspace_fragment>
   }`,
      });
  const mesh = new T.Mesh(geometry, material);
  mesh.name = 'volume-solved-ponds';
  mesh.userData.ponds = true;
  mesh.renderOrder = 2;
  scene.add(mesh);
  let world = state.world,
    capacity = 0;
  function update(force = false) {
    const changed = model.update(force);
    if (!changed && world === state.world && !force) return;
    world = state.world;
    headTexture.image.data = model.levelBuffer;
    headTexture.needsUpdate = true;
    const data = model.meshData(world),
      count = data.depths.length;
    if (count > capacity) {
      capacity = Math.max(count, Math.ceil(capacity * 1.5), 512);
      geometry.setAttribute(
        'position',
        new T.BufferAttribute(new Float32Array(capacity * 3), 3).setUsage(T.DynamicDrawUsage),
      );
      geometry.setAttribute(
        'normal',
        new T.BufferAttribute(new Float32Array(capacity * 3), 3).setUsage(T.DynamicDrawUsage),
      );
      geometry.setAttribute(
        'pondDepth',
        new T.BufferAttribute(new Float32Array(capacity), 1).setUsage(T.DynamicDrawUsage),
      );
    }
    if (count) {
      for (const [name, values] of [
        ['position', data.positions],
        ['normal', data.normals],
        ['pondDepth', data.depths],
      ]) {
        const a = geometry.attributes[name];
        a.array.set(values);
        a.clearUpdateRanges();
        a.addUpdateRange(0, values.length);
        a.needsUpdate = true;
      }
      geometry.computeBoundingSphere();
    }
    geometry.setDrawRange(0, count);
    mesh.visible = count > 0;
    model.summary.renderedTriangles = data.triangles;
  }
  update(true);
  return {
    model,
    mesh,
    update,
    dispose() {
      scene.remove(mesh);
      geometry.dispose();
      material.dispose();
      headTexture.dispose();
    },
  };
}
OM.PondModel = PondModel;
OM.createPonds = createPonds;
export { PondModel, triangleVolume, clipTriangle };
