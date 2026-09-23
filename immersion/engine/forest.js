/* Trees beyond a world's detailed planting, out toward the horizon.
 *
 * The world supplies its planting rule for one cell of a tree lattice
 * (OM.world.vegetation = { lattice, treeAt(i, j) }). The cells in lattice.detailed
 * are drawn as the world's full tree models; this system draws every other cell in
 * range as a camera-facing impostor (engine/impostors.js). Tiles of 16 x 16 cells
 * are generated nearest first, in a worker when the experience provides one (else
 * on the main thread under a per-frame time budget), cached, and appended to one
 * instanced draw; trees thin out by dithering toward the range limit.
 *
 * Which trees stand where belongs to the world (for the development cove, its
 * authored placeholder rule). Drawing them as cards is a perceptual convention for
 * distance, not a claim about their form.
 */
import { OM } from './om.js';
import { buildImpostorAtlas, prototypeFor, PROTOTYPES } from './impostors.js';

export const TILE_CELLS = 16;
export const FLOATS_PER_TREE = 8; // aTree (x, surface height, z, scale), aCard (prototype, tone, yaw, 0)
const now = () => (typeof performance !== 'undefined' ? performance.now() : Date.now());

/** Tile bookkeeping and instance data, independent of any renderer. */
export class ForestField {
  constructor(vegetation, surfaceHeight, { range = 3000 } = {}) {
    this.vegetation = vegetation;
    this.surfaceHeight = surfaceHeight; // uncurved landscape height; curvature is applied when drawn
    this.lattice = vegetation.lattice;
    this.range = range;
    this.cache = new Map();
    this.desired = [];
    this.pending = [];
    this.centreKey = null;
    this.generatedMs = 0;
  }
  tileSize() {
    return TILE_CELLS * this.lattice.spacing;
  }
  /** Distance from (x, z) to the nearest point of tile (ti, tj). */
  tileDistance(ti, tj, x, z) {
    const s = this.lattice.spacing,
      x0 = this.lattice.originX + (ti * TILE_CELLS - 0.5) * s,
      z0 = this.lattice.originZ + (tj * TILE_CELLS - 0.5) * s,
      size = this.tileSize();
    const dx = Math.max(x0 - x, 0, x - (x0 + size)),
      dz = Math.max(z0 - z, 0, z - (z0 + size));
    return Math.hypot(dx, dz);
  }
  /** Choose the tiles in range of (x, z), nearest first. Returns true if the set changed. */
  select(x, z) {
    const size = this.tileSize(),
      ci = Math.floor((x - this.lattice.originX) / size),
      cj = Math.floor((z - this.lattice.originZ) / size),
      key = `${ci},${cj},${this.range}`;
    if (key === this.centreKey) return false;
    this.centreKey = key;
    const n = Math.ceil(this.range / size) + 1,
      tiles = [];
    for (let tj = cj - n; tj <= cj + n; tj++)
      for (let ti = ci - n; ti <= ci + n; ti++) {
        const d = this.tileDistance(ti, tj, x, z);
        if (d <= this.range) tiles.push({ key: `${ti},${tj}`, ti, tj, d });
      }
    tiles.sort((a, b) => a.d - b.d || a.tj - b.tj || a.ti - b.ti);
    this.desired = tiles;
    this.pending = tiles.filter((t) => !this.cache.has(t.key));
    return true;
  }
  /** Instance data for one tile: every lattice cell outside the detailed planting. */
  generate(ti, tj) {
    const { i0, i1, j0, j1 } = this.lattice.detailed,
      out = [];
    for (let j = tj * TILE_CELLS; j < (tj + 1) * TILE_CELLS; j++)
      for (let i = ti * TILE_CELLS; i < (ti + 1) * TILE_CELLS; i++) {
        if (i >= i0 && i <= i1 && j >= j0 && j <= j1) continue;
        const t = this.vegetation.treeAt(i, j);
        if (!t) continue;
        const p = prototypeFor(t);
        out.push(
          t.x,
          this.surfaceHeight(t.x, t.z),
          t.z,
          t.height / PROTOTYPES[p].height,
          p,
          (t.seed % 997) / 997,
          (((t.seed >>> 10) % 3600) / 3600) * 2 * Math.PI,
          0,
        );
      }
    return new Float32Array(out);
  }
  /** Generate pending tiles until the budget (ms) is spent. Returns the tiles made. */
  work(budgetMs = Infinity) {
    const start = now(),
      made = [];
    while (this.pending.length && (made.length === 0 || now() - start < budgetMs)) {
      const t = this.pending.shift();
      if (this.cache.has(t.key)) continue;
      const t0 = now();
      this.cache.set(t.key, this.generate(t.ti, t.tj));
      this.generatedMs += now() - t0;
      made.push(t);
    }
    return made;
  }
  /** Instance data of every generated tile in the current selection, nearest first. */
  *selected() {
    for (const t of this.desired) {
      const data = this.cache.get(t.key);
      if (data) yield data;
    }
  }
  get complete() {
    return this.pending.length === 0;
  }
}

// The card faces the viewer; its image is the prototype seen from the viewer's
// direction in the tree's own (yawed) frame, dithered between the two nearest views.
const VERTEX_HEAD = /* glsl */ `
attribute vec4 aTree, aCard;
uniform float uR;
uniform vec2 uCards[${PROTOTYPES.length}];
uniform vec2 uAtlas;
varying float vOMTone, vOMBlend;
varying vec2 vOMUvB;
`;
const FRAGMENT_HEAD = /* glsl */ `
uniform float uForestFar;
varying float vOMTone, vOMBlend;
varying vec2 vOMUvB;
float omDither(vec2 p) { return fract(52.9829189 * fract(dot(p, vec2(0.06711056, 0.00583715)))); }
`;

/** The far-field trees of the active world, drawn as impostor cards. */
export class Forest {
  constructor(T, scene, atm, state, renderer, { range = 3000, worker = null } = {}) {
    this.T = T;
    // A worker takes the whole nearest-first queue ({ queue: [[ti, tj], ...]}, replaced
    // whenever the selection changes) and answers each tile with { ti, tj, tile }.
    this.worker = worker;
    this.received = false;
    if (worker) worker.onmessage = ({ data }) => this.receive(data);
    this.field = new ForestField(OM.world.vegetation, (x, z) => OM.world.landscape.height(x, z), {
      range,
    });
    this.atlas = buildImpostorAtlas(T, renderer);
    this.base = new T.PlaneGeometry(1, 1).translate(0, 0.5, 0);
    this.geometry = null;
    this.capacity = 0;
    this.count = 0;
    this.grow(65536);
    const material = new T.MeshStandardMaterial({
      map: this.atlas.colour,
      normalMap: this.atlas.normal,
      roughness: 0.85,
      metalness: 0,
    });
    OM.prepareMaterial(T, material, atm, state, {});
    const prepared = material.onBeforeCompile,
      uniforms = {
        uCards: { value: this.atlas.cards.map((c) => new T.Vector2(c.size, c.bottom)) },
        uAtlas: { value: new T.Vector2(this.atlas.views, this.atlas.rows) },
        uForestFar: { value: range },
      };
    this.uniforms = uniforms;
    this.fieldUniforms = state.fieldUniforms;
    if (this.fieldUniforms?.uForestFar) this.fieldUniforms.uForestFar.value = range;
    material.onBeforeCompile = (shader) => {
      prepared(shader);
      Object.assign(shader.uniforms, uniforms);
      shader.vertexShader = VERTEX_HEAD + shader.vertexShader;
      shader.vertexShader = shader.vertexShader
        .replace(
          '#include <beginnormal_vertex>',
          `vec3 objectNormal = vec3(omF.x, 0.0, omF.y);
   #ifdef USE_TANGENT
   vec3 objectTangent = vec3(omF.y, 0.0, -omF.x);
   #endif`,
        )
        .replace(
          '#include <begin_vertex>',
          `vec3 transformed = omBase + vec3(omF.y, 0.0, -omF.x) * position.x * omCard.x * aTree.w
     + vec3(0.0, (omCard.y + position.y * omCard.x) * aTree.w, 0.0);`,
        )
        .replace(
          '#include <uv_vertex>',
          `#include <uv_vertex>
   vec2 omCard = uCards[int(aCard.x + 0.5)];
   vec3 omBase = vec3(aTree.x, aTree.y - dot(aTree.xz, aTree.xz) / (2.0 * uR), aTree.z);
   vec2 omTo = cameraPosition.xz - omBase.xz;
   vec2 omF = omTo / max(1e-4, length(omTo));
   float omC = cos(aCard.z), omS = sin(aCard.z);
   float omView = atan(omF.x * omC - omF.y * omS, omF.x * omS + omF.y * omC) / 6.28318531 * uAtlas.x;
   omView = mod(omView + uAtlas.x, uAtlas.x);
   float omV0 = floor(omView), omV1 = mod(omV0 + 1.0, uAtlas.x);
   vec2 omCell = clamp(uv, 0.004, 0.996);
   vec2 omUv = (vec2(omV0, aCard.x) + omCell) / uAtlas;
   vOMUvB = (vec2(omV1, aCard.x) + omCell) / uAtlas;
   vOMBlend = omView - omV0;
   #ifdef USE_MAP
   vMapUv = omUv;
   #endif
   #ifdef USE_NORMALMAP
   vNormalMapUv = omUv;
   #endif
   vOMTone = aCard.y;`,
        );
      shader.fragmentShader = FRAGMENT_HEAD + shader.fragmentShader;
      // Pick this pixel's view. Texture reads take their derivatives from the
      // continuous card coordinates, so the choice never changes the mipmap level.
      shader.fragmentShader = shader.fragmentShader
        .replace(
          'void main() {',
          `void main() {
   vec2 omUv = omDither(gl_FragCoord.yx) < vOMBlend ? vOMUvB : vMapUv;
   vec2 omDx = dFdx(vMapUv), omDy = dFdy(vMapUv);`,
        )
        // The atlases hold colour and normals premultiplied by coverage.
        .replace(
          '#include <map_fragment>',
          `vec4 omTexel = textureGrad(map, omUv, omDx, omDy);
   diffuseColor *= vec4(omTexel.rgb / max(omTexel.a, 1e-3), omTexel.a);`,
        )
        .replace(
          '#include <normal_fragment_maps>',
          `vec4 omN = textureGrad(normalMap, omUv, omDx, omDy);
   vec3 mapN = omN.xyz / max(omN.a, 1e-3) * 2.0 - 1.0;
   mapN.xy *= normalScale;
   normal = normalize(tbn * mapN);`,
        );
      shader.fragmentShader = shader.fragmentShader.replace(
        '#include <alphatest_fragment>',
        // Coverage by dithered discard: distant cards sample coarse mipmaps whose alpha is
        // the crown's mean coverage, which a fixed alpha test would reject outright.
        `float omFade = 1.0 - smoothstep(uForestFar * 0.82, uForestFar, length(vOMWorld - cameraPosition));
   if (diffuseColor.a * omFade <= omDither(gl_FragCoord.xy)) discard;
   diffuseColor.a = 1.0;
   diffuseColor.rgb *= mix(0.86, 1.08, vOMTone);`,
      );
    };
    material.customProgramCacheKey = () => 'open-moon-forest-impostor';
    this.material = material;
    this.mesh = new T.Mesh(this.geometry, material);
    this.mesh.name = 'far-forest';
    this.mesh.frustumCulled = false;
    this.mesh.castShadow = false;
    this.mesh.receiveShadow = true;
    scene.add(this.mesh);
  }
  /** Enlarge the instance buffers. A new geometry is made because three.js fixes an
   * instanced geometry's maximum instance count when it is first drawn. */
  grow(capacity) {
    const T = this.T,
      next = Math.max(capacity, this.capacity * 2),
      geometry = new T.InstancedBufferGeometry(),
      old = this.geometry;
    for (const name of ['position', 'normal', 'uv'])
      geometry.setAttribute(name, this.base.getAttribute(name));
    geometry.setIndex(this.base.getIndex());
    for (const name of ['aTree', 'aCard']) {
      const array = new Float32Array(next * 4),
        previous = old?.getAttribute(name);
      if (previous) array.set(previous.array.subarray(0, this.count * 4));
      geometry.setAttribute(name, new T.InstancedBufferAttribute(array, 4));
    }
    geometry.instanceCount = this.count;
    this.geometry = geometry;
    if (this.mesh) this.mesh.geometry = geometry;
    old?.dispose();
    this.capacity = next;
  }
  append(data, upload = true) {
    const n = data.length / FLOATS_PER_TREE;
    if (!n) return;
    if (this.count + n > this.capacity) this.grow(this.count + n);
    const tree = this.geometry.getAttribute('aTree'),
      card = this.geometry.getAttribute('aCard');
    for (let k = 0; k < n; k++) {
      const o = (this.count + k) * 4,
        s = k * FLOATS_PER_TREE;
      tree.array.set(data.subarray(s, s + 4), o);
      card.array.set(data.subarray(s + 4, s + 8), o);
    }
    if (upload)
      for (const a of [tree, card]) {
        a.addUpdateRange(this.count * 4, n * 4);
        a.needsUpdate = true;
      }
    this.count += n;
    this.geometry.instanceCount = this.count;
  }
  rebuild() {
    this.count = 0;
    for (const data of this.field.selected()) this.append(data, false);
    for (const name of ['aTree', 'aCard']) {
      const a = this.geometry.getAttribute(name);
      a.clearUpdateRanges();
      a.needsUpdate = true;
    }
    this.geometry.instanceCount = this.count;
  }
  /** A tile from the worker: cache it, and draw it if it is still in range. */
  receive({ ti, tj, tile }) {
    const key = `${ti},${tj}`;
    if (this.field.cache.has(key)) return;
    this.field.cache.set(key, tile);
    this.field.pending = this.field.pending.filter((t) => t.key !== key);
    if (this.field.desired.some((t) => t.key === key)) {
      this.append(tile);
      this.received = true;
    }
  }
  /** Stream tiles around the viewer. True if the draw changed. */
  update(position, budgetMs = 3) {
    const moved = this.field.select(position.x, position.z);
    if (moved) this.rebuild();
    let changed = moved || this.received;
    this.received = false;
    if (this.worker) {
      if (moved) this.worker.postMessage({ queue: this.field.pending.map((t) => [t.ti, t.tj]) });
    } else if (this.field.pending.length) {
      const made = this.field.work(budgetMs);
      for (const t of made) this.append(this.field.cache.get(t.key));
      changed ||= made.length > 0;
    }
    return changed;
  }
  /** Generate everything in range now, on this thread (deterministic captures). */
  flush(position) {
    if (this.field.select(position.x, position.z)) this.rebuild();
    for (const t of this.field.work(Infinity)) this.append(this.field.cache.get(t.key));
  }
  setRange(range) {
    this.field.range = range;
    this.uniforms.uForestFar.value = range;
    if (this.fieldUniforms?.uForestFar) this.fieldUniforms.uForestFar.value = range;
    this.field.centreKey = null;
  }
  get stats() {
    return {
      range_m: this.field.range,
      trees: this.count,
      tilesInRange: this.field.desired.length,
      tilesPending: this.field.pending.length,
      generationMs: Math.round(this.field.generatedMs),
    };
  }
}
OM.Forest = Forest;
