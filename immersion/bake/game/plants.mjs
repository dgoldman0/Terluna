/* Export the web experience's plant models for the Unreal game: tree templates grown by
 * engine/vegetation.js's skeleton (trunk, limbs and leaves) and its grass tufts, as
 * indexed triangle meshes. World placement stays in trees.csv; the game picks a template
 * by family, age class and seed and scales it to each tree's height.
 *
 * The models are the web experience's authored placeholders, not botany. A template is
 * one tree of the template height; scaling it to other heights also scales its leaves
 * (the web grows each tree at its own height, with leaves of fixed size).
 *
 * Mesh files (little-endian): "TMSH", uint32 version 1, vertex count, index count and
 * section count; per section uint32 first index and index count; then per vertex float32
 * position xyz (Unreal centimetres: X = x, Y = z, Z = y, the base at the origin), float32
 * normal xyz, float32 uv (bark: metres around and along the limb; leaves and blades:
 * across and along, 0 to 1) and uint8 RGBA colour (linear tone, 255 = 1); then uint32
 * indices, three per triangle, front faces as in three.js (the axis swap keeps them
 * front-facing in Unreal).
 */
import fs from 'node:fs';
import path from 'node:path';
import crypto from 'node:crypto';
import * as THREE from 'three';
import Vegetation from '../../engine/vegetation.js';

// Mid-range heights of the planting rule's classes (world/cove-flora.js treeAt): mature
// 5.5-10 m for family 0 and 7.7-12.2 m for families 1 and 2, juvenile 2.1-5.2 m.
const TEMPLATE_HEIGHT = { mature: [7.75, 9.95, 9.95], juvenile: [3.65, 3.65, 3.65] };
const VARIANTS = [18061, 42013, 77317, 90491]; // skeleton seeds, one tree shape each
// The web's colours (three.js sRGB hex, linear once converted): leaves 0x647e37, grass
// 0x85824a, 0x60753c and 0x688442 by community.
const LEAF_COLOUR = 0x647e37,
  GRASS_COLOURS = [0x85824a, 0x60753c, 0x688442];

class MeshBuilder {
  constructor() {
    this.p = [];
    this.n = [];
    this.uv = [];
    this.c = [];
    this.index = [];
    this.sections = [];
  }
  vertex([x, y, z], [nx, ny, nz], [u, v], [r, g, b]) {
    // three.js metres, y up -> Unreal centimetres, Z up.
    this.p.push(x * 100, z * 100, y * 100);
    this.n.push(nx, nz, ny);
    this.uv.push(u, v);
    this.c.push(...[r, g, b].map((k) => Math.round(255 * Math.min(1, Math.max(0, k)))), 255);
    return this.p.length / 3 - 1;
  }
  section(fill) {
    const first = this.index.length;
    fill();
    this.sections.push([first, this.index.length - first]);
  }
  write(file) {
    const nv = this.p.length / 3,
      ni = this.index.length,
      ns = this.sections.length,
      head = Buffer.alloc(20 + ns * 8);
    head.write('TMSH', 0, 'ascii');
    head.writeUInt32LE(1, 4);
    head.writeUInt32LE(nv, 8);
    head.writeUInt32LE(ni, 12);
    head.writeUInt32LE(ns, 16);
    this.sections.forEach(([first, count], k) => {
      head.writeUInt32LE(first, 20 + k * 8);
      head.writeUInt32LE(count, 24 + k * 8);
    });
    const floats = (a) => Buffer.from(new Float32Array(a).buffer);
    fs.writeFileSync(
      file,
      Buffer.concat([
        head,
        floats(this.p),
        floats(this.n),
        floats(this.uv),
        Buffer.from(Uint8Array.from(this.c)),
        Buffer.from(new Uint32Array(this.index).buffer),
      ]),
    );
    return { vertices: nv, triangles: ni / 3 };
  }
}

const linear = (hex) => new THREE.Color(hex).toArray(); // three.js converts sRGB hex to linear

/** An 8-sided tube from a (radius r0) to b (radius r1), UVs in metres. */
function limb(m, { a, b, r0, r1 }) {
  const axis = new THREE.Vector3(b[0] - a[0], b[1] - a[1], b[2] - a[2]),
    length = axis.length();
  axis.normalize();
  const side = Math.abs(axis.y) < 0.9 ? new THREE.Vector3(0, 1, 0) : new THREE.Vector3(1, 0, 0),
    u = new THREE.Vector3().crossVectors(axis, side).normalize(),
    w = new THREE.Vector3().crossVectors(axis, u);
  const rings = [];
  for (const [end, radius, along] of [
    [a, r0, 0],
    [b, r1, length],
  ]) {
    const ring = [];
    for (let k = 0; k <= 8; k++) {
      const t = (k / 8) * Math.PI * 2,
        d = u.clone().multiplyScalar(Math.cos(t)).addScaledVector(w, Math.sin(t));
      ring.push(
        m.vertex(
          [end[0] + d.x * radius, end[1] + d.y * radius, end[2] + d.z * radius],
          [d.x, d.y, d.z],
          [(k / 8) * Math.PI * 2 * r0, along],
          [1, 1, 1],
        ),
      );
    }
    rings.push(ring);
  }
  for (let k = 0; k < 8; k++) {
    const [i0, i1, j0, j1] = [rings[0][k], rings[0][k + 1], rings[1][k], rings[1][k + 1]];
    m.index.push(i0, i1, j0, i1, j1, j0); // counter-clockwise from outside: front faces out
  }
}

/** Copies of a three.js geometry, each placed by a matrix and tinted. */
function place(m, geometry, placements) {
  const pos = geometry.attributes.position,
    nor = geometry.attributes.normal,
    uv = geometry.attributes.uv,
    idx = geometry.index.array,
    p = new THREE.Vector3(),
    n = new THREE.Vector3(),
    normalMatrix = new THREE.Matrix3();
  for (const { matrix, tone } of placements) {
    normalMatrix.getNormalMatrix(matrix);
    const base = [];
    for (let i = 0; i < pos.count; i++) {
      p.fromBufferAttribute(pos, i).applyMatrix4(matrix);
      n.fromBufferAttribute(nor, i).applyMatrix3(normalMatrix).normalize();
      base.push(m.vertex([p.x, p.y, p.z], [n.x, n.y, n.z], [uv.getX(i), uv.getY(i)], tone));
    }
    for (const i of idx) m.index.push(base[i]);
  }
}

/** One tree template: the skeleton of a tree of this family, age and seed at the origin. */
function treeModel(family, ageClass, seed) {
  const height = TEMPLATE_HEIGHT[ageClass][family],
    { segments, leaves } = Vegetation.skeleton({ x: 0, y: 0, z: 0, height, family, ageClass, seed }),
    m = new MeshBuilder(),
    up = new THREE.Vector3(0, 1, 0),
    q = new THREE.Quaternion(),
    o = new THREE.Object3D();
  // Section 0: bark (the web draws limbs thinner than 3 mm as nothing).
  m.section(() => segments.filter((s) => s.r0 >= 0.003).forEach((s) => limb(m, s)));
  // Section 1: leaves, as the web places and tints them.
  m.section(() =>
    place(
      m,
      Vegetation.leafGeometry(THREE, family),
      leaves.map((l) => {
        o.position.set(...l.position);
        q.setFromUnitVectors(up, new THREE.Vector3(...l.dir).normalize());
        o.quaternion.copy(q);
        o.rotateY(l.roll);
        o.scale.setScalar(l.size);
        o.updateMatrix();
        return {
          matrix: o.matrix.clone(),
          tone: [0.76 + l.tone * 0.28, 0.83 + l.tone * 0.22, 0.58 + l.tone * 0.31],
        };
      }),
    ),
  );
  return { m, height };
}

/** A grass tuft of one community (0 grass, 1 woodland, 2 wet margin), 1 m tall, untinted. */
function tuftModel(family) {
  const m = new MeshBuilder();
  m.section(() => place(m, Vegetation.tuftGeometry(THREE, family), [{ matrix: new THREE.Matrix4(), tone: [1, 1, 1] }]));
  return m;
}

/** Write the plant models to <dir>/plants and return their manifest. */
export function writePlants(dir, immersion) {
  const out = path.join(dir, 'plants');
  fs.mkdirSync(out, { recursive: true });
  const trees = [];
  for (const ageClass of ['mature', 'juvenile'])
    for (let family = 0; family < 3; family++)
      VARIANTS.forEach((seed, variant) => {
        const { m, height } = treeModel(family, ageClass, seed),
          file = `tree_${family}_${ageClass}_${variant}.tmsh`,
          counts = m.write(path.join(out, file));
        trees.push({ file, family, age_class: ageClass, variant, seed, height_m: height, ...counts });
      });
  const tufts = [0, 1, 2].map((family) => {
    const file = `tuft_${family}.tmsh`,
      counts = tuftModel(family).write(path.join(out, file));
    return { file, community: family, colour_linear: linear(GRASS_COLOURS[family]), ...counts };
  });
  const sha = (f) => crypto.createHash('sha256').update(fs.readFileSync(f)).digest('hex');
  const manifest = {
    schema: 'terluna.game.plants/1',
    evidence:
      "The web experience's placeholder plant models (engine/vegetation.js), not botany: " +
      'authored shapes for building the engine.',
    format:
      'TMSH meshes, described in immersion/bake/game/plants.mjs; Unreal centimetres, base at the origin.',
    trees: {
      templates: trees,
      sections: ['bark', 'leaves'],
      choice:
        'By family, age class and seed % variants (trees.csv); scaled by height_m / template height_m, ' +
        'turned by a yaw drawn from the seed.',
      bark: { surface_layer: 'bark', tile_m: 1.1 },
      leaves: { colour_linear: linear(LEAF_COLOUR), note: 'Times the vertex tone; two-sided.' },
    },
    tufts: {
      templates: tufts,
      note: 'One metre tall; the web scales each clump by its own height and tints it.',
    },
    sources: ['engine/vegetation.js', 'bake/game/plants.mjs'].map((f) => ({
      file: `immersion/${f}`,
      sha256: sha(path.join(immersion, f)),
    })),
  };
  fs.writeFileSync(path.join(out, 'manifest.json'), JSON.stringify(manifest, null, 1));
  return manifest;
}
