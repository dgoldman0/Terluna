/* Tree impostors: images of the procedural tree models (engine/vegetation.js),
 * rendered once at start-up into an atlas for drawing distant trees as
 * camera-facing cards.
 *
 * Each prototype is drawn from eight directions around it into two multisampled
 * atlases on the GPU: albedo (the leaf and bark colours of the full models, with an
 * occlusion term that darkens the crown's interior, its underside and the shaded
 * trunk) and normals, which blend each leaf's facing with the crown's rounded form so
 * that a card shades like a crown under the scene's Sun, Earth and sky. Multisample
 * resolve makes coverage the alpha and leaves colour premultiplied by it, so mipmaps
 * average correctly; the card shader divides the alpha back out.
 *
 * The occlusion and normal blend are display choices (perceptual convention), not a
 * canopy radiative-transfer model.
 */
import { OM } from './om.js';

// Families 0-2 of the vegetation models; two mature forms and one juvenile form each.
export const PROTOTYPES = [0, 1, 2].flatMap((family) => [
  ...[0, 1].map((v) => ({
    family,
    ageClass: 'mature',
    height: family === 0 ? 7.75 : 9.95,
    seed: 7919 * (v + 1) + 104729 * family,
  })),
  { family, ageClass: 'juvenile', height: 3.65, seed: 15485863 + 104729 * family },
]);
// One row per prototype, one column per viewing direction.
export const ATLAS = { views: 8, cell: 192 };
const LEAF_COLOUR = [0.1295, 0.2122, 0.0382]; // 0x647e37 in linear sRGB, as the near leaves
const LEAF_DETAIL = 0.9; // mean of the near-leaf vein and midrib shading
const BARK_SHADE = 0.5; // the trunk stands in its crown's shade

/** The prototype index for a tree: same family and age class, variant from its seed. */
export function prototypeFor(tree) {
  let first = -1,
    count = 0;
  PROTOTYPES.forEach((p, i) => {
    if (p.family === tree.family && p.ageClass === tree.ageClass) {
      if (first < 0) first = i;
      count++;
    }
  });
  return first < 0 ? 0 : first + (tree.seed % count);
}

/** Crown centre and radii from the leaves, and the card framing (metres, prototype scale). */
export function frame(model) {
  const n = model.leaves.length,
    c = [0, 0, 0];
  for (const l of model.leaves) for (let k = 0; k < 3; k++) c[k] += l.position[k] / n;
  let rx = 0.5,
    ry = 0.5,
    reach = 0.5,
    top = 0;
  for (const l of model.leaves) {
    const [x, y, z] = l.position;
    rx = Math.max(rx, Math.hypot(x - c[0], z - c[2]));
    ry = Math.max(ry, Math.abs(y - c[1]));
    reach = Math.max(reach, Math.hypot(x, z) + 0.45 * l.size);
    top = Math.max(top, y + 0.45 * l.size);
  }
  for (const s of model.segments) {
    reach = Math.max(reach, Math.hypot(s.a[0], s.a[2]) + s.r0, Math.hypot(s.b[0], s.b[2]) + s.r1);
    top = Math.max(top, s.a[1], s.b[1]);
  }
  const bottom = -0.02 * top,
    size = 1.03 * Math.max(2 * reach, top - bottom);
  return { centre: c, radii: [rx, ry, rx], size, bottom };
}

// three.js declares instanceMatrix and instanceColor for instanced meshes.
const VERTEX = /* glsl */ `
uniform float uLeaf;
varying vec3 vWorld, vNormal, vColour;
void main() {
  vec4 p = instanceMatrix * vec4(position, 1.0);
  vWorld = p.xyz;
  vNormal = normalize(mat3(instanceMatrix) * normal);
  vColour = vec3(1.0);
  #ifdef USE_INSTANCING_COLOR
  if (uLeaf > 0.5) vColour = instanceColor;
  #endif
  gl_Position = projectionMatrix * viewMatrix * p;
}`;
const COLOUR = /* glsl */ `
uniform vec3 uAlbedo, uCentre, uRadii;
uniform float uLeaf;
varying vec3 vWorld, vNormal, vColour;
void main() {
  vec3 q = (vWorld - uCentre) / uRadii;
  float ao = mix(0.52, 1.0, smoothstep(0.2, 1.0, length(q))) * mix(0.8, 1.0, smoothstep(-1.0, 0.8, q.y));
  if (uLeaf < 0.5) ao *= ${BARK_SHADE.toFixed(2)};
  gl_FragColor = vec4(uAlbedo * vColour * ao, 1.0);
}`;
const NORMAL = /* glsl */ `
uniform vec3 uCentre, uRadii;
uniform float uLeaf;
varying vec3 vWorld, vNormal, vColour;
void main() {
  // Normals in the capture camera's frame are the card's tangent-space normals.
  mat3 view = mat3(viewMatrix);
  vec3 n = view * normalize(vNormal) * (gl_FrontFacing ? 1.0 : -1.0);
  if (uLeaf > 0.5) {
    if (n.z < 0.0) n = -n;
    vec3 crown = normalize(view * normalize((vWorld - uCentre) / (uRadii * uRadii)));
    n = normalize(mix(n, crown, 0.65));
  }
  n.z = max(n.z, 0.1);
  gl_FragColor = vec4(normalize(n) * 0.5 + 0.5, 1.0);
}`;

function prototypeScene(T, spec, bark) {
  const V = OM.Vegetation,
    model = V.skeleton({ x: 0, y: 0, z: 0, ...spec }),
    f = frame(model),
    up = new T.Vector3(0, 1, 0),
    v = new T.Vector3(),
    o = new T.Object3D(),
    colour = new T.Color();
  const material = (leaf, albedo, fragmentShader) =>
    new T.ShaderMaterial({
      uniforms: {
        uLeaf: { value: leaf },
        uAlbedo: { value: new T.Vector3(...albedo) },
        uCentre: { value: new T.Vector3(...f.centre) },
        uRadii: { value: new T.Vector3(...f.radii) },
      },
      vertexShader: VERTEX,
      fragmentShader,
      side: T.DoubleSide,
    });
  const branches = new T.InstancedMesh(
    new T.CylinderGeometry(0.8, 1, 1, 6, 1),
    null,
    model.segments.length,
  );
  model.segments.forEach((s, i) => {
    v.set(s.b[0] - s.a[0], s.b[1] - s.a[1], s.b[2] - s.a[2]);
    const len = v.length();
    o.position.set(...s.a.map((x, k) => (x + s.b[k]) * 0.5));
    o.quaternion.setFromUnitVectors(up, v.normalize());
    o.scale.set(s.r0, len, s.r0);
    o.updateMatrix();
    branches.setMatrixAt(i, o.matrix);
    branches.setColorAt(i, colour.setRGB(1, 1, 1));
  });
  const leaves = new T.InstancedMesh(V.leafGeometry(T, spec.family), null, model.leaves.length);
  model.leaves.forEach((l, i) => {
    o.position.set(...l.position);
    o.quaternion.setFromUnitVectors(up, v.set(...l.dir).normalize());
    o.rotateY(l.roll);
    o.scale.setScalar(l.size);
    o.updateMatrix();
    leaves.setMatrixAt(i, o.matrix);
    leaves.setColorAt(
      i,
      colour.setRGB(0.76 + l.tone * 0.28, 0.83 + l.tone * 0.22, 0.58 + l.tone * 0.31),
    );
  });
  const leafAlbedo = LEAF_COLOUR.map((c) => c * LEAF_DETAIL);
  const passes = {
    colour: [material(0, bark, COLOUR), material(1, leafAlbedo, COLOUR)],
    normal: [material(0, bark, NORMAL), material(1, leafAlbedo, NORMAL)],
  };
  const scene = new T.Scene();
  scene.add(branches, leaves);
  const half = f.size / 2,
    camera = new T.OrthographicCamera(-half, half, f.bottom + f.size, f.bottom, 0.1, 400);
  const use = (pass) => {
    branches.material = passes[pass][0];
    leaves.material = passes[pass][1];
  };
  /** Look at the tree from azimuth theta (the direction from the tree to the viewer). */
  const view = (theta) => {
    camera.position.set(200 * Math.sin(theta), 0, 200 * Math.cos(theta));
    camera.lookAt(0, 0, 0);
    camera.updateMatrixWorld();
  };
  const dispose = () => {
    branches.geometry.dispose();
    leaves.geometry.dispose();
    for (const list of Object.values(passes)) for (const m of list) m.dispose();
    branches.dispose();
    leaves.dispose();
  };
  return { scene, camera, use, view, frame: f, dispose };
}

/** Render every prototype from every direction into colour and normal atlases (GPU only). */
export function buildImpostorAtlas(T, renderer) {
  const { views, cell } = ATLAS,
    rows = PROTOTYPES.length,
    width = views * cell,
    height = rows * cell;
  const layer = OM.surfaceSpec?.manifest?.layers?.find((l) => l.name === 'bark'),
    bark = layer?.mean_albedo_linear || [0.098, 0.065, 0.034];
  // Half floats keep dark premultiplied albedo free of banding; bytes are the fallback.
  const half = renderer.extensions.has('EXT_color_buffer_float'),
    samples = Math.min(4, renderer.capabilities.maxSamples || 0);
  const target = () => {
    const t = new T.WebGLRenderTarget(width, height, {
      type: half ? T.HalfFloatType : T.UnsignedByteType,
      samples,
      depthBuffer: true,
      generateMipmaps: true,
      minFilter: T.LinearMipmapLinearFilter,
      magFilter: T.LinearFilter,
    });
    t.texture.colorSpace = T.NoColorSpace;
    t.texture.anisotropy = 4;
    t.scissorTest = true;
    // Allocate the full mipmap chain now; the cells below then skip regeneration.
    renderer.initRenderTarget(t);
    return t;
  };
  const targets = { colour: target(), normal: target() };
  const saved = {
    target: renderer.getRenderTarget(),
    tone: renderer.toneMapping,
    clear: renderer.getClearColor(new T.Color()),
    alpha: renderer.getClearAlpha(),
    shadows: renderer.shadowMap.enabled,
  };
  const cards = [];
  renderer.toneMapping = T.NoToneMapping;
  renderer.shadowMap.enabled = false;
  renderer.setClearColor(0x000000, 0);
  try {
    PROTOTYPES.forEach((spec, row) => {
      const p = prototypeScene(T, spec, bark),
        last = row === rows - 1;
      cards.push({ size: p.frame.size, bottom: p.frame.bottom, height: spec.height });
      for (let col = 0; col < views; col++) {
        p.view((2 * Math.PI * col) / views);
        for (const pass of ['colour', 'normal']) {
          const t = targets[pass];
          t.viewport.set(col * cell, row * cell, cell, cell);
          t.scissor.set(col * cell, row * cell, cell, cell);
          // Mipmaps are built once, after the final cell.
          t.texture.generateMipmaps = last && col === views - 1;
          p.use(pass);
          renderer.setRenderTarget(t);
          renderer.clear();
          renderer.render(p.scene, p.camera);
        }
      }
      p.dispose();
    });
  } finally {
    renderer.setRenderTarget(saved.target);
    renderer.toneMapping = saved.tone;
    renderer.setClearColor(saved.clear, saved.alpha);
    renderer.shadowMap.enabled = saved.shadows;
  }
  for (const t of Object.values(targets)) t.scissorTest = false;
  return {
    colour: targets.colour.texture,
    normal: targets.normal.texture,
    targets,
    cards,
    columns: views,
    rows,
    views,
    samples,
  };
}
