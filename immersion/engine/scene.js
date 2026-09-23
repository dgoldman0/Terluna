import { OM } from './om.js';
import C from './core.js';
function mesh(T, geo, mat, group, pos, scale) {
  const o = new T.Mesh(geo, mat);
  if (pos) o.position.set(...pos);
  if (scale) o.scale.set(...scale);
  o.castShadow = o.receiveShadow = true;
  group.add(o);
  return o;
}
function cylinderBetween(T, a, b, r0, r1, mat, group) {
  const p = new T.Vector3(...a),
    q = new T.Vector3(...b),
    v = q.clone().sub(p);
  const o = mesh(T, new T.CylinderGeometry(r1, r0, v.length(), 7), mat, group);
  o.position.copy(p.add(q).multiplyScalar(0.5));
  o.quaternion.setFromUnitVectors(new T.Vector3(0, 1, 0), v.normalize());
  return o;
}
function prepareMaterial(T, mat, atm, state, { sway = 0, sheltered = false, ground = false } = {}) {
  const uniforms = {
    ...atm.uniforms,
    uWet: { value: 0 },
    uCanopyWet: { value: 0 },
    uLeafWet: { value: 0 },
    uSway: { value: sway },
    uProtected: { value: sheltered ? 1 : 0 },
    uLeafSurface: { value: sway ? 1 : 0 },
  };
  mat.userData.omUniforms = uniforms;
  state.materials.push(mat);
  const coverDecl = ground ? 'varying float vOMCanopy;\n' : '';
  const extra =
    coverDecl +
    `varying vec3 vOMWorld;uniform float uWet,uCanopyWet,uLeafWet,uSway,uProtected,uLeafSurface;\n`;
  mat.onBeforeCompile = (shader) => {
    Object.assign(shader.uniforms, uniforms);
    shader.vertexShader =
      (ground ? 'attribute float omCanopy;varying float vOMCanopy;\n' : '') +
      'varying vec3 vOMWorld;uniform float uTime,uWind,uSway;\n' +
      shader.vertexShader;
    shader.vertexShader = shader.vertexShader.replace(
      '#include <begin_vertex>',
      `#include <begin_vertex>
   ${ground ? 'vOMCanopy=omCanopy;' : ''}
   vec4 omP=vec4(position,1.);\n#ifdef USE_INSTANCING\nomP=instanceMatrix*omP;\n#endif
   omP=modelMatrix*omP;
   transformed.x+=min(1.,max(0.,position.y*position.y))*uSway*sin(omP.x*.36+omP.z*.16+uTime*(.65+uWind*.035))*min(1.5,uWind*.18);
   transformed.z+=min(1.,max(0.,position.y*position.y))*uSway*.45*sin(omP.x*.18-omP.z*.28+uTime*.77)*min(1.5,uWind*.18);`,
    );
    shader.vertexShader = shader.vertexShader.replace(
      '#include <project_vertex>',
      `vec4 omW=vec4(transformed,1.);\n#ifdef USE_INSTANCING\nomW=instanceMatrix*omW;\n#endif
   vOMWorld=(modelMatrix*omW).xyz;\n#include <project_vertex>`,
    );
    shader.fragmentShader =
      extra + OM.SKY_UNIFORMS + OM.GLSL_NOISE + OM.CLOUDS + OM.CLEAR_LOOKUP + shader.fragmentShader;
    shader.fragmentShader = shader.fragmentShader.replace(
      '#include <roughnessmap_fragment>',
      `#include <roughnessmap_fragment>
   float omRoof=(1.-step(5.7,abs(vOMWorld.x-20.)))*(1.-step(4.2,abs(vOMWorld.z-33.)))*(1.-step(8.0,vOMWorld.y));
   float omWet=mix(${ground ? 'mix(uWet,uCanopyWet,vOMCanopy)' : 'mix(uWet,uLeafWet,uLeafSurface)'},0.,max(omRoof,uProtected));
   roughnessFactor=mix(roughnessFactor,max(.12,roughnessFactor*.34),omWet);
   diffuseColor.rgb*=mix(1.,.70,omWet);`,
    );
    const lightChunk = T.ShaderChunk.lights_fragment_begin;
    const needle = 'getDirectionalLightInfo( directionalLight, directLight );';
    if (!lightChunk.includes(needle))
      throw Error('Pinned Three.js directional-light shader hook has changed.');
    shader.fragmentShader = shader.fragmentShader.replace(
      '#include <lights_fragment_begin>',
      lightChunk.replace(needle, needle + '\n directLight.color *= omCloudShadow(vOMWorld);'),
    );
    shader.fragmentShader = shader.fragmentShader.replace(
      '#include <lights_fragment_end>',
      `#include <lights_fragment_end>
   reflectedLight.indirectDiffuse*=mix(1.,.24,max(omRoof,uProtected*.35))*${ground ? 'mix(1.,.58,vOMCanopy)' : '1.'};
   reflectedLight.indirectSpecular*=mix(1.,.35,max(omRoof,uProtected*.35));`,
    );
    shader.fragmentShader = shader.fragmentShader.replace(
      '#include <tonemapping_fragment>',
      `vec3 omView=vOMWorld-cameraPosition;vec3 omFogTr=exp(-uLocalExtinction*length(omView));gl_FragColor.rgb=mix(omAirColour(normalize(omView)),gl_FragColor.rgb,omFogTr);\n#ifdef TONE_MAPPING\n gl_FragColor.rgb*=uWhiteBalance;\n#endif\n#include <tonemapping_fragment>`,
    );
  };
  mat.customProgramCacheKey = () => `open-moon-material-${sway}-${sheltered}-${ground}`;
  return mat;
}
function animateDepth(T, o, uniforms, sway) {
  if (!sway) return;
  const d = new T.MeshDepthMaterial({ depthPacking: T.RGBADepthPacking, side: T.DoubleSide });
  d.onBeforeCompile = (s) => {
    Object.assign(s.uniforms, uniforms);
    s.vertexShader = 'uniform float uTime,uWind,uSway;\n' + s.vertexShader;
    s.vertexShader = s.vertexShader.replace(
      '#include <begin_vertex>',
      `#include <begin_vertex>
 vec4 p=vec4(position,1.);\n#ifdef USE_INSTANCING\np=instanceMatrix*p;\n#endif
 p=modelMatrix*p;transformed.x+=min(1.,max(0.,position.y*position.y))*uSway*sin(p.x*.36+p.z*.16+uTime*(.65+uWind*.035))*min(1.5,uWind*.18);transformed.z+=min(1.,max(0.,position.y*position.y))*uSway*.45*sin(p.x*.18-p.z*.28+uTime*.77)*min(1.5,uWind*.18);`,
    );
  };
  o.customDepthMaterial = d;
}
function createScene(T, scene, atm, world = 'moon') {
  const inheritedObjects = new Set(scene.children);
  const state = {
      materials: [],
      obstacles: [],
      trees: [],
      puddles: [],
      groundMeshes: [],
      animated: [],
      world,
    },
    obj = new T.Object3D();
  const layout = OM.world.flora();
  state.surfaceWater = new OM.SurfaceWater();
  state.fieldUniforms = OM.createFieldUniforms(T, state.surfaceWater);
  const groundMat = OM.material(T, atm, state, 'terrain');
  state.terrain = new OM.TerrainSystem(T, scene, groundMat, world);
  state.groundMeshes = state.terrain.levels.map((l) => l.mesh);
  const ecologyCounts = OM.Vegetation.populate(T, scene, atm, state, layout);
  const S = OM.world.shelter;
  let deckY = C.groundHeight(S.x, S.z, world) + 0.22;
  state.pavilionY = deckY;
  const deckMat = OM.material(T, atm, state, 'timber'),
    frameMat = new T.MeshStandardMaterial({ color: 0x313e39, metalness: 0.3, roughness: 0.5 }),
    roofMat = new T.MeshStandardMaterial({ color: 0x616960, metalness: 0.35, roughness: 0.45 });
  OM.prepareMaterial(T, frameMat, atm, state);
  OM.prepareMaterial(T, roofMat, atm, state);
  mesh(T, new T.BoxGeometry(11, 0.24, 8), deckMat, scene, [S.x, deckY - 0.12, S.z]);
  for (const dx of [-4.8, 4.8])
    for (const dz of [-3.3, 3.3]) {
      cylinderBetween(
        T,
        [S.x + dx, deckY, S.z + dz],
        [S.x + dx, deckY + 4.1, S.z + dz],
        0.115,
        0.115,
        frameMat,
        scene,
      );
      state.obstacles.push({ x: S.x + dx, z: S.z + dz, r: 0.25 });
    }
  const roof = mesh(T, new T.BoxGeometry(12, 0.16, 9.2), roofMat, scene, [S.x, deckY + 4.2, S.z]);
  roof.rotation.x = 0.045;
  const lamp = new T.PointLight(0xffd29b, 0, 22, 2);
  lamp.position.set(S.x, deckY + 3.5, S.z);
  scene.add(lamp);
  state.lamp = lamp;
  state.bulb = mesh(
    T,
    new T.CylinderGeometry(0.14, 0.14, 0.06, 16),
    new T.MeshStandardMaterial({ color: 0xffdfa7, emissive: 0xffdfa7, emissiveIntensity: 0 }),
    scene,
    [S.x, deckY + 3.92, S.z],
  );
  const ramp = new T.PlaneGeometry(4, 7, 1, 12);
  ramp.rotateX(-Math.PI / 2);
  const rp = ramp.attributes.position;
  for (let i = 0; i < rp.count; i++) {
    const x = rp.getX(i) + S.x,
      z = rp.getZ(i) + S.z - 7.5;
    rp.setXYZ(
      i,
      x,
      C.mix(C.groundHeight(x, z, world), deckY, C.smooth(0, 1, (z - S.z + 11) / 7)),
      z,
    );
  }
  ramp.computeVertexNormals();
  const rampMesh = mesh(T, ramp, deckMat, scene);
  rampMesh.userData.surfaceConforming = true;
  state.height = (x, z) => {
    let y = C.groundHeight(x, z, state.world);
    if (OM.world.roofMask(x, z)) y = Math.max(y, deckY);
    if (Math.abs(x - S.x) < 2 && z > S.z - 11 && z < S.z - 4)
      y = C.mix(y, deckY, C.smooth(0, 1, (z - S.z + 11) / 7));
    return y;
  };
  state.canMove = (x, z) => {
    if (Math.abs(x) > 130 || z > 140 || C.surfaceHeight(x, z) < 0.12) return false;
    return !state.obstacles.some((o) => Math.hypot(x - o.x, z - o.z) < o.r + 0.28);
  };
  state.update = (l, time = 0, kind = 'clear') => {
    const before = state.surfaceWater.version;
    state.surfaceWater.seek(time, kind);
    if (
      before !== state.surfaceWater.version ||
      state.waterTextureVersion !== state.surfaceWater.version
    ) {
      state.fieldUniforms.uSurfaceState.value.image.data = state.surfaceWater.textureData();
      state.fieldUniforms.uSurfaceState.value.needsUpdate = true;
      state.waterTextureVersion = state.surfaceWater.version;
    }
    if (state.ponds) state.ponds.update();
    for (const m of state.materials) {
      const u = m.userData.omUniforms;
      if (u) {
        u.uWet.value = C.clamp(l.exposed / 0.9);
        u.uCanopyWet.value = C.clamp(l.canopy / 0.7);
        u.uLeafWet.value = C.clamp(l.leaf / 0.35);
      }
    }
  };
  const decorations = scene.children.filter(
    (o) => !inheritedObjects.has(o) && !o.userData.landscape,
  );
  // What stands on the land (vegetation, stones, the shelter), for line-of-sight tests.
  state.occluders = decorations;
  const baseY = new Map(),
    treeY = state.trees.map((t) => t.y),
    baseDeckY = deckY;
  for (const o of decorations) {
    if (o.isInstancedMesh)
      baseY.set(
        o,
        Float32Array.from({ length: o.count }, (_, i) => o.instanceMatrix.array[i * 16 + 13]),
      );
    else if (o.userData.surfaceConforming)
      baseY.set(
        o,
        Float32Array.from({ length: o.geometry.attributes.position.count }, (_, i) =>
          o.geometry.attributes.position.getY(i),
        ),
      );
    else baseY.set(o, o.position.y);
  }
  state.updateTerrain = (x, z, worldKey) => {
    if (worldKey !== state.world) {
      const delta = (x, z) => C.curvatureSag(x, z, world) - C.curvatureSag(x, z, worldKey),
        m = new T.Matrix4();
      for (const o of decorations) {
        if (o.isInstancedMesh) {
          for (let i = 0; i < o.count; i++) {
            o.getMatrixAt(i, m);
            m.elements[13] = baseY.get(o)[i] + delta(m.elements[12], m.elements[14]);
            o.setMatrixAt(i, m);
          }
          o.instanceMatrix.needsUpdate = true;
          o.computeBoundingSphere();
        } else if (o.userData.surfaceConforming) {
          const p = o.geometry.attributes.position;
          for (let i = 0; i < p.count; i++)
            p.setY(i, baseY.get(o)[i] + delta(p.getX(i), p.getZ(i)));
          p.needsUpdate = true;
          o.geometry.computeVertexNormals();
          o.geometry.computeBoundingSphere();
        } else o.position.y = baseY.get(o) + delta(o.position.x, o.position.z);
      }
      state.trees.forEach((t, i) => (t.y = treeY[i] + delta(t.x, t.z)));
      deckY = baseDeckY + delta(S.x, S.z);
      state.pavilionY = deckY;
      state.world = worldKey;
    }
    state.terrain.setWorld(worldKey);
    if (state.ponds) state.ponds.update();
    return state.terrain.update(x, z);
  };
  state.debug = (mode) => {
    state.fieldUniforms.uLandscapeDebug.value = mode;
  };
  state.ponds = OM.createPonds(T, scene, atm, state);
  state.counts = {
    ...ecologyCounts,
    terrainVertices: state.terrain.stats.vertices,
    terrainLevels: state.terrain.stats.levels,
  };
  return state;
}
OM.createScene = createScene;
OM.prepareMaterial = prepareMaterial;
OM.animateDepth = animateDepth;
