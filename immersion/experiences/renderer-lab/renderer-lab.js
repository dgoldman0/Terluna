/* Optional r186 TSL renderer using the production landscape and ecology model.
 * No automatic promotion: actual backend, omissions and synchronization mode
 * are part of every exported record. WebGL reference remains separately usable.
 */
import { OM } from '../../engine/om.js';
import C from '../../engine/core.js';
import { LEGACY_MOON_GRAVITY } from '../../../shared/constants.js';
const $ = (id) => document.getElementById(id);
OM.bootLab = async function (modules) {
  const T = { ...modules.webgl, ...modules.webgpu },
    N = modules.tsl;
  if (!N || !T.WebGPURenderer) throw Error('The complete r186 node dependency graph is required.');
  const requested =
    globalThis.OM_LAB_BACKEND || new URLSearchParams(location.search).get('backend') || 'auto';
  let adapter = null,
    device = null,
    adapterInfo = null,
    reason = null;
  if (requested !== 'webgl') {
    if (!navigator.gpu)
      reason =
        'WebGPU is unavailable in this browser/context. Native WebGPU requires a secure context and a compatible adapter.';
    else
      try {
        adapter = await navigator.gpu.requestAdapter({ powerPreference: 'high-performance' });
        if (!adapter) throw Error('No WebGPU adapter returned.');
        adapterInfo = {
          vendor: adapter.info?.vendor,
          architecture: adapter.info?.architecture,
          device: adapter.info?.device,
          description: adapter.info?.description,
          isFallbackAdapter: adapter.info?.isFallbackAdapter ?? adapter.isFallbackAdapter ?? null,
        };
        if (!adapter.features.has('float32-filterable'))
          throw Error(
            'This r186 experiment requires float32-filterable field textures. Use its TSL/WebGL2 mode on this adapter.',
          );
        const features = [
          'timestamp-query',
          'float32-filterable',
          'shader-f16',
          'core-features-and-limits',
        ].filter((f) => adapter.features.has(f));
        device = await adapter.requestDevice({ requiredFeatures: features });
      } catch (e) {
        reason = e.message;
      }
  }
  if (requested === 'webgpu' && !device)
    throw Error(reason || 'Native WebGPU could not be initialized. Select the TSL/WebGL2 mode.');
  const renderer = new T.WebGPURenderer({
    canvas: $('world'),
    forceWebGL: !device,
    device: device || undefined,
    antialias: false,
    alpha: false,
  });
  await renderer.init();
  renderer.info.autoReset = false;
  renderer.setPixelRatio(1);
  renderer.setSize(innerWidth, innerHeight);
  renderer.outputColorSpace = T.SRGBColorSpace;
  renderer.toneMapping = T.ACESFilmicToneMapping;
  renderer.toneMappingExposure = 0.244647046;
  renderer.shadowMap.enabled = true;
  renderer.shadowMap.type = T.PCFShadowMap;
  const backend = renderer.backend.isWebGPUBackend ? 'webgpu' : 'webgl2',
    gl = renderer.backend.gl;
  if (backend !== 'webgpu' && gl) {
    const ext = gl.getExtension('WEBGL_debug_renderer_info');
    adapterInfo = {
      renderer: ext ? gl.getParameter(ext.UNMASKED_RENDERER_WEBGL) : gl.getParameter(gl.RENDERER),
    };
  }
  const software =
    /swiftshader|llvmpipe|softpipe|software/i.test(JSON.stringify(adapterInfo)) ||
    adapterInfo?.isFallbackAdapter === true;
  $('backend-choice').value = requested;
  $('backend-choice').onchange = () => {
    const u = new URL(location.href);
    u.searchParams.set('backend', $('backend-choice').value);
    location.href = u.href;
  };
  const report = {
    schema: 'open-moon-renderer-lab/1',
    revision: T.REVISION,
    requested,
    backend,
    nativeAPIPresent: !!navigator.gpu,
    adapter: adapterInfo,
    software,
    secureContext: isSecureContext,
    fallbackReason: reason,
    hardwareBenchmark: false,
    featuresImplemented: {
      sharedLandscape: true,
      sharedEcology: true,
      sharedWaterStores: true,
      sharedPondGeometry: true,
      terrainGeometry: true,
      materialMaps: true,
      clearSkyLUT: true,
      localAerialPerspective: true,
      directShadows: true,
      depthBasedSeaAbsorption: true,
      planarReflection: true,
      cloudTransport: false,
      animatedRain: false,
      exactLeafLighting: false,
      leafSway: false,
    },
    runs: [],
  };
  $('backend').textContent =
    `Three.js r186 · ${backend === 'webgpu' ? 'native WebGPU' : 'TSL / WebGL2'}${software ? ' · software adapter' : ''}`;
  const scene = new T.Scene(),
    camera = new T.PerspectiveCamera(58, innerWidth / innerHeight, 0.12, 24000);
  camera.rotation.order = 'YXZ';
  const atm = new OM.Atmosphere(T, OM.skyData);
  await atm.init(null, true);
  await OM.loadSurfaceAssets(T);
  const node = OM.createNodeSystem(T, N, atm, renderer);
  OM.material = node.factory;
  OM.prepareMaterial = (T, mat, a, state) => {
    state.materials.push(mat);
    return mat;
  };
  OM.animateDepth = () => {};
  OM.makePondMaterial = node.pondFactory;
  const geography = OM.createScene(T, scene, atm, 'moon');
  for (const obj of scene.children) {
    if (obj.geometry?.attributes.position && !obj.geometry.attributes.uv) {
      const p = obj.geometry.attributes.position,
        uvs = new Float32Array(p.count * 2);
      for (let i = 0; i < p.count; i++) {
        uvs[i * 2] = p.getX(i);
        uvs[i * 2 + 1] = p.getZ(i);
      }
      obj.geometry.setAttribute('uv', new T.BufferAttribute(uvs, 2));
    }
    if (obj.isInstancedMesh && obj.material.userData.nodeKind === 'wood') {
      const a = new Float32Array(obj.count * 2),
        m = new T.Matrix4(),
        x = new T.Vector3(),
        y = new T.Vector3();
      for (let i = 0; i < obj.count; i++) {
        obj.getMatrixAt(i, m);
        x.setFromMatrixColumn(m, 0);
        y.setFromMatrixColumn(m, 1);
        a[i * 2] = x.length() * C.TAU;
        a[i * 2 + 1] = y.length();
      }
      obj.geometry = obj.geometry.clone();
      obj.geometry.setAttribute('omBarkScale', new T.InstancedBufferAttribute(a, 2));
    }
    if (obj.material && !obj.material.isNodeMaterial) {
      const old = obj.material,
        mat = new T.MeshStandardNodeMaterial({
          color: old.color,
          roughness: old.roughness ?? 0.8,
          metalness: old.metalness ?? 0,
          emissive: old.emissive,
          emissiveIntensity: old.emissiveIntensity ?? 0,
        });
      mat.outputNode = N.vec4(node.display(node.fog(N.output.rgb, N.positionWorld)), N.output.a);
      mat.maskNode = node.U.mirror.lessThan(0.5).or(N.positionWorld.y.greaterThan(0.015));
      obj.material = mat;
    }
  }
  const sun = new T.DirectionalLight(0xffffff, 1);
  sun.castShadow = true;
  sun.shadow.mapSize.set(1024, 1024);
  Object.assign(sun.shadow.camera, {
    left: -64,
    right: 64,
    top: 64,
    bottom: -64,
    near: 1,
    far: 900,
  });
  sun.shadow.bias = -0.0001;
  sun.shadow.normalBias = 0.007;
  sun.shadow.radius = 2;
  scene.add(sun, sun.target);
  const s = { world: 'moon', phase: 0, time: 0, yaw: 0, pitch: -0.04 };
  camera.position.set(0, geography.height(0, 9) + 1.7, 9);
  camera.rotation.set(s.pitch, s.yaw, 0);
  scene.backgroundNode = node.sky;
  const envScene = new T.Scene();
  envScene.backgroundNode = node.clear(N.positionWorldDirection);
  const pmrem = new T.PMREMGenerator(renderer);
  let envRT = null;
  const ocean = createNodeWater(T, N, renderer, scene, camera, atm, node, geography);
  function sync(environment = false) {
    const w = C.weatherAt(0, 'clear');
    geography.updateTerrain(camera.position.x, camera.position.z, s.world);
    const solar = atm.update(s.world, s.phase, w, 0, camera.position, 0);
    node.sync();
    node.refresh(geography);
    const max = Math.max(...atm.directRaw, 1e-15);
    sun.color.setRGB(...atm.directRaw.map((c) => Math.max(0, c) / max));
    sun.intensity = solar.y > -0.0047 ? max / 8500 : 0;
    sun.position.set(
      camera.position.x + solar.x * 350,
      camera.position.y + solar.y * 350,
      camera.position.z,
    );
    sun.target.position.copy(camera.position);
    sun.target.updateMatrixWorld();
    if (environment || !envRT) {
      const previousTM = renderer.toneMapping,
        previousDisplay = node.U.display.value;
      renderer.toneMapping = T.NoToneMapping;
      node.U.display.value = 0;
      const next = pmrem.fromScene(envScene, 0, 0.1, 10, { size: 256 });
      if (envRT) envRT.dispose();
      envRT = next;
      scene.environment = next.texture;
      renderer.toneMapping = previousTM;
      node.U.display.value = previousDisplay;
    }
  }
  function render() {
    renderer.info.reset();
    sync();
    camera.rotation.set(s.pitch, s.yaw, 0);
    camera.updateMatrixWorld();
    ocean.render();
    renderer.render(scene, camera);
  }
  const completionPixel = new Uint8Array(4);
  async function complete() {
    if (backend === 'webgpu') await renderer.backend.device.queue.onSubmittedWorkDone();
    else {
      const g = renderer.backend.gl;
      g.readPixels(1, 1, 1, 1, g.RGBA, g.UNSIGNED_BYTE, completionPixel);
      g.finish();
      if (g.isContextLost()) throw Error('TSL WebGL2 context lost.');
    }
  }
  function snapshot() {
    return {
      world: s.world,
      phase: s.phase,
      weather: 'clear',
      waterReplaySeconds: s.time,
      position: camera.position.toArray(),
      yaw: s.yaw,
      pitch: s.pitch,
      fov: camera.fov,
      exposure: renderer.toneMappingExposure,
      counts: geography.counts,
      terrain: geography.terrain.stats,
      water: geography.surfaceWater.ledger,
      ponds: geography.ponds.model.summary,
      renderPixels: [innerWidth, innerHeight],
      samples: 0,
      shadowMap: 1024,
      reflection: [ocean.reflectionTarget.width, ocean.reflectionTarget.height],
    };
  }
  async function setWater(t) {
    s.time = t;
    geography.update(C.ledgerAt(t, 'episode'), t, 'episode');
    node.refresh(geography);
    render();
    await complete();
    status();
  }
  function status() {
    const r = geography.ponds.model.summary;
    $('readout').textContent =
      `${geography.counts.trees} trees · ${geography.counts.leaves.toLocaleString()} leaf instances · ${r.activeBasins} wet basins · ${r.reconstructedVolume_m3.toFixed(2)} m³ reconstructed ponding`;
  }
  $('dry').onclick = () => setWater(0);
  $('wet').onclick = () => setWater(14400);
  $('shore').onclick = () => {
    camera.position.set(0, geography.height(0, 9) + 1.7, 9);
    s.yaw = 0;
    s.pitch = -0.04;
    render();
  };
  $('woods').onclick = () => {
    const p = C.WAYPOINTS.find((p) => p.id === 'path');
    camera.position.set(p.x, geography.height(p.x, p.z) + 1.7, p.z);
    s.yaw = p.yaw;
    s.pitch = p.pitch;
    render();
  };
  $('pools').onclick = () => {
    const p = geography.ponds.model.best();
    if (!p) return;
    camera.position.set(p.x + 4, geography.height(p.x + 4, p.z + 7) + 1.7, p.z + 7);
    s.yaw = Math.atan2(4, 7);
    s.pitch = -0.35;
    render();
  };
  $('phase').oninput = () => {
    s.phase = +$('phase').value;
    sync(true);
    render();
  };
  let drag = null;
  $('world').onpointerdown = (e) => {
    drag = [e.clientX, e.clientY];
    $('world').setPointerCapture(e.pointerId);
  };
  $('world').onpointermove = (e) => {
    if (!drag) return;
    s.yaw -= (e.clientX - drag[0]) * 0.003;
    s.pitch = C.clamp(s.pitch - (e.clientY - drag[1]) * 0.003, -1.5, 1.5);
    drag = [e.clientX, e.clientY];
    render();
  };
  $('world').onpointerup = () => (drag = null);
  addEventListener('keydown', (e) => {
    if (!['KeyW', 'KeyA', 'KeyS', 'KeyD'].includes(e.code)) return;
    const f = e.code === 'KeyW' ? 1 : e.code === 'KeyS' ? -1 : 0,
      side = e.code === 'KeyD' ? 1 : e.code === 'KeyA' ? -1 : 0,
      x = camera.position.x + (-Math.sin(s.yaw) * f + Math.cos(s.yaw) * side) * 0.5,
      z = camera.position.z + (-Math.cos(s.yaw) * f - Math.sin(s.yaw) * side) * 0.5;
    if (geography.canMove(x, z)) camera.position.set(x, geography.height(x, z) + 1.7, z);
    render();
  });
  async function benchmark(frames = 24) {
    if (!Number.isInteger(frames) || frames < 1 || frames > 120)
      throw Error('Expected 1–120 frames.');
    $('bench').disabled = true;
    const before = camera.position.clone(),
      yaw = s.yaw,
      pitch = s.pitch;
    const times = [],
      cpu = [],
      draws = [];
    let startState;
    try {
      camera.position.set(0, geography.height(0, 9) + 1.7, 9);
      s.yaw = 0;
      s.pitch = -0.04;
      for (let i = 0; i < 3; i++) {
        render();
        await complete();
      }
      startState = snapshot();
      for (let i = 0; i < frames; i++) {
        const x = i * 0.15,
          z = 9 + i * 0.025;
        camera.position.set(x, geography.height(x, z) + 1.7, z);
        const start = performance.now();
        render();
        cpu.push(performance.now() - start);
        await complete();
        times.push(performance.now() - start);
        draws.push({ ...renderer.info.render });
        await new Promise((r) => setTimeout(r, 0));
      }
      const quant = (a, q) =>
        [...a].sort((a, b) => a - b)[Math.min(a.length - 1, Math.floor(a.length * q))];
      const run = {
        frames,
        method:
          'Synchronized frame latency: WebGL2 uses synchronous pixel readback plus finish; WebGPU uses queue.onSubmittedWorkDone. Includes CPU and all scene passes; excludes sleep',
        state: startState,
        median_ms: quant(times, 0.5),
        p95_ms: quant(times, 0.95),
        max_ms: Math.max(...times),
        submissionMedian_ms: quant(cpu, 0.5),
        times_ms: times,
        draws,
        memory: { ...renderer.info.memory },
        completionPixel: Array.from(completionPixel),
      };
      report.runs.push(run);
      report.hardwareBenchmark = !software;
      $('result').textContent =
        `${backend}: median ${run.median_ms.toFixed(1)} ms · p95 ${run.p95_ms.toFixed(1)} ms${software ? ' (software adapter; no consumer-GPU claim)' : ''}`;
      return run;
    } finally {
      camera.position.copy(before);
      s.yaw = yaw;
      s.pitch = pitch;
      $('bench').disabled = false;
    }
  }
  $('bench').onclick = () =>
    benchmark().catch((e) => {
      $('result').textContent = e.message;
    });
  $('export').onclick = () =>
    OM.downloadBlob(
      new Blob([JSON.stringify({ ...report, state: snapshot() }, null, 2)], {
        type: 'application/json',
      }),
      'Open_Moon_Renderer_Lab.json',
    );
  addEventListener('resize', () => {
    camera.aspect = innerWidth / innerHeight;
    camera.updateProjectionMatrix();
    renderer.setSize(innerWidth, innerHeight);
    render();
  });
  const started = performance.now();
  sync(true);
  render();
  await complete();
  report.firstRender_ms = performance.now() - started;
  report.state = snapshot();
  status();
  $('boot').hidden = true;
  globalThis.rendererLab = {
    ready: true,
    renderer,
    scene,
    camera,
    geography,
    atm,
    node,
    ocean,
    state: s,
    report,
    snapshot,
    render,
    setWater,
    benchmark,
    complete,
  };
};
function createNodeWater(T, N, renderer, scene, camera, atm, system, geography) {
  const {
      Fn,
      vec2,
      vec3,
      vec4,
      float,
      uniform,
      texture,
      positionWorld,
      positionGeometry,
      positionView,
      cameraPosition,
      cameraViewMatrix,
      screenUV,
      If,
    } = N,
    U = system.U;
  const opaque = new T.RenderTarget(1, 1, { type: T.HalfFloatType, depthBuffer: true });
  opaque.depthTexture = new T.DepthTexture(1, 1, T.UnsignedIntType);
  const mirror = new T.RenderTarget(1, 1, {
    type: T.HalfFloatType,
    depthBuffer: true,
    generateMipmaps: true,
    minFilter: T.LinearMipmapLinearFilter,
  });
  const opaqueNode = texture(opaque.texture),
    depthNode = texture(opaque.depthTexture),
    mirrorNode = texture(mirror.texture),
    mirrorMatrix = uniform(new T.Matrix4()),
    near = uniform(camera.near),
    far = uniform(camera.far);
  const waves = Fn(([p, micro]) => {
    const sum = vec3(0).toVar();
    for (let i = 0; i < C.WATER_BANDS.length; i++) {
      const k = C.WATER_BANDS[i],
        theta = 0.95 + Math.sin(i * 2.399) * 0.71,
        dir = vec2(Math.cos(theta), Math.sin(theta));
      const a = U.uWind
          .mul(0.18)
          .add(0.4)
          .mul(0.095 * Math.pow(0.095 / k, 1.12)),
        omega = float(Math.sqrt(LEGACY_MOON_GRAVITY * k * Math.tanh(k * 12)));
      const q = p
          .dot(dir)
          .mul(k)
          .sub(U.uTime.mul(omega))
          .add(i * 2.721),
        weight = i < 6 ? float(1) : micro;
      sum.addAssign(
        vec3(q.sin(), dir.x.mul(q.cos()).mul(-k), dir.y.mul(q.cos()).mul(-k)).mul(a).mul(weight),
      );
    }
    const env = system.envelope(p, U.uWind, geography),
      e = 0.001,
      dx = system
        .envelope(p.add(vec2(e, 0)), U.uWind, geography)
        .sub(system.envelope(p.sub(vec2(e, 0)), U.uWind, geography))
        .div(2 * e),
      dz = system
        .envelope(p.add(vec2(0, e)), U.uWind, geography)
        .sub(system.envelope(p.sub(vec2(0, e)), U.uWind, geography))
        .div(2 * e);
    return vec3(
      sum.x.mul(env),
      sum.y.mul(env).sub(sum.x.mul(dx)),
      sum.z.mul(env).sub(sum.x.mul(dz)),
    );
  });
  const mat = new T.NodeMaterial();
  mat.positionNode = Fn(() => {
    const p = positionGeometry,
      w = waves(p.xz, float(0)),
      r2 = p.xz.dot(p.xz),
      sag = r2.div(U.uR.add(U.uR.mul(U.uR).sub(r2).max(1).sqrt()));
    return vec3(p.x, w.x.sub(sag), p.z);
  })();
  mat.fragmentNode = Fn(() => {
    const p = positionWorld,
      w = waves(p.xz, float(1)),
      v = cameraPosition.sub(p).normalize();
    const footprint = N.max(N.dFdx(p.xz).length(), N.dFdy(p.xz).length()),
      ripple = float(1).sub(footprint.smoothstep(0.08, 0.6));
    const fine = vec2(
      p.x.mul(38).add(p.z.mul(17)).add(U.uTime.mul(0.7)).sin(),
      p.z.mul(41).sub(p.x.mul(13)).add(U.uTime.mul(0.5)).sin(),
    )
      .mul(0.015)
      .mul(ripple);
    const n = vec3(
        w.y.add(fine.x).add(p.x.div(U.uR)),
        1,
        w.z.add(fine.y).add(p.z.div(U.uR)),
      ).normalize(),
      nv = n.dot(v).max(0.001),
      F = float(1).sub(nv).pow(5).mul(0.9796).add(0.0204);
    const uv = screenUV.add(n.xz.mul(0.01)).toVar(),
      viewZ = positionView.z.negate();
    const eye = (d) =>
      near
        .mul(far)
        .mul(2)
        .div(far.add(near).sub(d.mul(2).sub(1).mul(far.sub(near))));
    const opaqueZ = eye(depthNode.sample(uv).r).toVar();
    If(opaqueZ.lessThan(viewZ.add(0.005)), () => {
      uv.assign(screenUV);
      opaqueZ.assign(eye(depthNode.sample(uv).r));
    });
    const distance = cameraPosition.sub(p).length(),
      path = opaqueZ.sub(viewZ).mul(distance).div(viewZ.max(0.05)).clamp(0, 150),
      trans = vec3(0.29, 0.071, 0.036).mul(path).negate().exp();
    const irradiance = U.uDiffuse.add(U.uDirect.mul(U.uSun.y.max(0))).div(Math.PI),
      scatter = irradiance.mul(vec3(0.009, 0.037, 0.042));
    const body = opaqueNode
      .sample(uv)
      .rgb.mul(trans)
      .add(scatter.mul(vec3(1).sub(trans)));
    const reflected = system.clear(v.negate().reflect(n)),
      clip = mirrorMatrix.mul(vec4(p, 1)),
      muv = clip.xy.div(clip.w).add(n.xz.mul(0.006));
    const margin = N.min(N.min(muv.x, muv.y), N.min(float(1).sub(muv.x), float(1).sub(muv.y))),
      blend = margin.smoothstep(0.002, 0.02).mul(float(1).sub(distance.smoothstep(550, 1800)));
    const ref = N.mix(
      reflected,
      mirrorNode.sample(vec2(muv.x, float(1).sub(muv.y)).clamp(0.001, 0.999)).bias(0.65).rgb,
      blend,
    );
    let colour = N.mix(body, ref, F);
    const h = v.add(U.uSun).normalize(),
      nl = n.dot(U.uSun).max(0),
      nh = n.dot(h).max(0),
      vh = v.dot(h).max(0),
      rough = U.uWind.mul(0.006).add(0.042),
      alpha2 = rough.pow(4).max(0.00465421 ** 2),
      den = nh.mul(nh).mul(alpha2.sub(1)).add(1),
      D = alpha2.div(den.mul(den).mul(Math.PI)),
      k = rough.add(1).pow(2).div(8),
      G = nv.div(nv.mul(float(1).sub(k)).add(k)).mul(nl.div(nl.mul(float(1).sub(k)).add(k))),
      spec = float(1).sub(vh).pow(5).mul(0.9796).add(0.0204);
    colour = colour.add(U.uDirect.mul(D.mul(G).mul(spec).div(nv.mul(4).add(0.0001))));
    const contact = float(1).sub(path.smoothstep(0.02, 0.09)).mul(path.smoothstep(0, 0.012));
    colour = N.mix(colour, irradiance.mul(0.18), contact.mul(0.14));
    return vec4(system.display(system.fog(colour, p)), 1);
  })();
  const mesh = new T.Mesh(OM.waterGeometry(T), mat);
  mesh.name = 'tsl-refractive-water';
  mesh.frustumCulled = false;
  scene.add(mesh);
  const mc = new T.PerspectiveCamera(),
    bias = new T.Matrix4().set(0.5, 0, 0, 0.5, 0, 0.5, 0, 0.5, 0, 0, 0.5, 0.5, 0, 0, 0, 1),
    dir = new T.Vector3();
  function render() {
    const size = renderer.getDrawingBufferSize(new T.Vector2());
    opaque.setSize(size.x, size.y);
    mirror.setSize(
      Math.max(256, Math.round(size.x * 0.5)),
      Math.max(256, Math.round(size.y * 0.5)),
    );
    near.value = camera.near;
    far.value = camera.far;
    mc.copy(camera);
    mc.position.y = -camera.position.y;
    camera.getWorldDirection(dir);
    dir.y = -dir.y;
    mc.up.set(0, -1, 0);
    mc.lookAt(mc.position.clone().add(dir));
    mc.updateMatrixWorld();
    mirrorMatrix.value.copy(bias).multiply(mc.projectionMatrix).multiply(mc.matrixWorldInverse);
    const old = renderer.getRenderTarget(),
      tm = renderer.toneMapping;
    mesh.visible = false;
    U.display.value = 0;
    renderer.toneMapping = T.NoToneMapping;
    try {
      U.mirror.value = 0;
      renderer.setRenderTarget(opaque);
      renderer.render(scene, camera);
      U.mirror.value = 1;
      renderer.setRenderTarget(mirror);
      renderer.render(scene, mc);
    } finally {
      U.mirror.value = 0;
      U.display.value = 1;
      renderer.setRenderTarget(old);
      renderer.toneMapping = tm;
      mesh.visible = true;
    }
  }
  return { mesh, render, opaqueTarget: opaque, reflectionTarget: mirror };
}
