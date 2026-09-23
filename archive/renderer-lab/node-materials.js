/* TSL implementation for the renderer experiment. Geometry, field textures,
 * vegetation identities and pond hypsometry come from the production model.
 * This first path evaluates clear-sky lighting; cloud transport/rain and exact
 * leaf transmission parity remain explicit gates before promotion.
 */
import { OM } from '../../engine/om.js';
import C from '../../engine/core.js';
OM.createNodeSystem = function (T, N, atm, renderer) {
  const {
    Fn,
    float,
    int,
    vec2,
    vec3,
    vec4,
    uniform,
    texture,
    positionWorld,
    positionWorldDirection,
    positionLocal,
    positionGeometry,
    normalViewGeometry,
    normalWorldGeometry,
    cameraPosition,
    cameraViewMatrix,
    modelWorldMatrix,
    uv,
    attribute,
    output,
    If,
  } = N;
  const U = {};
  for (const [key, u] of Object.entries(atm.uniforms)) {
    if (key === 'uSkyA' || key === 'uSkyB') U[key] = texture(u.value);
    else if (u.value && typeof u.value === 'object' && 'x' in u.value) U[key] = uniform(u.value);
    else if (typeof u.value === 'number' || typeof u.value === 'boolean') U[key] = uniform(u.value);
  }
  U.display = uniform(1);
  U.mirror = uniform(0);
  function sync() {
    for (const [key, u] of Object.entries(atm.uniforms))
      if (U[key] && typeof u.value !== 'object') U[key].value = u.value;
  }
  const clear = Fn(([dir]) => {
    const d = dir.normalize(),
      e = d.y.clamp(-1, 1).asin();
    const v = float(0.5).add(
      e
        .sign()
        .mul(
          e
            .abs()
            .div(Math.PI * 0.5)
            .sqrt(),
        )
        .mul(0.5),
    );
    const az = d.x.mul(U.uSide).div(d.xz.length().max(0.000001)).clamp(-1, 1).acos().div(Math.PI);
    const st = vec2(az.mul(32).add(0.5).div(33), v.mul(40).add(0.5).div(41));
    return N.mix(U.uSkyA.sample(st).rgb, U.uSkyB.sample(st).rgb, U.uBlend);
  });
  const fog = Fn(([colour, p]) => {
    const ray = p.sub(cameraPosition),
      tr = U.uLocalExtinction.mul(ray.length()).negate().exp();
    return N.mix(clear(ray.normalize()), colour, tr);
  });
  const display = (colour) => colour.mul(N.mix(vec3(1), U.uWhiteBalance, U.display));
  const sky = Fn(() => {
    const d = positionWorldDirection.normalize(),
      angle = d.dot(U.uSun).clamp(-1, 1).acos(),
      pixel = angle.fwidth().max(0.00001);
    const disk = float(1).sub(
      N.smoothstep(float(0.00465421).sub(pixel), float(0.00465421).add(pixel), angle),
    );
    const visibility = U.uSun.y.smoothstep(-0.00465421, 0.00465421);
    const solar = U.uDirect.div(Math.PI * 0.00465421 ** 2);
    const glare = angle
      .mul(angle)
      .div(0.00013)
      .negate()
      .exp()
      .mul(0.000055)
      .add(angle.mul(angle).div(0.0018).negate().exp().mul(0.000006));
    return display(clear(d).add(solar.mul(disk.add(glare)).mul(visibility)));
  });
  const noise = Fn(([p]) => {
    const i = p.floor(),
      f = p.fract().toVar();
    f.assign(f.mul(f).mul(float(3).sub(f.mul(2))));
    const hash = (q) => q.dot(vec2(127.1, 311.7)).sin().mul(43758.5453).fract();
    return N.mix(
      N.mix(hash(i), hash(i.add(vec2(1, 0))), f.x),
      N.mix(hash(i.add(vec2(0, 1))), hash(i.add(1)), f.x),
      f.y,
    );
  });
  const gridUV = (p, b) => p.sub(b.xy).div(b.z).add(0.5).div(b.w);
  const inside = (p, b) => {
    const q = p.sub(b.xy).div(b.z),
      a = q.smoothstep(0, 6),
      c = vec2(1).sub(q.smoothstep(b.w.sub(7), b.w.sub(1)));
    return a.x.mul(a.y).mul(c.x).mul(c.y);
  };
  let fieldCache = null;
  function fields(state) {
    if (fieldCache) return fieldCache;
    const f = state.fieldUniforms;
    fieldCache = {
      bathy: texture(f.uBathymetry.value),
      mat: texture(f.uMaterialFields.value),
      env: texture(f.uEnvironmentFields.value),
      water: texture(f.uSurfaceState.value),
      pond: texture(f.uPondLevels.value),
      bounds: uniform(f.uFieldBounds.value),
      waterBounds: uniform(f.uStateBounds.value),
      pondBounds: uniform(f.uPondBounds.value),
      debug: uniform(0),
      albedo: texture(OM.surfaceAssets.albedo),
      packed: texture(OM.surfaceAssets.packed),
    };
    return fieldCache;
  }
  const layer = (F, p, id) => {
    const a = F.albedo.sample(p).depth(int(id)),
      m = F.packed.sample(p).depth(int(id)),
      xy = m.rg.mul(2).sub(1);
    return {
      a: a.rgb,
      ao: a.a,
      r: m.b,
      s: xy.div(float(1).sub(xy.dot(xy)).max(0.02).sqrt().max(0.2)),
    };
  };
  function factory(T, unused, state, kind, options = {}) {
    const F = fields(state),
      mat = new T.MeshStandardNodeMaterial({ color: 0xffffff, roughness: 0.9, ...options }),
      p = positionWorld,
      coord = p.xz;
    const env = F.env.sample(gridUV(coord, F.bounds)),
      water = F.water.sample(gridUV(coord, F.waterBounds)).mul(inside(coord, F.waterBounds));
    const elevation = p.y.add(coord.dot(coord).div(U.uR.mul(2)));
    const roof = p.x
      .sub(20)
      .abs()
      .lessThan(5.7)
      .and(p.z.sub(33).abs().lessThan(4.2))
      .and(p.y.lessThan(8));
    let colour = uniform(mat.color),
      rough = float(mat.roughness),
      ao = float(1),
      gradient = vec3(0),
      wet = water.r;
    if (kind === 'terrain') {
      const w = N.mix(
        attribute('omWeights', 'vec4'),
        F.mat.sample(gridUV(coord, F.bounds)),
        inside(coord, F.bounds),
      )
        .max(0)
        .toVar();
      const weights = w.div(w.dot(vec4(1)).max(0.001));
      const layers = [
          layer(F, coord.div(0.75), 0),
          layer(F, coord.div(1.5), 1),
          layer(F, coord.div(2.8), 2),
          layer(F, coord.div(1.1), 3),
        ],
        litter = layer(F, coord.div(1.1), 4);
      const cover = env.g.mul(inside(coord, F.bounds)).mul(weights.a.smoothstep(0.08, 0.42));
      layers[3] = {
        a: N.mix(layers[3].a, litter.a, cover),
        ao: N.mix(layers[3].ao, litter.ao, cover),
        r: N.mix(layers[3].r, litter.r, cover),
        s: N.mix(layers[3].s, litter.s, cover),
      };
      const herbs = coord
        .sub(cameraPosition.xz)
        .length()
        .smoothstep(28, 85)
        .mul(float(1).sub(cover))
        .mul(elevation.smoothstep(1.3, 4))
        .mul(noise(coord.mul(0.05)).mul(0.38).add(0.45));
      layers[3].a = N.mix(
        layers[3].a,
        vec3(0.102, 0.153, 0.035).mul(noise(coord.mul(0.13)).mul(0.38).add(0.78)),
        herbs,
      );
      colour = vec3(0);
      rough = float(0);
      ao = float(0);
      let slope = vec2(0);
      for (let i = 0; i < 4; i++) {
        const weight = weights.element(i);
        colour = colour.add(layers[i].a.mul(weight));
        rough = rough.add(layers[i].r.mul(weight));
        ao = ao.add(layers[i].ao.mul(weight));
        slope = slope.add(layers[i].s.mul(weight));
      }
      colour = colour.mul(noise(coord.mul(0.035)).mul(0.2).add(0.9));
      const damp = float(1).sub(elevation.smoothstep(-0.03, 0.24)),
        head = F.pond.sample(gridUV(coord, F.pondBounds));
      const basinWet = head.a
        .mul(inside(coord, F.pondBounds))
        .mul(float(1).sub(elevation.sub(head.r).smoothstep(-0.001, 0.025)));
      wet = water.r.max(damp).max(basinWet);
      colour = colour.mul(
        N.mix(float(1), float(0.66), water.g.mul(weights.a.add(weights.r.mul(0.65)))),
      );
      const footprint = N.max(N.dFdx(coord).length(), N.dFdy(coord).length());
      gradient = vec3(slope.x, 0, slope.y)
        .mul(0.72)
        .mul(float(1).sub(footprint.smoothstep(0.035, 0.25)));
      ao = ao.mul(N.mix(float(1), float(0.68), env.g)); // Canopy-light proxy, applied to indirect lighting.
    } else if (kind === 'rock') {
      const wn = normalWorldGeometry.abs().pow(4).toVar(),
        weight = wn.div(wn.x.add(wn.y).add(wn.z).max(0.001));
      const a = layer(F, p.xz.div(2.8), 2),
        b = layer(F, p.zy.div(2.8), 2),
        c = layer(F, p.xy.div(2.8), 2);
      colour = a.a.mul(weight.y).add(b.a.mul(weight.x)).add(c.a.mul(weight.z));
      const lichen = noise(coord.mul(0.6))
        .smoothstep(0.64, 0.83)
        .mul(weight.y.smoothstep(0.35, 0.9))
        .mul(0.16);
      colour = N.mix(colour, vec3(0.19, 0.21, 0.09), lichen);
      rough = a.r.mul(weight.y).add(b.r.mul(weight.x)).add(c.r.mul(weight.z));
      ao = a.ao.mul(weight.y).add(b.ao.mul(weight.x)).add(c.ao.mul(weight.z));
      gradient = vec3(a.s.x, 0, a.s.y)
        .mul(weight.y)
        .add(vec3(0, b.s.y, b.s.x).mul(weight.x))
        .add(vec3(c.s.x, c.s.y, 0).mul(weight.z))
        .mul(0.7);
      wet = water.r.max(float(1).sub(elevation.smoothstep(-0.03, 0.16)));
    } else if (kind === 'wood') {
      // Bark metric UVs are prepared from the same branch radii and lengths below.
      const a = layer(F, uv().mul(attribute('omBarkScale', 'vec2')).div(1.1), 5);
      colour = a.a;
      rough = a.r;
      ao = a.ao;
      gradient = vec3(a.s.x.mul(0.25), 0, a.s.y.mul(0.1));
      wet = water.r.mul(0.7);
    } else if (kind === 'timber') {
      colour = N.mix(
        vec3(0.085, 0.055, 0.029),
        vec3(0.2, 0.135, 0.071),
        noise(vec2(p.x.mul(12).add(p.z.mul(7)), p.y.mul(0.65))),
      );
      rough = float(0.85);
    } else {
      const tex = uv(),
        mid = tex.x.sub(0.5).abs().mul(-85).exp(),
        vein = tex.y
          .mul(13)
          .add(tex.x.sub(0.5).abs().mul(8))
          .mul(2 * Math.PI)
          .cos()
          .mul(0.5)
          .add(0.5)
          .pow(18)
          .mul(tex.x.sub(0.5).abs().smoothstep(0.015, 0.18));
      colour = colour
        .mul(N.mix(float(0.83), float(1.05), tex.y))
        .mul(float(1).sub(mid.mul(0.11)).sub(vein.mul(0.08)));
      rough = float(0.76);
      wet = water.a;
    }
    const wetness = N.select(roof, float(0), wet);
    mat.colorNode = colour.mul(N.mix(float(1), float(0.72), wetness));
    mat.roughnessNode = N.mix(rough, rough.mul(0.2).max(0.055), wetness);
    mat.aoNode = ao;
    const viewGradient = cameraViewMatrix.mul(vec4(gradient, 0)).xyz,
      baseNormal = N.normalView;
    mat.normalNode = baseNormal
      .add(viewGradient)
      .sub(baseNormal.mul(baseNormal.dot(viewGradient)))
      .normalize();
    mat.outputNode = vec4(display(fog(output.rgb, p)), output.a);
    mat.maskNode = U.mirror.lessThan(0.5).or(positionWorld.y.greaterThan(0.015));
    state.materials.push(mat);
    mat.userData.nodeKind = kind;
    return mat;
  }
  function refresh(state) {
    const F = fields(state);
    F.pond.value = state.fieldUniforms.uPondLevels.value;
    F.debug.value = state.fieldUniforms.uLandscapeDebug.value;
  }
  function pondFactory() {
    const mat = new T.NodeMaterial({ transparent: true, depthWrite: false, side: T.DoubleSide }),
      p = positionWorld,
      d = attribute('pondDepth', 'float');
    mat.fragmentNode = Fn(() => {
      const n = vec3(p.x.div(U.uR), 1, p.z.div(U.uR)).normalize(),
        v = cameraPosition.sub(p).normalize(),
        nv = n.dot(v).max(0.03),
        F = float(1).sub(nv).pow(5).mul(0.97963).add(0.02037),
        absorb = d.max(0).mul(-0.28).div(nv).exp().oneMinus(),
        alpha = F.add(F.oneMinus().mul(absorb));
      const refl = clear(v.negate().reflect(n)),
        scatter = U.uDiffuse.div(Math.PI).mul(vec3(0.011, 0.039, 0.042));
      let colour = refl.mul(F).add(scatter.mul(F.oneMinus()).mul(absorb)).div(alpha.max(0.0001));
      const h = v.add(U.uSun).normalize();
      colour = colour.add(U.uDirect.mul(n.dot(h).max(0).pow(4000)).mul(0.12));
      return vec4(display(fog(colour, p)), alpha.mul(d.smoothstep(0, 0.0015)));
    })();
    return mat;
  }
  function envelope(p, wind, state) {
    const F = fields(state),
      depth = N.mix(
        float(12),
        F.bathy.sample(gridUV(p, F.bounds)).r.negate().max(0),
        inside(p, F.bounds),
      ),
      amplitude = wind
        .mul(0.18)
        .add(0.4)
        .mul(C.WATER_BANDS.reduce((sum, k) => sum + 0.095 * Math.pow(0.095 / k, 1.12), 0));
    return depth.div(depth.add(amplitude.div(0.45)));
  }
  return {
    U,
    sync,
    clear,
    fog,
    display,
    sky: sky(),
    fields,
    factory,
    refresh,
    pondFactory,
    envelope,
  };
};
