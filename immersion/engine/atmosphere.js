import * as THREE from 'three';
import { OM } from './om.js';
import { MOON_RADIUS } from '../../shared/constants.js';
import C from './core.js';
import Column from './cloud-renderer.js';
const GLSL_NOISE = `
float omHash(vec3 p){return fract(sin(dot(p,vec3(127.1,311.7,74.7)))*43758.5453);}
float omNoise(vec3 p){vec3 i=floor(p),f=fract(p);f=f*f*(3.-2.*f);return mix(mix(mix(omHash(i),omHash(i+vec3(1,0,0)),f.x),mix(omHash(i+vec3(0,1,0)),omHash(i+vec3(1,1,0)),f.x),f.y),mix(mix(omHash(i+vec3(0,0,1)),omHash(i+vec3(1,0,1)),f.x),mix(omHash(i+vec3(0,1,1)),omHash(i+vec3(1,1,1)),f.x),f.y),f.z);}
float omFBM(vec3 p){return .58*omNoise(p)+.28*omNoise(p*2.03+7.3)+.14*omNoise(p*4.07+11.2);}
`;
const SKY_UNIFORMS =
  `
uniform sampler2D uSkyA,uSkyB;uniform float uBlend,uTime,uCover,uTau,uCloudBase,uCloudDepth,uWind,uVisibility,uR,uSide,uSteps,uEyeHeight,uDrift;
uniform vec3 uSun,uDirect,uDiffuse,uCloudDirect,uCloudDiffuse,uLocalExtinction,uWhiteBalance;
uniform vec2 uObserver;uniform bool uShowDisk;
uniform sampler2D uEarthSkyA,uEarthSkyB,uEarthDay,uEarthClouds;uniform float uEarthBlend,uEarthScale,uEarthShow,uEarthCos,uEarthSin,uEarthGain,uEarthCloudAlbedo,uEarthLod,uEarthDisplay;
uniform vec2 uEarthSide;uniform vec3 uEarth,uEarthSunlit,uEarthHaze;uniform mat3 uEarthBasis;
const float OM_PI=3.141592653589793;
` + Column.UNIFORMS;
const CLOUDS = `
float omShell(vec3 ro,vec3 rd,float r){float b=dot(ro,rd),c=dot(ro,ro)-r*r;return -b+sqrt(max(0.,b*b-c));}
${Column.SHADOW_GLSL}
float omDensity(vec3 p){float h=(p.y-uCloudBase)/uCloudDepth;float shape=smoothstep(0.,.16,h)*(1.-smoothstep(.63,1.,h));vec3 q=vec3(p.x+uDrift,p.y*.68,p.z+uDrift*.23)*.0011;float n=omFBM(q);float threshold=mix(.81,.12,uCover);float mass=smoothstep(threshold-.045,threshold+.20,n);float detail=.72+.28*omNoise(q*5.);return max(0.,mass*shape*detail);}
float omCloudShadow(vec3 w){if(uColumnMode>.5)return omColumnShadow(w);if(uCover<.005||uTau<.001||uSun.y<-.01)return 1.;vec3 ro=vec3(w.x,uR+max(w.y,0.),w.z);float l0=omShell(ro,uSun,uR+uCloudBase),l1=omShell(ro,uSun,uR+uCloudBase+uCloudDepth);float depth=0.;for(int i=0;i<4;i++){float s=mix(l0,l1,(float(i)+.5)/4.);vec3 p=ro+uSun*s;p.y=length(p)-uR;depth+=omDensity(p)*(l1-l0)*.25/uCloudDepth;}return exp(-min(30.,uTau*depth));}
`;
const CLEAR_LOOKUP = `vec3 omAtlas(vec3 d,sampler2D a,sampler2D b,float blend,vec2 side){float e=asin(clamp(d.y,-1.,1.));float v=.5+.5*sign(e)*sqrt(abs(e)/(OM_PI*.5));float az=acos(clamp(dot(d.xz,side)/max(length(d.xz),.000001),-1.,1.))/OM_PI;vec2 uv=vec2((az*32.+.5)/33.,(v*40.+.5)/41.);return mix(texture2D(a,uv).rgb,texture2D(b,uv).rgb,blend);}
vec3 omEarthDisk(vec3 d){float c=dot(d,uEarth);if(uEarthShow<.5||c<uEarthCos)return vec3(0.);vec3 perp=d-uEarth*c;float r=length(perp)/uEarthSin;float aa=max(fwidth(r),.0005);float cover=1.-smoothstep(1.-aa,1.+aa,r);float s=min(r,1.);vec3 p=r>1e-6?normalize(perp):vec3(0.);vec3 n=-uEarth*sqrt(max(0.,1.-s*s))+p*s;vec3 q=uEarthBasis*n;vec2 uv=vec2(atan(q.y,q.x)/(2.*OM_PI)+.5,asin(clamp(q.z,-1.,1.))/OM_PI+.5);vec3 ground=textureLod(uEarthDay,uv,uEarthLod).rgb;float cloud=textureLod(uEarthClouds,uv,uEarthLod).r;vec3 albedo=min(vec3(1.),uEarthGain*mix(ground,vec3(uEarthCloudAlbedo),cloud)+uEarthHaze);return albedo/OM_PI*uEarthSunlit*max(0.,dot(n,uSun))*cover*uEarthDisplay;}
vec3 omEarthLight(vec3 d){return uEarthScale>0.?uEarthScale*omAtlas(d,uEarthSkyA,uEarthSkyB,uEarthBlend,uEarthSide)+omEarthDisk(d):vec3(0.);}
vec3 omSunSky(vec3 d){float e=asin(clamp(d.y,-1.,1.));float v=.5+.5*sign(e)*sqrt(abs(e)/(OM_PI*.5));float az=acos(clamp(d.x*uSide/max(length(d.xz),.000001),-1.,1.))/OM_PI;vec2 uv=vec2((az*32.+.5)/33.,(v*40.+.5)/41.);return mix(texture2D(uSkyA,uv).rgb,texture2D(uSkyB,uv).rgb,uBlend);}
vec3 omClear(vec3 d){return omSunSky(d)+omEarthLight(d);}
vec3 omAirColour(vec3 d){return uColumnMode>.5&&uColumnMode<1.5?uColumnFogColour:omClear(d);}
`;
const SKY_FRAGMENT =
  SKY_UNIFORMS +
  GLSL_NOISE +
  CLOUDS +
  CLEAR_LOOKUP +
  Column.CACHED_SKY_GLSL +
  `
varying vec3 vDir;
vec3 omSky(vec3 d){if(uColumnMode>.5)return omColumnSky(d);vec3 clear=omClear(d);if(d.y<-.002)return clear;
 vec3 ro=vec3(uObserver.x,uR+uEyeHeight,uObserver.y);float l0=omShell(ro,d,uR+uCloudBase),l1=omShell(ro,d,uR+uCloudBase+uCloudDepth);float ds=(l1-l0)/uSteps;
 vec3 cloud=vec3(0.);float tr=1.;float jitter=.5;
 for(int i=0;i<32;i++){if(uCover<.001||float(i)>=uSteps||tr<.004)break;float s=l0+(float(i)+jitter)*ds;vec3 p=ro+d*s;p.y=length(p)-uR;float density=omDensity(p),dt=min(20.,density*uTau*ds/uCloudDepth);float a=1.-exp(-dt);float h=clamp((p.y-uCloudBase)/uCloudDepth,0.,1.);float topPath=(1.-h)*uTau*density/max(.10,uSun.y);
 float nu=dot(d,uSun),g=.65;float hg=(1.-g*g)/pow(max(.05,1.+g*g-2.*g*nu),1.5);
 vec3 source=uCloudDiffuse/OM_PI*(.30+.35*h)+uCloudDirect*(.025+.085*exp(-topPath*.26))*(.75+.18*hg)*smoothstep(-.15,.03,uSun.y);
 cloud+=tr*a*source;tr*=1.-a;}
 vec3 light=clear*tr+cloud;
 float k=3.912/max(80.,uVisibility),fogPath=(uVisibility<500.?55.:350.)/max(.03,d.y);float fogTr=exp(-k*fogPath);vec3 fogColour=(uDiffuse+uDirect*max(0.,uSun.y)*.14)/OM_PI;
 if(uVisibility<35000.)light=mix(fogColour,light,fogTr);else fogTr=1.;
 if(uShowDisk){float angular=acos(clamp(dot(d,uSun),-1.,1.));float pixel=max(fwidth(angular),.00001);float disk=1.-smoothstep(.00465421-pixel,.00465421+pixel,angular);float visibility=smoothstep(-.00465421,.00465421,uSun.y);
 vec3 solar=uDirect/(OM_PI*.00465421*.00465421);light+=solar*disk*visibility*tr*fogTr;
 // Display-only glare. The physical disk stays at 0.5333 degrees.
 light+=solar*(.000055*exp(-angular*angular/.00013)+.000006*exp(-angular*angular/.0018))*visibility*tr*fogTr;}
 return max(light,vec3(0.));}
void main(){vec3 d=normalize(vDir);gl_FragColor=vec4(omSky(d),1.);\n#ifdef TONE_MAPPING\n gl_FragColor.rgb*=uWhiteBalance;\n#endif\n#include <tonemapping_fragment>\n#include <colorspace_fragment>
}
`;
// Put preprocessor directives at the start of a line.
const SKY_FS = SKY_FRAGMENT.replace('1.);#include', '1.);\n#include');
// A light meter: the sky's irradiance on level ground, integrated over the upper
// hemisphere with the same sky function (clouds, fog and earth-glow included).
const METER_FS =
  SKY_FS.split('void main(){')[0].replace('varying vec3 vDir;', '') +
  `void main(){vec3 e=vec3(0.);for(int i=0;i<12;i++){float mu=(float(i)+.5)/12.,s=sqrt(1.-mu*mu);for(int j=0;j<24;j++){float phi=(float(j)+.5)/24.*6.283185307;e+=omSky(vec3(s*cos(phi),mu,s*sin(phi)))*mu;}}gl_FragColor=vec4(e*6.283185307/288.,1.);}`;
function half(v) {
  const s = v & 0x8000 ? -1 : 1,
    e = (v >> 10) & 31,
    f = v & 1023;
  return s * (e ? Math.pow(2, e - 15) * (1 + f / 1024) : (Math.pow(2, -14) * f) / 1024);
}
async function unpack(w) {
  const bytes = Uint8Array.from(atob(w.data), (c) => c.charCodeAt(0));
  if (!globalThis.DecompressionStream)
    throw Error('This browser needs DecompressionStream support to open the embedded sky data.');
  const raw = new Uint8Array(
    await new Response(
      new Blob([bytes]).stream().pipeThrough(new DecompressionStream('gzip')),
    ).arrayBuffer(),
  );
  const u = new DataView(raw.buffer),
    a = new Float32Array(raw.length / 2),
    frame = w.width * w.height * 3;
  for (let i = 0; i < a.length; i++)
    a[i] = (half(u.getUint16(i * 2, true)) * w.scale[Math.floor(i / frame)]) / 8500;
  return a;
}
class Atmosphere {
  constructor(THREE, data) {
    this.T = THREE;
    this.data = data;
    this.world = 'moon';
    this.whiteBalanceStrength = 1;
    this.arrays = {};
    this.index = -99;
    this.envTime = -1e9;
    this.envRT = null;
    this.clearLux = 0;
    this.cacheKey = '';
    this.uniforms = {
      uSkyA: { value: null },
      uSkyB: { value: null },
      uBlend: { value: 0 },
      uTime: { value: 0 },
      uCover: { value: 0 },
      uTau: { value: 0 },
      uCloudBase: { value: 1000 },
      uCloudDepth: { value: 1000 },
      uWind: { value: 2 },
      uVisibility: { value: 35000 },
      uR: { value: MOON_RADIUS },
      uSide: { value: 1 },
      uSteps: { value: 14 },
      uSun: { value: new THREE.Vector3(0, 1, 0) },
      uDirect: { value: new THREE.Vector3() },
      uDiffuse: { value: new THREE.Vector3() },
      uWhiteBalance: { value: new THREE.Vector3(1, 1, 1) },
      uLocalExtinction: { value: new THREE.Vector3(0.000007, 0.000016, 0.000037) },
      uCloudDirect: { value: new THREE.Vector3() },
      uCloudDiffuse: { value: new THREE.Vector3() },
      uObserver: { value: new THREE.Vector2() },
      uShowDisk: { value: true },
      uEyeHeight: { value: 1.7 },
      uDrift: { value: 0 },
      uEarthSkyA: { value: null },
      uEarthSkyB: { value: null },
      uEarthDay: { value: null },
      uEarthClouds: { value: null },
      uEarthBlend: { value: 0 },
      uEarthScale: { value: 0 },
      uEarthShow: { value: 0 },
      uEarthCos: { value: 1 },
      uEarthSin: { value: 0.0166 },
      uEarthGain: { value: 1 },
      uEarthCloudAlbedo: { value: 0.75 },
      uEarthLod: { value: 0 },
      uEarthDisplay: { value: 1 },
      uEarthSide: { value: new THREE.Vector2(-1, 0) },
      uEarth: { value: new THREE.Vector3(-1, 0, 0) },
      uEarthSunlit: { value: new THREE.Vector3() },
      uEarthHaze: { value: new THREE.Vector3() },
      uEarthBasis: { value: new THREE.Matrix3() },
    };
    this.bodies = null;
    this.earth = null;
    this.earthDirectRaw = [0, 0, 0];
    this.earthLux = 0;
    this.columnClouds = new Column.Controller(THREE, this);
  }
  async init(renderer, dataOnly = false) {
    const T = this.T;
    for (const [key, w] of Object.entries(this.data.worlds)) this.arrays[key] = await unpack(w);
    for (const k of ['uSkyA', 'uSkyB', 'uEarthSkyA', 'uEarthSkyB']) {
      const tex = new T.DataTexture(
        new Float32Array(33 * 41 * 4),
        33,
        41,
        T.RGBAFormat,
        T.FloatType,
      );
      tex.minFilter = tex.magFilter = T.LinearFilter;
      tex.colorSpace = T.LinearSRGBColorSpace;
      tex.generateMipmaps = false;
      tex.needsUpdate = true;
      this.uniforms[k].value = tex;
    }
    if (dataOnly) return;
    this.material = new T.ShaderMaterial({
      uniforms: this.uniforms,
      vertexShader:
        'varying vec3 vDir;void main(){vDir=position;vec4 p=projectionMatrix*mat4(mat3(viewMatrix))*vec4(position,1.);gl_Position=p.xyww;}',
      fragmentShader: SKY_FS,
      side: T.BackSide,
      depthWrite: false,
      depthTest: true,
      toneMapped: true,
    });
    this.mesh = new T.Mesh(new T.SphereGeometry(1, 24, 16), this.material);
    this.mesh.frustumCulled = false;
    this.mesh.renderOrder = -100;
    this.envScene = new T.Scene();
    this.envScene.add(this.mesh.clone());
    this.cubeRT = new T.WebGLCubeRenderTarget(256, {
      type: T.HalfFloatType,
      generateMipmaps: true,
      minFilter: T.LinearMipmapLinearFilter,
    });
    this.cubeCamera = new T.CubeCamera(0.1, 10, this.cubeRT);
    this.pmrem = new T.PMREMGenerator(renderer);
    // Without float render targets the meter is off and exposure uses the clear-sky light.
    this.meterTarget = renderer.extensions.has('EXT_color_buffer_float')
      ? new T.WebGLRenderTarget(1, 1, { type: T.FloatType, depthBuffer: false })
      : null;
    this.meterMaterial = new T.ShaderMaterial({
      uniforms: this.uniforms,
      vertexShader: 'void main(){gl_Position=vec4(position.xy,0.,1.);}',
      fragmentShader: METER_FS,
      depthTest: false,
      depthWrite: false,
      toneMapped: false,
    });
    const meterQuad = new T.Mesh(new T.PlaneGeometry(2, 2), this.meterMaterial);
    meterQuad.frustumCulled = false;
    this.meterScene = new T.Scene();
    this.meterScene.add(meterQuad);
    this.meterCamera = new T.OrthographicCamera(-1, 1, 1, -1, 0, 1);
    this.meterPixel = new Float32Array(4);
    this.meterTime = -1e9;
    this.columnClouds.initRenderer(renderer, SKY_UNIFORMS + GLSL_NOISE + CLOUDS + CLEAR_LOOKUP);
  }
  update(world, phase, w, time, observer, drift = 0) {
    const T = this.T;
    this.world = world;
    const d = this.data.worlds[world],
      s = C.sunAt(phase),
      arr = this.arrays[world];
    let lo = 0,
      hi = d.suns.length - 1;
    while (hi - lo > 1) {
      const m = (hi + lo) >> 1;
      if (d.suns[m] <= s.elevation) lo = m;
      else hi = m;
    }
    const blend = C.clamp((s.elevation - d.suns[lo]) / (d.suns[hi] - d.suns[lo]));
    const key = world + ':' + lo;
    if (key !== this.cacheKey) {
      this.uploadFrames(arr, lo, hi, 'uSkyA', 'uSkyB');
      this.cacheKey = key;
    }
    const interp = (k) => d[k][lo].map((v, i) => C.mix(v, d[k][hi][i], blend));
    this.directRaw = interp('direct');
    this.diffuseRaw = interp('diffuse');
    const horiz = interp('direct_horizontal');
    this.clearLux = horiz
      .map((v, i) => v + this.diffuseRaw[i])
      .reduce((s, v, i) => s + v * [0.2126, 0.7152, 0.0722][i], 0);
    for (const [key, field] of [
      ['uDirect', 'direct'],
      ['uDiffuse', 'diffuse'],
      ['uCloudDirect', 'cloud_direct'],
      ['uCloudDiffuse', 'cloud_diffuse'],
    ]) {
      const v = interp(field);
      this.uniforms[key].value.set(
        Math.max(0, v[0]) / 8500,
        Math.max(0, v[1]) / 8500,
        Math.max(0, v[2]) / 8500,
      );
    }
    const u = this.uniforms;
    const E = d.direct.at(-1).map((v, i) => v + d.diffuse.at(-1)[i]),
      Y = E.reduce((sum, v, i) => sum + v * [0.2126, 0.7152, 0.0722][i], 0);
    u.uWhiteBalance.value.set(
      ...E.map((v) => C.mix(1, Y / Math.max(1e-8, v), this.whiteBalanceStrength)),
    );
    const density = world === 'earth' ? 1 : 1.2;
    // RGB-band local molecular extinction plus separately prescribed haze.
    // Selected proxy, not a new spectral transport solve.
    const haze = Math.max(0, 3.912 / w.visibility - 3.912 / 180000);
    u.uLocalExtinction.value.set(
      density * 0.0000058 + haze,
      density * 0.0000135 + haze,
      density * 0.0000331 + haze,
    );
    u.uBlend.value = blend;
    u.uSun.value.set(s.x, s.y, s.z);
    u.uSide.value = s.side;
    u.uTime.value = time;
    u.uDrift.value = drift;
    u.uCover.value = w.coverage;
    u.uTau.value = w.tau;
    u.uCloudBase.value = w.cloudBase;
    u.uCloudDepth.value = w.thickness;
    u.uWind.value = w.wind;
    u.uVisibility.value = w.visibility;
    u.uR.value = d.radius;
    u.uObserver.value.set(observer.x, observer.z);
    u.uEyeHeight.value = Math.max(0.1, observer.y);
    this.updateEarth(world, phase);
    this.columnClouds.update(world);
    return s;
  }
  uploadFrames(arr, lo, hi, a, b) {
    for (const [name, idx] of [
      [a, lo],
      [b, hi],
    ]) {
      const dest = this.uniforms[name].value.image.data,
        offset = idx * 33 * 41 * 3;
      for (let i = 0; i < 33 * 41; i++) {
        dest[i * 4] = arr[offset + i * 3];
        dest[i * 4 + 1] = arr[offset + i * 3 + 1];
        dest[i * 4 + 2] = arr[offset + i * 3 + 2];
        dest[i * 4 + 3] = 1;
      }
      this.uniforms[name].value.needsUpdate = true;
    }
  }
  /** Bracketing atlas frames and blend for a source at the given elevation (degrees). */
  bracket(d, elevation) {
    let lo = 0,
      hi = d.suns.length - 1;
    while (hi - lo > 1) {
      const m = (hi + lo) >> 1;
      if (d.suns[m] <= elevation) lo = m;
      else hi = m;
    }
    return { lo, hi, blend: C.clamp((elevation - d.suns[lo]) / (d.suns[hi] - d.suns[lo])) };
  }
  /** Sky bodies (engine/sky-bodies.js) and the Earth imagery with its calibration. */
  async setSkyBodies(bodies, earth) {
    const T = this.T,
      u = this.uniforms;
    this.bodies = bodies;
    if (!earth) return;
    const loader = new T.TextureLoader();
    const [day, clouds] = await Promise.all([
      loader.loadAsync(earth.day),
      loader.loadAsync(earth.clouds),
    ]);
    day.colorSpace = T.SRGBColorSpace;
    clouds.colorSpace = T.NoColorSpace;
    for (const t of [day, clouds]) {
      t.wrapS = T.RepeatWrapping;
      t.minFilter = T.LinearMipmapLinearFilter;
      t.generateMipmaps = true;
      t.needsUpdate = true;
    }
    u.uEarthDay.value = day;
    u.uEarthClouds.value = clouds;
    u.uEarthGain.value = earth.calibration.gain;
    u.uEarthCloudAlbedo.value = earth.calibration.cloud_albedo;
    u.uEarthHaze.value.set(...earth.calibration.haze);
    const radius = bodies.k.earth_angular_radius_rad;
    u.uEarthSin.value = Math.sin(radius);
    u.uEarthCos.value = Math.cos(radius * 1.04);
    this.earthTextureWidth = day.image.width;
  }
  /** Earth's direction, phase, glow and light for this phase of the lunar day. */
  updateEarth(world, phase) {
    const u = this.uniforms;
    if (!this.bodies || world === 'earth') {
      u.uEarthScale.value = 0;
      u.uEarthShow.value = 0;
      this.earth = null;
      this.earthDirectRaw = [0, 0, 0];
      this.earthLux = 0;
      return;
    }
    const st = this.bodies.scene(phase * this.bodies.period),
      d = this.data.worlds[world],
      [x, y, z] = st.earth,
      elevation = (Math.asin(C.clamp(y, -1, 1)) * 180) / Math.PI,
      { lo, hi, blend } = this.bracket(d, elevation);
    const key = world + ':' + lo;
    if (key !== this.earthCacheKey) {
      this.uploadFrames(this.arrays[world], lo, hi, 'uEarthSkyA', 'uEarthSkyB');
      this.earthCacheKey = key;
    }
    const interp = (k) => d[k][lo].map((v, i) => C.mix(v, d[k][hi][i], blend));
    const direct = interp('direct').map((v) => Math.max(0, v)),
      ratio = st.earthlightRatio;
    const horizontal = interp('direct_horizontal'),
      diffuse = interp('diffuse');
    this.earth = st;
    this.earthDirectRaw = direct.map((v) => v * ratio);
    this.earthLux =
      ratio *
      horizontal.reduce((sum, v, i) => sum + (v + diffuse[i]) * [0.2126, 0.7152, 0.0722][i], 0);
    this.clearLux += this.earthLux;
    u.uEarthScale.value = ratio;
    u.uEarthShow.value = u.uEarthDay.value ? 1 : 0;
    u.uEarthBlend.value = blend;
    u.uEarth.value.set(x, y, z);
    const h = Math.hypot(x, z) || 1;
    u.uEarthSide.value.set(x / h, z / h);
    u.uEarthSunlit.value.set(...direct.map((v) => v / 8500));
    // Scene frame -> Earth-fixed: undo the celestial rotation, then the Earth's spin.
    const m = st.celestialToScene,
      g = st.earthRotation,
      c = Math.cos(g),
      sn = Math.sin(g);
    const col = (j) => [m[0][j], m[1][j], m[2][j]]; // rows of the transpose
    const [r0, r1, r2] = [col(0), col(1), col(2)];
    u.uEarthBasis.value.set(
      c * r0[0] + sn * r1[0],
      c * r0[1] + sn * r1[1],
      c * r0[2] + sn * r1[2],
      -sn * r0[0] + c * r1[0],
      -sn * r0[1] + c * r1[1],
      -sn * r0[2] + c * r1[2],
      r2[0],
      r2[1],
      r2[2],
    );
  }
  /** Light for exposure (lux): metered sky irradiance on level ground, plus the direct
   * Sun or Earth beam after cloud and fog. Unlike clearLux, it sees the actual sky. */
  meter(renderer, stamp, force = false) {
    if (!this.meterTarget) return this.clearLux;
    if (!force && stamp - this.meterTime < 0.5 && Number.isFinite(this.meterLux))
      return this.meterLux;
    this.meterTime = stamp;
    this.columnClouds.prepare(renderer);
    const u = this.uniforms,
      disk = u.uShowDisk.value,
      tm = renderer.toneMapping,
      exp = renderer.toneMappingExposure,
      target = renderer.getRenderTarget();
    let sky = NaN;
    u.uShowDisk.value = false;
    renderer.toneMapping = this.T.NoToneMapping;
    renderer.toneMappingExposure = 1;
    try {
      renderer.setRenderTarget(this.meterTarget);
      renderer.render(this.meterScene, this.meterCamera);
      renderer.readRenderTargetPixels(this.meterTarget, 0, 0, 1, 1, this.meterPixel);
      const p = this.meterPixel;
      sky = (0.2126 * p[0] + 0.7152 * p[1] + 0.0722 * p[2]) * 8500;
    } finally {
      renderer.setRenderTarget(target);
      renderer.toneMapping = tm;
      renderer.toneMappingExposure = exp;
      u.uShowDisk.value = disk;
    }
    const lum = (v) =>
      0.2126 * Math.max(0, v[0]) + 0.7152 * Math.max(0, v[1]) + 0.0722 * Math.max(0, v[2]);
    const through = (mu) => {
      const m = Math.max(0.1, mu),
        v = u.uVisibility.value;
      let t = 1 - u.uCover.value * (1 - Math.exp(-u.uTau.value / m));
      if (v < 35000)
        t *= Math.exp(((-3.912 / Math.max(80, v)) * (v < 500 ? 55 : 350)) / Math.max(0.03, mu));
      return Math.max(0, t);
    };
    const sunUp = u.uSun.value.y,
      sun = lum(this.directRaw || [0, 0, 0]) * Math.max(0, sunUp) * through(sunUp);
    const earthUp = this.earth ? this.earth.earth[1] : 0,
      earth = this.earth ? lum(this.earthDirectRaw) * Math.max(0, earthUp) * through(earthUp) : 0;
    this.meterSky = sky;
    this.meterLux = Number.isFinite(sky) && sky >= 0 ? sky + sun + earth : this.clearLux;
    return this.meterLux;
  }
  /** Display-only scaling of the Earth disk so its brightest point sits near display
   * white at this exposure, as the eye's local adaptation keeps a bright disk readable.
   * Earthlight on the ground and in the sky glow is not affected. */
  setEarthDisplay(exposure) {
    const peak = (Math.max(...this.uniforms.uEarthSunlit.value.toArray()) / Math.PI) * exposure;
    this.uniforms.uEarthDisplay.value = peak > 0 ? Math.min(1, 2.5 / peak) : 1;
  }
  /** Mip level for the Earth disk given its size on screen in pixels. */
  setEarthPixels(diameterPixels) {
    const texels = (this.earthTextureWidth || 2048) / 2; // texels across the visible hemisphere
    this.uniforms.uEarthLod.value = Math.max(0, Math.log2(texels / Math.max(1, diameterPixels)));
  }
  clearEnvironmentKey() {
    const u = this.uniforms;
    if (u.uColumnMode.value > 0.5) return this.columnClouds.renderKey;
    if (u.uCover.value > 0 && u.uTau.value > 0) return null;
    return [
      this.world,
      this.cacheKey,
      u.uBlend.value,
      u.uSide.value,
      u.uVisibility.value,
      u.uR.value,
      ...u.uSun.value.toArray(),
      ...u.uDirect.value.toArray(),
      ...u.uDiffuse.value.toArray(),
      u.uEarthScale.value,
      ...u.uEarth.value.toArray(),
    ].join(':');
  }
  environment(renderer, scene, stamp, force = false) {
    this.columnClouds.prepare(renderer);
    const key = this.clearEnvironmentKey();
    if (key !== null && key === this.clearEnvKey && this.envRT) {
      scene.environment = this.envRT.texture;
      return;
    }
    if (!force && stamp - this.envTime < 2) return;
    this.envTime = stamp;
    const disk = this.uniforms.uShowDisk.value,
      tm = renderer.toneMapping,
      exp = renderer.toneMappingExposure,
      target = renderer.getRenderTarget();
    const columnSteps = this.uniforms.uColumnSteps.value;
    this.uniforms.uColumnSteps.value = Math.min(16, columnSteps);
    this.uniforms.uShowDisk.value = false;
    renderer.toneMapping = this.T.NoToneMapping;
    renderer.toneMappingExposure = 1;
    try {
      this.cubeCamera.update(renderer, this.envScene);
      const next = this.pmrem.fromCubemap(this.cubeRT.texture);
      scene.environment = next.texture;
      if (this.envRT) this.envRT.dispose();
      this.envRT = next;
      this.clearEnvKey = key;
      this.environmentUpdates = (this.environmentUpdates || 0) + 1;
    } finally {
      this.uniforms.uColumnSteps.value = columnSteps;
      renderer.setRenderTarget(target);
      renderer.toneMapping = tm;
      renderer.toneMappingExposure = exp;
      this.uniforms.uShowDisk.value = disk;
    }
  }
}
OM.Atmosphere = Atmosphere;
OM.SKY_UNIFORMS = SKY_UNIFORMS;
OM.GLSL_NOISE = GLSL_NOISE;
OM.CLOUDS = CLOUDS;
OM.CLEAR_LOOKUP = CLEAR_LOOKUP;
