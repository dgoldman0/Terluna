/* Column cloud transfer into the existing WebGL scene. Field values and lighting
 * are sampled at physical radial altitude. Subpixel morphology, multiple cloud
 * scattering and RGB aerial in-scattering are explicit rendering approximations.
 */
// The column model belongs to the atmosphere domain. It is still a browser script
// that registers itself globally; baked column products will replace this import.
import '../../atmosphere/column/weather-column.js';
import O from './cloud-optics.js';
const W = globalThis.OpenMoonWeatherColumn;
const UNIFORMS = `
uniform float uColumnMode,uColumnBase,uColumnTop,uColumnScale,uColumnCoverage,uColumnSteps,uColumnTime,uOpticalH,uOpticalDensity,uLightTop,uColumnCacheReady,uColumnShadowReady,uColumnAverageShadow;
uniform sampler2D uColumnProfile,uColumnSun,uColumnSkyMap,uColumnShadowMap;
uniform vec4 uColumnShadowBounds;
uniform highp sampler3D uColumnNoiseTexture;
uniform vec3 uColumnSolar,uOpticalBeta,uColumnFogColour;
`;
// One immutable periodic lattice replaces repeated hash/trigonometric noise in
// each view/shadow ray step. Interpolation remains a world-space density field.
function noiseLattice(n = 64) {
  if (n !== 64) throw new RangeError('The shader noise period is pinned to 64');
  const bytes = new Uint8Array(n * n * n);
  for (let z = 0; z < n; z++)
    for (let y = 0; y < n; y++)
      for (let x = 0; x < n; x++) {
        let h =
          Math.imul(x + 17, 374761393) ^
          Math.imul(y + 31, 668265263) ^
          Math.imul(z + 79, 2147483647);
        h = Math.imul(h ^ (h >>> 13), 1274126177);
        bytes[(z * n + y) * n + x] = (h ^ (h >>> 16)) >>> 24;
      }
  return bytes;
}
const DENSITY_GLSL = `
float omColumnNoise(vec3 p){vec3 a=floor(p),f=fract(p);f=f*f*(3.-2.*f);return texture(uColumnNoiseTexture,(a+f+.5)/64.).r;}
float omColumnFBM(vec3 p){return .58*omColumnNoise(p)+.28*omColumnNoise(p*2.03+7.3)+.14*omColumnNoise(p*4.07+11.2);}

// Stable radial altitude; all directions, including the reflection camera.
float omRadialHeight(vec3 ro){float q=dot(ro.xz,ro.xz)+(ro.y-uR)*(ro.y+uR);return q/(length(ro)+uR);}
vec4 omColumnAt(float h){float t=(h-uColumnBase)/max(1.,uColumnTop-uColumnBase);return texture2D(uColumnProfile,vec2((clamp(t,0.,1.)*255.+.5)/256.,.5));}
vec2 omCellHash(vec2 p){return fract(sin(vec2(dot(p,vec2(127.1,311.7)),dot(p,vec2(269.5,183.3))))*43758.5453);}
float omColumnExt(vec3 ro){
 float h=omRadialHeight(ro);if(h<uColumnBase||h>uColumnTop)return 0.;
 vec4 profile=omColumnAt(h);if(uColumnMode<1.5)return profile.r; // Fog is locally stratified.
 float y=(h-uColumnBase)/max(1.,uColumnTop-uColumnBase);
 vec2 p=ro.xz-vec2(profile.b,profile.a)*uColumnTime;
 float scale=uColumnScale;
 if(uColumnMode<2.5){
  // Continuous multiscale field gives clustered, unequal cloud bodies, with
  // no grid of identical cloud sprites. Vertical support still comes from the column.
  vec3 q=vec3(p/scale*3.8,y*1.7);
  float broad=omColumnFBM(vec3(p/scale*.57,1.7));
  float n=.70*omColumnFBM(q)+.20*omColumnNoise(q*2.91+13.)+.10*omColumnNoise(q*7.03+6.);
  float threshold=mix(.72,.36,uColumnCoverage)+.17*smoothstep(.48,1.,y)-.12*(broad-.5);
  float cloud=smoothstep(threshold-.045,threshold+.095,n);
  return profile.r*cloud;
 }
 vec2 grid=p/scale,cell=floor(grid);float mass=0.;
 float breakup=omColumnFBM(vec3(p/scale*15.,h/scale*15.));
 for(int j=-1;j<=1;j++)for(int i=-1;i<=1;i++){
  vec2 index=cell+vec2(float(i),float(j)),r=omCellHash(index);
  if(r.y>clamp(uColumnCoverage*2.1,.05,1.))continue;
  vec2 centre=index+vec2(.2)+r*.6,offset=grid-centre;
  float radius=mix(.20,.33,r.x),height=.78+.22*r.y,yy=y/height;
  // Overlapping buoyant lobes and a capped spreading head. Individual towers
  // terminate below the shared neutral-buoyancy ceiling; no periodic sine ribs.
  vec2 lean=vec2(r.x-.5,r.y-.5)*.16;
  float body=length(vec3(offset/(radius*.70),(yy-.22)/.39));
  float middle=length(vec3((offset-lean*.5)/(radius*.94),(yy-.49)/.32));
  float head=length(vec3((offset-lean)/(radius*1.08),(yy-.73)/.23));
  float anvil=length(vec3((offset-lean*1.5)/(radius*1.80),(yy-.87)/.10));
  float d=min(min(body,middle),min(head,anvil));
  float shape=1.-smoothstep(.73,1.08,d+(breakup-.5)*.72);
  mass=max(mass,shape);
 }
 float detail=mix(.73,1.,omColumnFBM(vec3(p/(scale*.19),h/(scale*.19))));
 return profile.r*mass*detail;
}
float omColumnSunDepth(vec3 ro){
 float end=max(0.,omShell(ro,uSun,uR+uColumnTop));float start=omRadialHeight(ro)<uColumnBase?max(0.,omShell(ro,uSun,uR+uColumnBase)):0.;float tau=0.;
 for(int j=0;j<8;j++){float a=pow(float(j)/8.,2.),b=pow(float(j+1)/8.,2.);tau+=omColumnExt(ro+uSun*mix(start,end,(a+b)*.5))*(end-start)*(b-a);}
 return min(40.,tau);
}

`;
const SHADOW_GLSL = `
float omColumnShadow(vec3 w){
 if(uSun.y<0.)return 1.;
 if(uColumnShadowReady>.5){
  vec2 uv=(w.xz-uColumnShadowBounds.xy)/uColumnShadowBounds.z;
  vec2 edge=smoothstep(vec2(0.),vec2(.04),uv)*(1.-smoothstep(vec2(.96),vec2(1.),uv));
  return mix(uColumnAverageShadow,texture2D(uColumnShadowMap,clamp(uv,0.,1.)).r,edge.x*edge.y);
 }
 return 1.; // Cache is prepared before scene lighting; no per-object raymarch.
}
`;
const SKY_GLSL = `
vec3 omColumnLight(vec3 ro){
 float h=max(0.,omRadialHeight(ro)),mu=dot(normalize(ro),uSun),v=.5+.5*sign(mu)*sqrt(abs(mu));
 if(mu<0.&&pow(uR+h,2.)*(1.-mu*mu)<uR*uR)return vec3(0.);
 vec2 st=vec2((v*127.+.5)/128.,(clamp(h/uLightTop,0.,1.)*63.+.5)/64.);
 return uColumnSolar*texture2D(uColumnSun,st).rgb;
}
vec3 omColumnEyeTrans(vec3 ro,vec3 dir,float s){
 float column=0.;for(int j=0;j<6;j++){float t=(float(j)+.5)/6.;column+=exp(-max(0.,omRadialHeight(ro+dir*s*t))/uOpticalH);}
 return exp(-uOpticalBeta*column*(s/6.)*uOpticalDensity);
}
vec4 omColumnTransfer(vec3 d){
 vec3 clear=omClear(d),ro=vec3(uObserver.x,uR+uEyeHeight,uObserver.y);
 float roHeight=omRadialHeight(ro),mu=dot(normalize(ro),d);
 // Stop at the physical solid horizon; below-horizon sky remains behind terrain.
 if(mu<0.&&pow(uR+roHeight,2.)*(1.-mu*mu)<uR*uR)return vec4(clear,1.);
 float entry=roHeight>=uColumnBase?0.:max(0.,omShell(ro,d,uR+uColumnBase));
 float end=max(0.,omShell(ro,d,uR+uColumnTop));
 float count=uColumnSteps,tr=1.;vec3 cloud=vec3(0.);
 for(int i=0;i<192;i++){
  if(float(i)>=count||tr<.002||end<=entry)break;
  // Quadratic sampling resolves foreground fog and cloud entry surfaces.
  float a=pow(float(i)/count,1.35),b=pow(float(i+1)/count,1.35);
  float jitter=.12+.76*omHash(vec3(gl_FragCoord.xy,float(i)+7.));
  float s=mix(entry,end,mix(a,b,jitter)),ds=(end-entry)*(b-a);vec3 pos=ro+d*s;
  float ext=omColumnExt(pos),tau=ext*ds,opacity=1.-exp(-min(30.,tau));
  if(opacity>.00002){
   float h=max(0.,omRadialHeight(pos)),sunDepth=omColumnSunDepth(pos),nu=dot(d,uSun);
   float iceFraction=omColumnAt(h).g,g=mix(.78,.82,iceFraction),phase=(1.-g*g)/(12.56637061*pow(max(.015,1.+g*g-2.*g*nu),1.5));
   float backG=-.25;phase=.85*phase+.15*(1.-backG*backG)/(12.56637061*pow(1.+backG*backG-2.*backG*nu,1.5));
   vec3 direct=omColumnLight(pos);
   // Diffuse field and softened internal illumination are selected closures.
   // Two radiative boundaries: blue downwelling sky and surface-reflected
   // upwelling light. A selected regional Lambertian albedo of 0.18 supplies
   // lower cloud faces; the old sky-only closure made those faces too blue.
   float layerY=(h-uColumnBase)/max(1.,uColumnTop-uColumnBase);
   vec3 skyDiffuse=uDiffuse/OM_PI*exp(-h/(2.*uOpticalH))*(.24+.34*exp(-sunDepth*.05));
   vec3 groundTrans=exp(-uOpticalBeta*uOpticalDensity*uOpticalH*(1.-exp(-h/uOpticalH)));
   vec3 groundDiffuse=(uDiffuse+uDirect*max(0.,uSun.y)*uColumnAverageShadow)*(.18/OM_PI)*groundTrans;
   vec3 diffuse=skyDiffuse+groundDiffuse*mix(1.,.25,smoothstep(.0,.8,layerY));
   vec3 source=diffuse+direct*(phase*exp(-sunDepth)+.09*exp(-sunDepth*.16));
   if(uColumnMode<1.5)source=uColumnFogColour;
   vec3 gasTr=omColumnEyeTrans(ro,d,s);
   cloud+=tr*opacity*(source*gasTr+clear*(1.-gasTr));tr*=1.-opacity;
  }
 }
 return vec4(max(clear*tr+cloud,vec3(0.)),tr);
}

`;
const CACHED_SKY_GLSL = `
vec3 omColumnSky(vec3 d){
 vec4 transfer;
 if(uColumnCacheReady>.5){
  float e=asin(clamp(d.y,-1.,1.)),v=.5+.5*sign(e)*sqrt(abs(e)/(OM_PI*.5));
  float az=atan(d.z,d.x)/6.283185307+.5;
  transfer=texture2D(uColumnSkyMap,vec2(az,v));
 }else transfer=vec4(omClear(d),1.);
 vec3 light=transfer.rgb;float tr=transfer.a;
 if(uShowDisk){float angular=acos(clamp(dot(d,uSun),-1.,1.)),pixel=max(fwidth(angular),.00001);
  float disk=1.-smoothstep(.00465421-pixel,.00465421+pixel,angular),visibility=smoothstep(-.00465421,.00465421,uSun.y);
  light+=uDirect/(OM_PI*.00465421*.00465421)*(disk+.000055*exp(-angular*angular/.00013)+.000006*exp(-angular*angular/.0018))*visibility*tr;
 }
 return max(light,vec3(0.));
}
`;
class Controller {
  constructor(T, atm) {
    this.T = T;
    this.atm = atm;
    this.key = 'reference';
    this.cache = new Map();
    this.current = null;
    this.world = null;
    this.renderKey = null;
    this.bakeCount = 0;
    this.lastBakeMs = 0;
    const texture = (data, w, h) => {
      const t = new T.DataTexture(data, w, h, T.RGBAFormat, T.FloatType);
      t.minFilter = t.magFilter = T.LinearFilter;
      t.generateMipmaps = false;
      t.needsUpdate = true;
      return t;
    };
    this.noiseTexture = new T.Data3DTexture(noiseLattice(), 64, 64, 64);
    this.noiseTexture.format = T.RedFormat;
    this.noiseTexture.type = T.UnsignedByteType;
    this.noiseTexture.minFilter = this.noiseTexture.magFilter = T.LinearFilter;
    this.noiseTexture.wrapS = this.noiseTexture.wrapT = this.noiseTexture.wrapR = T.RepeatWrapping;
    this.noiseTexture.unpackAlignment = 1;
    this.noiseTexture.needsUpdate = true;
    this.profileTexture = texture(new Float32Array(256 * 4), 256, 1);
    this.sunTexture = texture(new Float32Array(128 * 64 * 4), 128, 64);
    Object.assign(atm.uniforms, {
      uColumnNoiseTexture: { value: this.noiseTexture },
      uColumnCacheReady: { value: 0 },
      uColumnShadowReady: { value: 0 },
      uColumnAverageShadow: { value: 1 },
      uColumnSkyMap: { value: null },
      uColumnShadowMap: { value: null },
      uColumnShadowBounds: { value: new T.Vector4(-16384, -16384, 32768, 128) },
      uColumnMode: { value: 0 },
      uColumnBase: { value: 0 },
      uColumnTop: { value: 1 },
      uColumnScale: { value: 10000 },
      uColumnCoverage: { value: 0.4 },
      uColumnSteps: { value: 128 },
      uColumnTime: { value: 0 },
      uOpticalH: { value: 48400 },
      uOpticalDensity: { value: 1.2 },
      uOpticalBeta: { value: new T.Vector3(...O.BETA) },
      uLightTop: { value: 65000 },
      uColumnSolar: { value: new T.Vector3(1, 1, 1) },
      uColumnFogColour: { value: new T.Vector3() },
      uColumnProfile: { value: this.profileTexture },
      uColumnSun: { value: this.sunTexture },
    });
  }
  select(key) {
    if (key !== 'reference' && !W.PRESETS[key]) throw new RangeError('Unknown atmospheric study');
    this.key = key;
    this.world = null;
    this.current = null;
    this.atm.envTime = -1e9;
    this.atm.clearEnvKey = null;
    this.atm.uniforms.uColumnMode.value = 0;
    this.atm.uniforms.uColumnCacheReady.value = 0;
    this.atm.uniforms.uColumnShadowReady.value = 0;
    this.renderKey = null;
  }
  update(world) {
    const u = this.atm.uniforms;
    if (this.key === 'reference') {
      u.uColumnMode.value = 0;
      return;
    }
    if (this.world !== world || !this.current) {
      const key = world + ':' + this.key;
      let record = this.cache.get(key);
      if (!record) {
        const model = W.create(this.key, world),
          profile = W.textureData(model),
          light = O.buildTable(world, profile.top + 1000);
        record = { model, profile, light };
        this.cache.set(key, record);
      }
      this.current = record;
      this.world = world;
      this.profileTexture.image.data = record.profile.data;
      this.profileTexture.needsUpdate = true;
      this.sunTexture.image.data = record.light.data;
      this.sunTexture.needsUpdate = true;
      u.uColumnBase.value = record.profile.base;
      u.uColumnTop.value = record.profile.top;
      u.uLightTop.value = record.light.maxHeight;
      const m = record.model.inputs.morphology;
      u.uColumnScale.value = Math.max(
        300,
        (record.profile.top - record.profile.base) * m.aspect * 3.2,
      );
      u.uColumnCoverage.value = m.coverage;
      const optical = O.profile(world),
        t = O.transmission(world, 0, 1);
      u.uOpticalH.value = optical.H;
      u.uOpticalDensity.value = optical.density;
      // Anchor the RGB proxy to the inherited noon ground irradiance. Twilight
      // then follows altitude-dependent spherical ray transmittance, not a fixed height.
      u.uColumnSolar.value.set(
        ...this.atm.data.worlds[world].direct.at(-1).map((e, i) => Math.max(0, e) / 8500 / t[i]),
      );
    }
    const c = this.current.model;
    u.uColumnMode.value = c.inputs.morphology.type + 1;
    u.uColumnTime.value = u.uTime.value;
    u.uCover.value = c.inputs.morphology.coverage;
    u.uTau.value = c.summary.inCloudOpticalDepth;
    if (c.key === 'fog') {
      u.uLocalExtinction.value.addScalar(c.sample(u.uEyeHeight.value).extinction);
      u.uColumnFogColour.value
        .copy(u.uDiffuse.value)
        .multiplyScalar(0.65 / Math.PI)
        .addScaledVector(
          u.uDirect.value,
          (0.025 + 0.075 * Math.exp(-c.summary.inCloudOpticalDepth * 0.25)) / Math.PI,
        );
    }
  }
  snapshot() {
    return this.current
      ? {
          regime: this.key,
          ...this.current.model.summary,
          optics: {
            kind: 'inherited exponential molecular profile with altitude-dependent RGB cloud transport',
            ...O.profile(this.world),
            ozoneEnabled: O.profile(this.world).ozone === 1,
            ozoneProfile: 'selected 10–40 km triangle; effective RGB coefficients',
            diffuse: 'selected sky/surface two-boundary closure; regional Lambertian albedo 0.18',
          },
          time_s: this.atm.uniforms.uColumnTime.value,
          renderCache: {
            skyPixels: this.skyTarget ? [this.skyTarget.width, this.skyTarget.height] : null,
            shadowPixels: this.shadowTarget
              ? [this.shadowTarget.width, this.shadowTarget.height]
              : null,
            quality: this.quality || 'economy',
            bakes: this.bakeCount,
            lastSubmission_ms: this.lastBakeMs,
            windCadence_s: 2,
            positionThreshold_m: 32,
            shadowPlaneAltitude_m: 0,
            shadowExtent_m: 32768,
          },
        }
      : null;
  }
  initRenderer(renderer, prelude) {
    const T = this.T;
    this.renderer = renderer;
    const options = {
      type: T.HalfFloatType,
      depthBuffer: false,
      minFilter: T.LinearFilter,
      magFilter: T.LinearFilter,
    };
    // The full sphere uses a signed-square-root elevation coordinate. Half of
    // the angular rows cover the visible upper sky; the lower rows support reflections.
    this.skyTarget = new T.WebGLRenderTarget(768, 384, options);
    this.skyTarget.texture.wrapS = T.RepeatWrapping;
    this.shadowTarget = new T.WebGLRenderTarget(128, 128, options);
    const u = this.atm.uniforms;
    u.uColumnSkyMap.value = this.skyTarget.texture;
    u.uColumnShadowMap.value = this.shadowTarget.texture;
    const vertex = 'varying vec2 vUV;void main(){vUV=uv;gl_Position=vec4(position.xy,0.,1.);}';
    this.bakeScene = new T.Scene();
    this.bakeCamera = new T.OrthographicCamera(-1, 1, 1, -1, 0, 1);
    this.skyBakeMaterial = new T.ShaderMaterial({
      uniforms: u,
      vertexShader: vertex,
      fragmentShader:
        prelude +
        DENSITY_GLSL +
        SKY_GLSL +
        `
   varying vec2 vUV;
   void main(){float q=vUV.y*2.-1.,e=sign(q)*q*q*1.57079632679,az=(vUV.x-.5)*6.283185307;
    vec3 d=vec3(cos(e)*cos(az),sin(e),cos(e)*sin(az));gl_FragColor=omColumnTransfer(d);}
  `,
      toneMapped: false,
      depthTest: false,
      depthWrite: false,
    });
    this.shadowBakeMaterial = new T.ShaderMaterial({
      uniforms: u,
      vertexShader: vertex,
      fragmentShader:
        prelude +
        DENSITY_GLSL +
        `
   varying vec2 vUV;
   void main(){vec2 p=uColumnShadowBounds.xy+vUV*uColumnShadowBounds.z;float sag=dot(p,p)/(uR+sqrt(max(1.,uR*uR-dot(p,p))));
    vec3 ro=vec3(p.x,uR-sag,p.y);float t=uSun.y<0.?1.:exp(-omColumnSunDepth(ro));gl_FragColor=vec4(t,t,t,1.);}
  `,
      toneMapped: false,
      depthTest: false,
      depthWrite: false,
    });
    this.bakeMesh = new T.Mesh(new T.PlaneGeometry(2, 2), this.skyBakeMaterial);
    this.bakeScene.add(this.bakeMesh);
  }
  setQuality(quality) {
    if (!['economy', 'balanced', 'high'].includes(quality))
      throw new RangeError('Unknown cloud quality');
    const width = quality === 'high' ? 3072 : quality === 'balanced' ? 1536 : 768;
    this.quality = quality;
    if (this.skyTarget && this.skyTarget.width !== width) {
      this.skyTarget.setSize(width, width / 2);
      this.renderKey = null;
      this.atm.clearEnvKey = null;
      this.atm.envTime = -1e9;
    }
  }

  prepare(renderer) {
    const u = this.atm.uniforms;
    if (this.key === 'reference' || !this.current || !this.skyTarget) return false;
    const time = Math.floor(u.uTime.value / 2) * 2,
      observer = u.uObserver.value;
    const key = [
      this.world,
      this.key,
      u.uBlend.value,
      u.uSide.value,
      ...u.uSun.value.toArray(),
      time,
      Math.floor(observer.x / 32),
      Math.floor(observer.y / 32),
      Math.round(u.uEyeHeight.value),
    ].join(':');
    if (key === this.renderKey) return false;
    const started = performance.now(),
      oldTarget = renderer.getRenderTarget(),
      tm = renderer.toneMapping,
      exposure = renderer.toneMappingExposure,
      oldTime = u.uColumnTime.value;
    const centreX = Math.round(observer.x / 512) * 512,
      centreZ = Math.round(observer.y / 512) * 512;
    u.uColumnShadowBounds.value.set(centreX - 16384, centreZ - 16384, 32768, 128);
    u.uColumnAverageShadow.value =
      1 -
      this.current.model.inputs.morphology.coverage *
        (1 -
          Math.exp(
            -this.current.model.summary.inCloudOpticalDepth / Math.max(0.1, u.uSun.value.y),
          ));
    u.uColumnTime.value = time;
    u.uColumnCacheReady.value = 0;
    u.uColumnShadowReady.value = 0;
    try {
      renderer.toneMapping = this.T.NoToneMapping;
      renderer.toneMappingExposure = 1;
      this.bakeMesh.material = this.skyBakeMaterial;
      // Tile expensive ray integration so higher angular resolution does not become
      // one unbounded GPU draw. Global fragment coordinates preserve the same jitter.
      const sky = this.skyTarget;
      sky.scissorTest = true;
      for (let y = 0; y < sky.height; y += 128)
        for (let x = 0; x < sky.width; x += 256) {
          sky.scissor.set(x, y, Math.min(256, sky.width - x), Math.min(128, sky.height - y));
          renderer.setRenderTarget(sky);
          renderer.render(this.bakeScene, this.bakeCamera);
        }
      sky.scissorTest = false;
      this.bakeMesh.material = this.shadowBakeMaterial;
      renderer.setRenderTarget(this.shadowTarget);
      renderer.render(this.bakeScene, this.bakeCamera);
      u.uColumnCacheReady.value = 1;
      u.uColumnShadowReady.value = 1;
      this.renderKey = key;
      this.bakeCount++;
    } finally {
      renderer.setRenderTarget(oldTarget);
      renderer.toneMapping = tm;
      renderer.toneMappingExposure = exposure;
      u.uColumnTime.value = oldTime;
      this.lastBakeMs = performance.now() - started;
    }
    return true;
  }

  dispose() {
    this.noiseTexture.dispose();
    this.profileTexture.dispose();
    this.sunTexture.dispose();
    this.cache.clear();
    this.skyTarget?.dispose();
    this.shadowTarget?.dispose();
    this.skyBakeMaterial?.dispose();
    this.shadowBakeMaterial?.dispose();
    this.bakeMesh?.geometry.dispose();
  }
}
const API = {
  noiseLattice,
  UNIFORMS,
  DENSITY_GLSL,
  SHADOW_GLSL,
  SKY_GLSL,
  CACHED_SKY_GLSL,
  Controller,
};
export default API;
