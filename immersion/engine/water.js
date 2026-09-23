import { OM } from './om.js';
import C from './core.js';
/* M1 water: a continuous curved surface, resolved planar near-field reflection,
 * depth-buffer refraction and RGB Beer–Lambert attenuation of the actual seabed.
 * The selected wave modes are a render model, with lunar gravity in dispersion.
 */
function waterGeometry(T) {
  const verts = [],
    indices = [],
    angular = 256,
    radii = [];
  for (const [a, b, step] of [
    [0.75, 96, 0.75],
    [98, 256, 2],
    [272, 1024, 16],
    [1088, 4096, 64],
    [4352, 16384, 256],
  ])
    for (let r = a; r <= b; r += step) radii.push(r);
  verts.push(0, 0, 0);
  radii.forEach((r, j) => {
    for (let i = 0; i < angular; i++) {
      const a = (i * C.TAU) / angular;
      verts.push(Math.cos(a) * r, 0, Math.sin(a) * r);
      const b = 1 + j * angular + i,
        c = 1 + j * angular + ((i + 1) % angular);
      if (j === 0) indices.push(0, c, b);
      else {
        const a = b - angular,
          d = c - angular;
        indices.push(a, d, b, b, d, c);
      }
    }
  });
  const geo = new T.BufferGeometry();
  geo.setAttribute('position', new T.Float32BufferAttribute(verts, 3));
  geo.setIndex(indices);
  return geo;
}
function createWater(T, scene, atm, renderer, fields) {
  const reflectionTarget = new T.WebGLRenderTarget(1280, 800, {
    type: T.HalfFloatType,
    depthBuffer: true,
    generateMipmaps: true,
    minFilter: T.LinearMipmapLinearFilter,
  });
  reflectionTarget.samples = OM.graphicsBackend?.software ? 0 : 4;
  reflectionTarget.texture.colorSpace = T.LinearSRGBColorSpace;
  const opaqueTarget = new T.WebGLRenderTarget(1280, 800, {
    type: T.HalfFloatType,
    depthBuffer: true,
  });
  opaqueTarget.samples = OM.graphicsBackend?.software ? 0 : 4;
  opaqueTarget.texture.colorSpace = T.LinearSRGBColorSpace;
  opaqueTarget.depthTexture = new T.DepthTexture(1280, 800, T.UnsignedIntType);
  const u = {
    ...atm.uniforms,
    ...fields,
    uReflection: { value: reflectionTarget.texture },
    uOpaque: { value: opaqueTarget.texture },
    uOpaqueDepth: { value: opaqueTarget.depthTexture },
    uMirrorMatrix: { value: new T.Matrix4() },
    uGravity: { value: 1.62 },
    uCamera: { value: new T.Vector3() },
    uRain: { value: 0 },
    uViewport: { value: new T.Vector2(1280, 800) },
    uNear: { value: 0.08 },
    uFar: { value: 30000 },
    uReflectValid: { value: 0 },
  };
  const wave = `
 uniform float uTime,uWind,uGravity,uR;

 uniform sampler2D uBathymetry;uniform vec4 uFieldBounds;
 float omEnvelope(vec2 p){
  vec2 q=(p-uFieldBounds.xy)/uFieldBounds.z,uv=(q+.5)/uFieldBounds.w;
  vec2 a=smoothstep(vec2(0),vec2(6),q),b=1.-smoothstep(vec2(uFieldBounds.w-7.),vec2(uFieldBounds.w-1.),q);
  float cover=a.x*a.y*b.x*b.y;
  float h=mix(12.,max(0.,-texture2D(uBathymetry,uv).r),cover);
  float amplitude=0.;
  ${C.WATER_BANDS.map((k) => `amplitude+=.095*pow(.095/${k.toFixed(6)},1.12)*(.4+uWind*.18);`).join('')}
  return h/(h+amplitude/.45);
 }
 vec3 omWaves(vec2 p,bool micro){
  vec3 sum=vec3(0.);float filterWidth=0.;
  float ks[12];${C.WATER_BANDS.map((k, i) => `ks[${i}]=${k.toFixed(6)};`).join('')}
  for(int i=0;i<12;i++){
   if(!micro&&i>5)break;
   float k=ks[i],theta=.95+sin(float(i)*2.399)*.71;vec2 dir=vec2(cos(theta),sin(theta));
   float a=.095*pow(.095/k,1.12)*(.4+uWind*.18);
   float depth=12.,e=exp(-2.*k*depth),omega=sqrt(uGravity*k*(1.-e)/(1.+e));
   float phase=k*dot(p,dir)-omega*uTime+float(i)*2.721;
   sum.x+=a*sin(phase);sum.yz-=a*k*dir*cos(phase);
  }
  float envelope=omEnvelope(p);
  vec2 gradient=vec2(omEnvelope(p+vec2(.001,0))-omEnvelope(p-vec2(.001,0)),omEnvelope(p+vec2(0,.001))-omEnvelope(p-vec2(0,.001)))/.002;
  sum.yz=sum.yz*envelope-sum.x*gradient;sum.x*=envelope;
  return sum;
 }
 `;
  const vs =
    wave +
    `
 uniform mat4 uMirrorMatrix;varying vec3 vW;varying vec4 vMirror;varying float vViewZ;
 void main(){vec3 p=position;vec3 w=omWaves(p.xz,false);float r2=dot(p.xz,p.xz);p.y=w.x-r2/(uR+sqrt(max(1.,uR*uR-r2)));
 vW=p;vMirror=uMirrorMatrix*vec4(p,1.);vec4 view=viewMatrix*vec4(p,1.);vViewZ=-view.z;gl_Position=projectionMatrix*view;}
 `;
  const fs =
    OM.SKY_UNIFORMS +
    OM.GLSL_NOISE +
    OM.CLOUDS +
    OM.CLEAR_LOOKUP +
    wave.replace('uniform float uTime,uWind,uGravity,uR;', 'uniform float uGravity;') +
    `
 uniform sampler2D uReflection,uOpaque,uOpaqueDepth;uniform vec3 uCamera;uniform vec2 uViewport;uniform float uRain,uNear,uFar,uReflectValid;
 varying vec3 vW;varying vec4 vMirror;varying float vViewZ;
 float eyeZ(float d){return 2.*uNear*uFar/(uFar+uNear-(2.*d-1.)*(uFar-uNear));}
 void main(){
  vec3 wave=omWaves(vW.xz,true),v=normalize(uCamera-vW);
  float footprint=max(length(dFdx(vW.xz)),length(dFdy(vW.xz)));
  float ripple=(1.-smoothstep(.08,.6,footprint));
  vec2 fine=vec2(sin(vW.x*38.+vW.z*17.+uTime*.7),sin(vW.z*41.-vW.x*13.+uTime*.5))*.015*ripple;
  vec3 n=normalize(vec3(wave.y+fine.x+vW.x/uR,1.,wave.z+fine.y+vW.z/uR));
  float nv=max(.001,dot(n,v)),fres=.0204+.9796*pow(1.-nv,5.);
  vec2 screen=gl_FragCoord.xy/uViewport,uv=screen+n.xz*.010;
  float opaqueZ=eyeZ(texture2D(uOpaqueDepth,uv).r);
  if(opaqueZ<vViewZ+.005){uv=screen;opaqueZ=eyeZ(texture2D(uOpaqueDepth,uv).r);}
  float path=clamp((opaqueZ-vViewZ)*length(uCamera-vW)/max(.05,vViewZ),0.,150.);
  vec3 extinction=vec3(.29,.071,.036),transmission=exp(-extinction*path);
  vec3 irradiance=(uDiffuse+uDirect*max(0.,uSun.y)*omCloudShadow(vW))/OM_PI;
  vec3 waterScatter=irradiance*vec3(.009,.037,.042);
  vec3 bottom=texture2D(uOpaque,uv).rgb;
  vec3 body=bottom*transmission+waterScatter*(1.-transmission);
  vec3 reflectedDir=reflect(-v,n);vec3 reflected=omClear(normalize(reflectedDir));
  vec2 mirrorUV=vMirror.xy/vMirror.w+n.xz*.006;
  float margin=min(min(mirrorUV.x,mirrorUV.y),min(1.-mirrorUV.x,1.-mirrorUV.y));
  float planarWeight=smoothstep(.002,.02,margin)*(1.-smoothstep(550.,1800.,length(uCamera-vW)))*uReflectValid;
  reflected=mix(reflected,texture2D(uReflection,clamp(mirrorUV,.001,.999),.65).rgb,planarWeight);
  vec3 colour=mix(body,reflected,fres);
  vec3 h=normalize(v+uSun);float nl=max(0.,dot(n,uSun)),nh=max(0.,dot(n,h)),vh=max(0.,dot(v,h));
  float rough=.042+uWind*.006,alpha2=max(pow(rough,4.),.00465421*.00465421);
  float den=nh*nh*(alpha2-1.)+1.,D=alpha2/(OM_PI*den*den);
  float k=pow(rough+1.,2.)/8.,G=(nv/(nv*(1.-k)+k))*(nl/(nl*(1.-k)+k));
  float F=.0204+.9796*pow(1.-vh,5.);colour+=uDirect*D*G*F/(4.*nv+.0001)*omCloudShadow(vW);
  float contact=(1.-smoothstep(.02,.09,path))*smoothstep(.0,.012,path);
  colour=mix(colour,irradiance*.18,contact*.14);
  vec3 fogT=exp(-uLocalExtinction*length(uCamera-vW));
  colour=mix(omAirColour(normalize(vW-uCamera)),colour,fogT);
  gl_FragColor=vec4(max(colour,vec3(0.)),1.);
  #ifdef TONE_MAPPING\n gl_FragColor.rgb*=uWhiteBalance;\n#endif\n#include <tonemapping_fragment>
  #include <colorspace_fragment>
 }
 `;
  const geo = waterGeometry(T);
  const mat = new T.ShaderMaterial({
    uniforms: u,
    vertexShader: vs,
    fragmentShader: fs,
    toneMapped: true,
  });
  const water = new T.Mesh(geo, mat);
  water.frustumCulled = false;
  water.name = 'refractive-water';
  scene.add(water);
  const mirrorCam = new T.PerspectiveCamera(),
    bias = new T.Matrix4().set(0.5, 0, 0, 0.5, 0, 0.5, 0, 0.5, 0, 0, 0.5, 0.5, 0, 0, 0, 1),
    dir = new T.Vector3();
  let width = 0,
    height = 0,
    quality = 'balanced',
    passes = 0;
  function resize(w, h, q = quality) {
    quality = q;
    if (
      w === width &&
      h === height &&
      reflectionTarget.width === Math.round(w * (q === 'economy' ? 0.5 : 1))
    )
      return;
    width = w;
    height = h;
    const scale = q === 'economy' ? 0.5 : 1;
    reflectionTarget.setSize(
      Math.max(256, Math.round(w * scale)),
      Math.max(256, Math.round(h * scale)),
    );
    opaqueTarget.setSize(w, h);
    u.uViewport.value.set(w, h);
  }
  function reflection(camera, now, force = false) {
    const size = renderer.getDrawingBufferSize(new T.Vector2());
    resize(size.x, size.y);
    u.uCamera.value.copy(camera.position);
    u.uNear.value = camera.near;
    u.uFar.value = camera.far;
    mirrorCam.copy(camera);
    mirrorCam.position.y = -camera.position.y;
    camera.getWorldDirection(dir);
    dir.y = -dir.y;
    mirrorCam.up.set(0, -1, 0);
    mirrorCam.lookAt(mirrorCam.position.clone().add(dir));
    mirrorCam.updateMatrixWorld();
    u.uMirrorMatrix.value
      .copy(bias)
      .multiply(mirrorCam.projectionMatrix)
      .multiply(mirrorCam.matrixWorldInverse);
    const old = {
      rt: renderer.getRenderTarget(),
      clip: renderer.clippingPlanes,
      tm: renderer.toneMapping,
      exposure: renderer.toneMappingExposure,
      shadow: renderer.shadowMap.autoUpdate,
      disk: atm.uniforms.uShowDisk.value,
    };
    water.visible = false;
    renderer.toneMapping = T.NoToneMapping;
    renderer.toneMappingExposure = 1;
    try {
      renderer.clippingPlanes = [];
      renderer.setRenderTarget(opaqueTarget);
      renderer.render(scene, camera);
      renderer.shadowMap.autoUpdate = false;
      renderer.clippingPlanes = [new T.Plane(new T.Vector3(0, 1, 0), -0.015)];
      atm.uniforms.uShowDisk.value = false;
      renderer.setRenderTarget(reflectionTarget);
      renderer.render(scene, mirrorCam);
      u.uReflectValid.value = 1;
      passes++;
    } finally {
      renderer.setRenderTarget(old.rt);
      renderer.clippingPlanes = old.clip;
      renderer.toneMapping = old.tm;
      renderer.toneMappingExposure = old.exposure;
      renderer.shadowMap.autoUpdate = old.shadow;
      atm.uniforms.uShowDisk.value = old.disk;
      water.visible = true;
    }
  }
  return {
    mesh: water,
    uniforms: u,
    reflection,
    resize,
    reflectionTarget,
    opaqueTarget,
    get passes() {
      return passes;
    },
  };
}
function createRain(T, scene, atm) {
  const random = C.rng(719),
    positions = [],
    seeds = [];
  for (let i = 0; i < 5000; i++) {
    const x = random(),
      y = random(),
      z = random();
    for (let j = 0; j < 2; j++) {
      positions.push(x, y, z);
      seeds.push(j);
    }
  }
  const g = new T.BufferGeometry();
  g.setAttribute('position', new T.Float32BufferAttribute(positions, 3));
  g.setAttribute('endPoint', new T.Float32BufferAttribute(seeds, 1));
  const u = {
    ...atm.uniforms,
    uCamera: { value: new T.Vector3() },
    uFall: { value: 1.56 },
    uRain: { value: 0 },
    uRoofY: { value: 8.24 },
  };
  const vs = `uniform float uTime,uWind,uFall,uRain,uRoofY,uDrift;uniform vec3 uCamera;attribute float endPoint;varying float vAlpha;void main(){vec3 p=position;float x=fract(p.x+uDrift*.018),z=fract(p.z+uDrift*.005);float y=mod(p.y*24.-uTime*uFall,24.);vec3 w=vec3((x-.5)*50.+uCamera.x,y+uCamera.y-10.,(z-.5)*50.+uCamera.z);w-=endPoint*vec3(uWind*.018,-uFall/35.,uWind*.005);float protectedArea=(1.-step(6.,abs(w.x-20.)))*(1.-step(4.6,abs(w.z-33.)))*(1.-step(uRoofY,w.y));vAlpha=step(position.x,clamp(uRain/10.,0.,1.))*(1.-protectedArea)*smoothstep(-1.,1.,w.y)*.30;gl_Position=projectionMatrix*viewMatrix*vec4(w,1.);}`;
  const fs = `uniform vec3 uDiffuse,uDirect,uWhiteBalance;uniform vec3 uSun;varying float vAlpha;void main(){vec3 c=(uDiffuse+uDirect*max(0.,uSun.y)*.1)*.18;gl_FragColor=vec4(c,vAlpha);\n#ifdef TONE_MAPPING\n gl_FragColor.rgb*=uWhiteBalance;\n#endif\n#include <tonemapping_fragment>\n#include <colorspace_fragment>\n}`;
  const m = new T.ShaderMaterial({
    uniforms: u,
    vertexShader: vs,
    fragmentShader: fs,
    transparent: true,
    depthWrite: false,
  });
  const rain = new T.LineSegments(g, m);
  rain.frustumCulled = false;
  scene.add(rain);
  return { mesh: rain, uniforms: u };
}
OM.waterGeometry = waterGeometry;
OM.createWater = createWater;
OM.createRain = createRain;
