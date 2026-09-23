/* Run in the exact revision-04 viewer. Read its production density GLSL on GPU,
 * then feed one immutable voxel field to both the GPU and offline reference.
 * The benchmark isolates cloud transport with unit normal solar irradiance,
 * black external boundaries and zero gas; it is not the full scene comparison. */
async ({regime='fair',world='moon',dimensions=[49,33,49]})=>{
 const app=openMoonShoreline,T=THREE,renderer=app.renderer,C=OpenMoonCloudRenderer,F=OpenMoonFrozenCloud;
 const ctrl=app.atmosphere.columnClouds;ctrl.select(regime);ctrl.update(world);
 const u=app.atmosphere.uniforms;u.uR.value=ctrl.current.model.planet.radius;u.uColumnTime.value=0;
 const base=u.uColumnBase.value,top=u.uColumnTop.value,extent=Math.max(2000,Math.min(50000,u.uColumnScale.value));
 const lo=[-extent,base,-extent],hi=[extent,top,extent],N=dimensions.reduce((a,b)=>a*b,1),width=512,height=Math.ceil(N/width);
 const old=renderer.getRenderTarget(),tm=renderer.toneMapping,exposure=renderer.toneMappingExposure;
 const scene=new T.Scene(),camera=new T.OrthographicCamera(-1,1,1,-1,0,1),geometry=new T.PlaneGeometry(2,2),mesh=new T.Mesh(geometry);scene.add(mesh);
 const vertex='void main(){gl_Position=vec4(position.xy,0.,1.);}';
 const target=new T.WebGLRenderTarget(width,height,{type:T.FloatType,depthBuffer:false});
 const uniforms={...u,a1BoundsMin:{value:new T.Vector3(...lo)},a1BoundsMax:{value:new T.Vector3(...hi)},a1Dims:{value:new T.Vector3(...dimensions)},a1Width:{value:width}};
 const production=new T.ShaderMaterial({uniforms,toneMapped:false,vertexShader:vertex,fragmentShader:OM.SKY_UNIFORMS+OM.GLSL_NOISE+OM.CLOUDS+C.DENSITY_GLSL+`
 uniform vec3 a1BoundsMin,a1BoundsMax,a1Dims;uniform float a1Width;
 void main(){float k=floor(gl_FragCoord.y)*a1Width+floor(gl_FragCoord.x);float x=mod(k,a1Dims.x),y=mod(floor(k/a1Dims.x),a1Dims.y),z=floor(k/(a1Dims.x*a1Dims.y));
 vec3 p=mix(a1BoundsMin,a1BoundsMax,vec3(x,y,z)/(a1Dims-1.));vec3 ro=vec3(p.x,uR+p.y,p.z);float h=omRadialHeight(ro);gl_FragColor=vec4(omColumnExt(ro),omColumnAt(h).g,0.,1.);}`});
 const resources=[target,production,geometry],results={};
 try{
  renderer.toneMapping=T.NoToneMapping;renderer.toneMappingExposure=1;mesh.material=production;renderer.setRenderTarget(target);renderer.render(scene,camera);
  const raw=new Float32Array(width*height*4);renderer.readRenderTargetPixels(target,0,0,width,height,raw);
  const ext=new Float32Array(N),ice=new Float32Array(N);for(let i=0;i<N;i++){ext[i]=raw[i*4];ice[i]=Math.max(0,Math.min(1,raw[i*4+1]));}
  const base64=a=>{const bytes=new Uint8Array(a.buffer);let s='';for(let i=0;i<bytes.length;i+=32768)s+=String.fromCharCode(...bytes.subarray(i,i+32768));return btoa(s);};
  const spec={schema:'open-moon-frozen-cloud/1',dimensions,bounds_min_m:lo,bounds_max_m:hi,extinction_float32le:base64(ext),ice_fraction_float32le:base64(ice),single_scattering_albedo:1,
   provenance:{regime,world,time_s:0,source:'Actual revision-04 production density GLSL float-target readback',radial_geometry:'p_world=(x,R+y,z); exported bounds are Cartesian, not radial heights',column:ctrl.current.model.summary},
   phase:'0.85 HG(0.78+0.04*iceFraction) + 0.15 HG(-0.25)',boundary:'Finite crop; density is zero outside the box. This truncation is a separate scene assumption.'};
  const field=F.create(spec);results.field=spec;
  const texture=data=>{const t=new T.Data3DTexture(data,...dimensions);t.format=T.RedFormat;t.type=T.FloatType;t.minFilter=t.magFilter=T.LinearFilter;t.unpackAlignment=1;t.generateMipmaps=false;t.needsUpdate=true;resources.push(t);return t;};
  const fu={a1FrozenExt:{value:texture(ext)},a1FrozenIce:{value:texture(ice)},a1FrozenMin:{value:new T.Vector3(...lo)},a1FrozenMax:{value:new T.Vector3(...hi)},a1FrozenSize:{value:new T.Vector3(...dimensions)}};
  const sampleTarget=new T.WebGLRenderTarget(1,1,{type:T.FloatType,depthBuffer:false});resources.push(sampleTarget);
  const sampleU={...fu,point:{value:new T.Vector3()}},material=new T.ShaderMaterial({uniforms:sampleU,toneMapped:false,vertexShader:vertex,fragmentShader:F.GLSL+'uniform vec3 point;void main(){vec2 v=a1FrozenSample(point);gl_FragColor=vec4(v,0.,1.);}'});resources.push(material);mesh.material=material;renderer.setRenderTarget(sampleTarget);
  const pixel=new Float32Array(4);let maxFieldError=0;const samples=[];
  for(let i=0;i<24;i++){const p=[.5+.43*Math.sin(i*1.73),.5+.47*Math.sin(i*2.37+.3),.5+.45*Math.cos(i*1.17)].map((v,k)=>lo[k]+v*(hi[k]-lo[k]));sampleU.point.value.set(...p);renderer.render(scene,camera);renderer.readRenderTargetPixels(sampleTarget,0,0,1,1,pixel);const cpu=field.sample(p);const error=Math.abs(pixel[0]-cpu.extinction);maxFieldError=Math.max(maxFieldError,error);samples.push({position:p,gpu:pixel[0],cpu:cpu.extinction,error});}
  results.field_parity={maximum_absolute_extinction_error_m1:maxFieldError,samples};
  // Off-grid comparisons quantify the separate discretization of the original
  // continuous density field into this solver-neutral voxel volume.
  const pu={...u,point:{value:new T.Vector3()}},continuous=new T.ShaderMaterial({uniforms:pu,toneMapped:false,vertexShader:vertex,fragmentShader:OM.SKY_UNIFORMS+OM.GLSL_NOISE+OM.CLOUDS+C.DENSITY_GLSL+'uniform vec3 point;void main(){vec3 p=vec3(point.x,uR+point.y,point.z);gl_FragColor=vec4(omColumnExt(p),0.,0.,1.);}'});resources.push(continuous);mesh.material=continuous;
  let maxVoxel=0,absoluteSum=0,sourceSum=0;const voxelSamples=[];
  for(const row of samples){pu.point.value.set(...row.position);renderer.render(scene,camera);renderer.readRenderTargetPixels(sampleTarget,0,0,1,1,pixel);const error=Math.abs(pixel[0]-row.cpu);maxVoxel=Math.max(maxVoxel,error);absoluteSum+=error;sourceSum+=Math.abs(pixel[0]);voxelSamples.push({position:row.position,continuous_production:pixel[0],frozen:row.cpu,error});}
  results.voxelization={maximum_absolute_extinction_error_m1:maxVoxel,sum_absolute_difference_over_sum_continuous:sourceSum?absoluteSum/sourceSum:null,samples:voxelSamples,scope:'24 deterministic off-grid probes. These are representation errors, separate from solver agreement; no whole-field bound.'};

  const rays=[];
  for(let iz=0;iz<3;iz++)for(let ix=0;ix<4;ix++)rays.push({origin:[lo[0]+(ix+.5)/4*(hi[0]-lo[0]),lo[1]-100,lo[2]+(iz+.5)/3*(hi[2]-lo[2])],direction:[0,1,0]});
  const sun=[Math.SQRT1_2,Math.SQRT1_2,0];
  const ru={...fu,origin:{value:new T.Vector3()},direction:{value:new T.Vector3()},sun:{value:new T.Vector3(...sun)},primary:{value:128},secondary:{value:32}};
  const transport=new T.ShaderMaterial({uniforms:ru,toneMapped:false,vertexShader:vertex,fragmentShader:F.GLSL+`
 uniform vec3 origin,direction,sun;uniform float primary,secondary;
 float sunTr(vec3 p){vec2 b=a1Box(p,sun);if(b.y<=b.x)return 1.;float ds=(b.y-b.x)/secondary,tau=0.;for(int j=0;j<1024;j++){if(float(j)>=secondary)break;tau+=a1FrozenSample(p+sun*(b.x+(float(j)+.5)*ds)).x*ds;}return exp(-tau);}
 void main(){vec2 b=a1Box(origin,direction);float tr=1.,radiance=0.;if(b.y>b.x){float ds=(b.y-b.x)/primary;for(int i=0;i<4096;i++){if(float(i)>=primary)break;vec3 p=origin+direction*(b.x+(float(i)+.5)*ds);vec2 s=a1FrozenSample(p);float a=1.-exp(-s.x*ds);radiance+=tr*a*a1Phase(dot(direction,sun),s.y)*sunTr(p);tr*=1.-a;}}gl_FragColor=vec4(radiance,tr,0.,1.);}`});resources.push(transport);mesh.material=transport;
  results.transport={model:'Gray cloud-only first order; same frozen voxel field; unit collimated solar normal irradiance; black boundaries; no gas, surface, tone mapping or cache',sun_direction:sun,rays,readbacks:[]};
  for(const budget of[[128,32],[512,128],[2048,512]]){ru.primary.value=budget[0];ru.secondary.value=budget[1];const values=[];for(const ray of rays){ru.origin.value.set(...ray.origin);ru.direction.value.set(...ray.direction);renderer.render(scene,camera);renderer.readRenderTargetPixels(sampleTarget,0,0,1,1,pixel);if(!Array.from(pixel).every(Number.isFinite))throw Error('Nonfinite cloud transport');values.push({single_scatter_radiance:pixel[0],transmittance:pixel[1]});}results.transport.readbacks.push({primary:budget[0],secondary:budget[1],values});}
  const gl=renderer.getContext(),e=gl.getExtension('WEBGL_debug_renderer_info');if(gl.isContextLost())throw Error('Context lost');results.gpu={renderer:e?gl.getParameter(e.UNMASKED_RENDERER_WEBGL):gl.getParameter(gl.RENDERER),context_lost:false};
  return results;
 }finally{renderer.setRenderTarget(old);renderer.toneMapping=tm;renderer.toneMappingExposure=exposure;for(const v of resources)v.dispose();}
}
