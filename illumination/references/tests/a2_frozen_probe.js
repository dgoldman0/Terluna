/* A2 cloud field convergence probe. It reads the exact production density GLSL
 * into a solver-neutral voxel field, then measures discretization and finite-crop
 * effects separately. The production renderer itself is unchanged. */
async ({regime='fair',world='moon',dimensions=[49,33,49],extentMultiplier=1,probeExtentMultiplier=1,probeCount=96})=>{
 const app=openMoonShoreline,T=THREE,renderer=app.renderer,C=OpenMoonCloudRenderer,F=OpenMoonFrozenCloud;
 const ctrl=app.atmosphere.columnClouds;ctrl.select(regime);ctrl.update(world);
 const u=app.atmosphere.uniforms;u.uR.value=ctrl.current.model.planet.radius;u.uColumnTime.value=0;
 const base=u.uColumnBase.value,top=u.uColumnTop.value,baseExtent=Math.max(2000,Math.min(50000,u.uColumnScale.value)),extent=baseExtent*extentMultiplier;
 if(!dimensions.every(Number.isInteger)||dimensions.some(v=>v<2)||!(extentMultiplier>0)||!(probeExtentMultiplier>0&&probeExtentMultiplier<=extentMultiplier)||!Number.isInteger(probeCount)||probeCount<8)throw Error('Invalid A2 convergence configuration');
 const lo=[-extent,base,-extent],hi=[extent,top,extent],N=dimensions.reduce((a,b)=>a*b,1),width=512,height=Math.ceil(N/width);
 const old=renderer.getRenderTarget(),tm=renderer.toneMapping,exposure=renderer.toneMappingExposure;
 const scene=new T.Scene(),camera=new T.OrthographicCamera(-1,1,1,-1,0,1),geometry=new T.PlaneGeometry(2,2),mesh=new T.Mesh(geometry);scene.add(mesh);
 const vertex='void main(){gl_Position=vec4(position.xy,0.,1.);}';
 const target=new T.WebGLRenderTarget(width,height,{type:T.FloatType,depthBuffer:false});
 const uniforms={...u,a2BoundsMin:{value:new T.Vector3(...lo)},a2BoundsMax:{value:new T.Vector3(...hi)},a2Dims:{value:new T.Vector3(...dimensions)},a2Width:{value:width}};
 const production=new T.ShaderMaterial({uniforms,toneMapped:false,vertexShader:vertex,fragmentShader:OM.SKY_UNIFORMS+OM.GLSL_NOISE+OM.CLOUDS+C.DENSITY_GLSL+`
 uniform vec3 a2BoundsMin,a2BoundsMax,a2Dims;uniform float a2Width;
 void main(){float k=floor(gl_FragCoord.y)*a2Width+floor(gl_FragCoord.x);float x=mod(k,a2Dims.x),y=mod(floor(k/a2Dims.x),a2Dims.y),z=floor(k/(a2Dims.x*a2Dims.y));
 vec3 p=mix(a2BoundsMin,a2BoundsMax,vec3(x,y,z)/(a2Dims-1.));vec3 ro=vec3(p.x,uR+p.y,p.z);float h=omRadialHeight(ro);gl_FragColor=vec4(omColumnExt(ro),omColumnAt(h).g,0.,1.);}`});
 const resources=[target,production,geometry],results={};
 function halton(index,base){let f=1,r=0,i=index;while(i>0){f/=base;r+=f*(i%base);i=Math.floor(i/base);}return r;}
 try{
  renderer.toneMapping=T.NoToneMapping;renderer.toneMappingExposure=1;mesh.material=production;renderer.setRenderTarget(target);renderer.render(scene,camera);
  const raw=new Float32Array(width*height*4);renderer.readRenderTargetPixels(target,0,0,width,height,raw);
  const ext=new Float32Array(N),ice=new Float32Array(N);for(let i=0;i<N;i++){ext[i]=raw[i*4];ice[i]=Math.max(0,Math.min(1,raw[i*4+1]));}
  const base64=a=>{const bytes=new Uint8Array(a.buffer);let s='';for(let i=0;i<bytes.length;i+=32768)s+=String.fromCharCode(...bytes.subarray(i,i+32768));return btoa(s);};
  const spec={schema:'open-moon-frozen-cloud/1',dimensions,bounds_min_m:lo,bounds_max_m:hi,extinction_float32le:base64(ext),ice_fraction_float32le:base64(ice),single_scattering_albedo:1,
   provenance:{regime,world,time_s:0,source:'Actual revision-04 production density GLSL float-target readback',radial_geometry:'p_world=(x,R+y,z); exported bounds are Cartesian, not radial heights',column:ctrl.current.model.summary,a2:{baseExtent_m:baseExtent,extentMultiplier,probeExtentMultiplier,probeCount}},
   phase:'0.85 HG(0.78+0.04*iceFraction) + 0.15 HG(-0.25)',boundary:'Finite crop; density is zero outside the box. Crop convergence is measured independently.'};
  const field=F.create(spec);results.field=spec;
  const texture=data=>{const t=new T.Data3DTexture(data,...dimensions);t.format=T.RedFormat;t.type=T.FloatType;t.minFilter=t.magFilter=T.LinearFilter;t.unpackAlignment=1;t.generateMipmaps=false;t.needsUpdate=true;resources.push(t);return t;};
  const fu={a1FrozenExt:{value:texture(ext)},a1FrozenIce:{value:texture(ice)},a1FrozenMin:{value:new T.Vector3(...lo)},a1FrozenMax:{value:new T.Vector3(...hi)},a1FrozenSize:{value:new T.Vector3(...dimensions)}};
  const sampleTarget=new T.WebGLRenderTarget(1,1,{type:T.FloatType,depthBuffer:false});resources.push(sampleTarget);
  const sampleU={...fu,point:{value:new T.Vector3()}},material=new T.ShaderMaterial({uniforms:sampleU,toneMapped:false,vertexShader:vertex,fragmentShader:F.GLSL+'uniform vec3 point;void main(){vec2 v=a1FrozenSample(point);gl_FragColor=vec4(v,0.,1.);}'});resources.push(material);mesh.material=material;renderer.setRenderTarget(sampleTarget);
  const pu={...u,point:{value:new T.Vector3()}},continuous=new T.ShaderMaterial({uniforms:pu,toneMapped:false,vertexShader:vertex,fragmentShader:OM.SKY_UNIFORMS+OM.GLSL_NOISE+OM.CLOUDS+C.DENSITY_GLSL+'uniform vec3 point;void main(){vec3 p=vec3(point.x,uR+point.y,point.z);gl_FragColor=vec4(omColumnExt(p),0.,0.,1.);}'});resources.push(continuous);
  const pixel=new Float32Array(4);let maxFieldError=0,maxVoxel=0,absoluteSum=0,sourceSum=0,squaredSum=0;const samples=[];
  const probeExtent=baseExtent*probeExtentMultiplier;
  for(let i=1;i<=probeCount;i++){
   const p=[(2*halton(i,2)-1)*probeExtent,base+(top-base)*halton(i,3),(2*halton(i,5)-1)*probeExtent];
   mesh.material=material;sampleU.point.value.set(...p);renderer.render(scene,camera);renderer.readRenderTargetPixels(sampleTarget,0,0,1,1,pixel);const cpu=field.sample(p),gpu=pixel[0],fieldError=Math.abs(gpu-cpu.extinction);maxFieldError=Math.max(maxFieldError,fieldError);
   mesh.material=continuous;pu.point.value.set(...p);renderer.render(scene,camera);renderer.readRenderTargetPixels(sampleTarget,0,0,1,1,pixel);const source=pixel[0],error=Math.abs(source-cpu.extinction);maxVoxel=Math.max(maxVoxel,error);absoluteSum+=error;sourceSum+=Math.abs(source);squaredSum+=error*error;
   samples.push({position:p,continuous_production:source,frozen:cpu.extinction,gpu_frozen:gpu,field_sampling_error:fieldError,voxel_error:error});
  }
  results.field_parity={maximum_absolute_extinction_error_m1:maxFieldError};
  results.voxelization={maximum_absolute_extinction_error_m1:maxVoxel,sum_absolute_difference_over_sum_continuous:sourceSum?absoluteSum/sourceSum:null,rms_absolute_extinction_error_m1:Math.sqrt(squaredSum/probeCount),probe_count:probeCount,probe_extent_m:probeExtent,samples};

  // Common central rays are defined from the unscaled production extent, so crop
  // comparisons change only the finite transport boundary, not the viewed points.
  const rays=[];for(let iz=0;iz<3;iz++)for(let ix=0;ix<4;ix++)rays.push({origin:[-baseExtent+(ix+.5)/4*(2*baseExtent),base-100,-baseExtent+(iz+.5)/3*(2*baseExtent)],direction:[0,1,0]});
  const sun=[Math.SQRT1_2,Math.SQRT1_2,0],ru={...fu,origin:{value:new T.Vector3()},direction:{value:new T.Vector3()},sun:{value:new T.Vector3(...sun)},primary:{value:1024},secondary:{value:256}};
  const transport=new T.ShaderMaterial({uniforms:ru,toneMapped:false,vertexShader:vertex,fragmentShader:F.GLSL+`
 uniform vec3 origin,direction,sun;uniform float primary,secondary;
 float sunTr(vec3 p){vec2 b=a1Box(p,sun);if(b.y<=b.x)return 1.;float ds=(b.y-b.x)/secondary,tau=0.;for(int j=0;j<1024;j++){if(float(j)>=secondary)break;tau+=a1FrozenSample(p+sun*(b.x+(float(j)+.5)*ds)).x*ds;}return exp(-tau);}
 void main(){vec2 b=a1Box(origin,direction);float tr=1.,radiance=0.;if(b.y>b.x){float ds=(b.y-b.x)/primary;for(int i=0;i<4096;i++){if(float(i)>=primary)break;vec3 p=origin+direction*(b.x+(float(i)+.5)*ds);vec2 s=a1FrozenSample(p);float a=1.-exp(-s.x*ds);radiance+=tr*a*a1Phase(dot(direction,sun),s.y)*sunTr(p);tr*=1.-a;}}gl_FragColor=vec4(radiance,tr,0.,1.);}`});resources.push(transport);mesh.material=transport;
  const values=[];for(const ray of rays){ru.origin.value.set(...ray.origin);ru.direction.value.set(...ray.direction);renderer.render(scene,camera);renderer.readRenderTargetPixels(sampleTarget,0,0,1,1,pixel);if(!Array.from(pixel).every(Number.isFinite))throw Error('Nonfinite cloud transport');values.push({single_scatter_radiance:pixel[0],transmittance:pixel[1]});}
  results.transport={model:'Gray cloud-only first order convergence probe; same production-derived field; unit collimated Sun at 45 degrees; black boundaries; zero gas',sun_direction:sun,primary:1024,secondary:256,rays,values};
  const gl=renderer.getContext(),e=gl.getExtension('WEBGL_debug_renderer_info');if(gl.isContextLost())throw Error('Context lost');results.gpu={renderer:e?gl.getParameter(e.UNMASKED_RENDERER_WEBGL):gl.getParameter(gl.RENDERER),context_lost:false};
  return results;
 }finally{renderer.setRenderTarget(old);renderer.toneMapping=tm;renderer.toneMappingExposure=exposure;for(const v of resources)v.dispose();}
}
