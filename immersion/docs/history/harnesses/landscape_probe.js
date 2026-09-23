/* Landscape GPU optics checks. Three submerged rays and one independently
 * raycast foreground-rock occlusion, using the actual rendered depth texture.
 * Flat-water path estimates retain the M1 diagnostic approximation. */
() => {
 const a=openMoonShoreline,T=THREE,r=a.renderer;
 const old=r.getRenderTarget(),oldTM=r.toneMapping,oldExp=r.toneMappingExposure;
 const rt=new T.WebGLRenderTarget(1,1,{type:T.FloatType,depthBuffer:false});
 const quadScene=new T.Scene(),quadCamera=new T.OrthographicCamera(-1,1,1,-1,0,1);
 const u={tex:{value:a.water.opaqueTarget.depthTexture},sampleUV:{value:new T.Vector2(.5,.44)},nearZ:{value:a.camera.near},farZ:{value:a.camera.far}};
 const shader=new T.ShaderMaterial({uniforms:u,vertexShader:'void main(){gl_Position=vec4(position.xy,0.,1.);}',fragmentShader:'uniform sampler2D tex;uniform vec2 sampleUV;uniform float nearZ,farZ;void main(){float d=texture2D(tex,sampleUV).r;float z=2.*nearZ*farZ/(farZ+nearZ-(2.*d-1.)*(farZ-nearZ));gl_FragColor=vec4(d,z,0.,1.);}',toneMapped:false});
 const q=new T.Mesh(new T.PlaneGeometry(2,2),shader);quadScene.add(q);
 const pixels=new Float32Array(4),depth=[];
 for(const uv of [[.50,.43],[.50,.45],[.50,.48],[.62,.43]]){
  u.sampleUV.value.set(...uv);r.setRenderTarget(rt);r.render(quadScene,quadCamera);r.readRenderTargetPixels(rt,0,0,1,1,pixels);
  const ray=new T.Raycaster();ray.setFromCamera(new T.Vector2(uv[0]*2-1,uv[1]*2-1),a.camera);
  const t=-a.camera.position.y/ray.ray.direction.y,hit=ray.ray.at(t,new T.Vector3());
  const waterZ=-hit.clone().applyMatrix4(a.camera.matrixWorldInverse).z;
  const rockHits=uv[0]>.6?ray.intersectObjects(a.scene.children.filter(o=>o.name.startsWith('substrate-rocks')),false):[];
  const blocker=rockHits.find(o=>o.distance<t);
  depth.push({uv,occluder:blocker?{name:blocker.object.name,distance_m:blocker.distance,flatWaterDistance_m:t}:null,deviceDepth:pixels[0],opaqueViewZ_m:pixels[1],flatWaterViewZ_m:waterZ,geometricDepth_m:OpenMoonCore.waterDepth(hit.x,hit.z),estimatedOpticalPath_m:Math.max(0,(pixels[1]-waterZ)*t/waterZ)});
 }
 const scene=new T.Scene();scene.environment=a.scene.environment;
 const cardMat=new T.MeshStandardMaterial({color:new T.Color(.18,.18,.18),roughness:1,metalness:0});
 const card=new T.Mesh(new T.PlaneGeometry(4,4),cardMat);card.rotation.x=-Math.PI/2;scene.add(card);
 const sun=a.scene.children.find(o=>o.isDirectionalLight),light=new T.DirectionalLight(sun.color,sun.intensity);light.position.set(0,10,0);scene.add(light);
 const camera=new T.OrthographicCamera(-1,1,1,-1,.1,20);camera.position.set(0,5,.0001);camera.lookAt(0,0,0);
 r.toneMapping=T.NoToneMapping;r.setRenderTarget(rt);r.render(scene,camera);r.readRenderTargetPixels(rt,0,0,1,1,pixels);
 const raw=Array.from(pixels).slice(0,3),gain=a.atmosphere.uniforms.uWhiteBalance.value.toArray();
 const toneUniforms={uRGB:{value:new T.Vector3(...raw)},uGain:{value:new T.Vector3(...gain)},toneMappingExposure:{value:.28}};
 q.material=new T.ShaderMaterial({uniforms:toneUniforms,vertexShader:shader.vertexShader,fragmentShader:'uniform vec3 uRGB,uGain;\n'+T.ShaderChunk.tonemapping_pars_fragment+'\nvoid main(){gl_FragColor=vec4(ACESFilmicToneMapping(uRGB*uGain),1.);}',toneMapped:false});
 function meter(e){toneUniforms.toneMappingExposure.value=e;r.render(quadScene,quadCamera);r.readRenderTargetPixels(rt,0,0,1,1,pixels);return .2126*pixels[0]+.7152*pixels[1]+.0722*pixels[2];}
 let lo=.001,hi=2;for(let i=0;i<22;i++){const e=(lo+hi)/2;if(meter(e)<.18)lo=e;else hi=e;}
 const calibrated=(lo+hi)/2,oldLuminance=meter(oldExp),calibratedLuminance=meter(calibrated);
 r.setRenderTarget(old);r.toneMapping=oldTM;r.toneMappingExposure=oldExp;
 const record={depth,grayCard:{reflectance:.18,rawLinearRGB:raw,displayWhiteBalanceGain:gain,referenceExposure:oldExp,referenceDisplayLuminance:oldLuminance,exposureFor18Percent:calibrated,calibratedDisplayLuminance:calibratedLuminance},renderTargets:{reflection:[a.water.reflectionTarget.width,a.water.reflectionTarget.height],samples:a.water.reflectionTarget.samples,opaque:[a.water.opaqueTarget.width,a.water.opaqueTarget.height]},inheritedNoonLux:a.atmosphere.clearLux};
 shader.dispose();q.material.dispose();q.geometry.dispose();rt.dispose();card.geometry.dispose();cardMat.dispose();
 return record;
}
