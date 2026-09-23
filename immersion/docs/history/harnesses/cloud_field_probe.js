/* Actual float-target readback of column extinction/phase/wind and solar transport. */
() => {
 const a=openMoonShoreline,T=THREE,r=a.renderer,c=a.atmosphere.columnClouds;
 if(!c.current)throw Error('Select a column study before the field probe');
 const old=r.getRenderTarget(),tm=r.toneMapping,exposure=r.toneMappingExposure;
 const target=new T.WebGLRenderTarget(1,1,{type:T.FloatType,depthBuffer:false});
 const scene=new T.Scene(),camera=new T.OrthographicCamera(-1,1,1,-1,0,1),u={field:{value:c.profileTexture},sampleUV:{value:new T.Vector2()}};
 const material=new T.ShaderMaterial({uniforms:u,toneMapped:false,vertexShader:'void main(){gl_Position=vec4(position.xy,0.,1.);}',fragmentShader:'uniform sampler2D field;uniform vec2 sampleUV;void main(){gl_FragColor=texture2D(field,sampleUV);}'});
 const geo=new T.PlaneGeometry(2,2);scene.add(new T.Mesh(geo,material));const rows=[],out=new Float32Array(4);let max=0,noiseMaterial=null;
 try{
  r.toneMapping=T.NoToneMapping;r.toneMappingExposure=1;r.setRenderTarget(target);
  for(const [name,tex]of [['column',c.profileTexture],['sunlight',c.sunTexture]]){
   u.field.value=tex;const {data,width,height}=tex.image;
   const indices=name==='column'?[0,64,128,192,255]:[32*128+64,32*128+80,63*128+110];
   for(const k of indices){u.sampleUV.value.set((k%width+.5)/width,(Math.floor(k/width)+.5)/height);r.render(scene,camera);r.readRenderTargetPixels(target,0,0,1,1,out);
    const expected=Array.from(data.slice(k*4,k*4+4)),actual=Array.from(out),error=Math.max(...actual.map((v,i)=>Math.abs(v-expected[i])));max=Math.max(max,error);rows.push({field:name,texel:k,expected,actual,error});}
  }
  // The cached sky and cloud-shadow outputs are sampled through the GPU too.
  const outputs=[];
  for(const [name,t]of [['sky',c.skyTarget],['cloudShadow',c.shadowTarget]]){
   u.field.value=t.texture;
   for(const uv of [[.12,.74],[.46,.81],[.75,.65]]){u.sampleUV.value.set(...uv);r.render(scene,camera);r.readRenderTargetPixels(target,0,0,1,1,out);const values=Array.from(out);if(!values.every(Number.isFinite))throw Error('Nonfinite '+name+' cache');if(name==='sky'&&(values[3]<0||values[3]>1))throw Error('Invalid cloud transmission');if(name==='cloudShadow'&&values.some(v=>v<0||v>1))throw Error('Invalid cloud shadow');outputs.push({field:name,uv,values});}
  }
  const noiseSamples=[],noiseU={noiseTex:{value:c.noiseTexture},noiseUV:{value:new T.Vector3()}};
  noiseMaterial=new T.ShaderMaterial({uniforms:noiseU,toneMapped:false,vertexShader:'void main(){gl_Position=vec4(position.xy,0.,1.);}',fragmentShader:'uniform highp sampler3D noiseTex;uniform vec3 noiseUV;void main(){float n=texture(noiseTex,noiseUV).r;gl_FragColor=vec4(n,n,n,1.);}'});
  scene.children[0].material=noiseMaterial;
  for(const xyz of [[0,0,0],[17,5,9],[39,48,21]]){noiseU.noiseUV.value.set(...xyz.map(v=>(v+.5)/64));r.render(scene,camera);r.readRenderTargetPixels(target,0,0,1,1,out);const expected=c.noiseTexture.image.data[(xyz[2]*64+xyz[1])*64+xyz[0]]/255,error=Math.abs(out[0]-expected);max=Math.max(max,error);noiseSamples.push({xyz,expected,actual:out[0],error});}
  return {maxAbsoluteError:max,samples:rows,cacheSamples:outputs,noiseSamples,backend:'WebGL2 float render-target readback'};
 }finally{r.setRenderTarget(old);r.toneMapping=tm;r.toneMappingExposure=exposure;geo.dispose();material.dispose();noiseMaterial?.dispose();target.dispose();}
}
