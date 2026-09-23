/* Read back actual uploaded state texels; no screenshot colour inference. */
() => {
 const a=openMoonShoreline,T=THREE,r=a.renderer,w=a.geography.surfaceWater;
 const old=r.getRenderTarget(),tm=r.toneMapping,exp=r.toneMappingExposure;
 const rt=new T.WebGLRenderTarget(1,1,{type:T.FloatType,depthBuffer:false});
 const scene=new T.Scene(),camera=new T.OrthographicCamera(-1,1,1,-1,0,1);
 const u={uStateTex:{value:a.geography.fieldUniforms.uSurfaceState.value},uSampleUV:{value:new T.Vector2()}};
 const m=new T.ShaderMaterial({uniforms:u,vertexShader:'void main(){gl_Position=vec4(position.xy,0.,1.);}',fragmentShader:'uniform sampler2D uStateTex;uniform vec2 uSampleUV;void main(){gl_FragColor=texture2D(uStateTex,uSampleUV);}',toneMapped:false});
 const g=new T.PlaneGeometry(2,2);scene.add(new T.Mesh(g,m));
 const bytes=w.textureData(),out=new Float32Array(4),samples=[];let maximumAbsoluteError=0;
 try {
  r.toneMapping=T.NoToneMapping;r.toneMappingExposure=1;r.setRenderTarget(rt);
  for(let channel=0;channel<4;channel++){
   let best=0;for(let i=0;i<w.n*w.n;i++)if(bytes[i*4+channel]>bytes[best*4+channel])best=i;
   u.uSampleUV.value.set((best%w.n+.5)/w.n,(Math.floor(best/w.n)+.5)/w.n);r.render(scene,camera);r.readRenderTargetPixels(rt,0,0,1,1,out);
   const expected=Array.from(bytes.slice(best*4,best*4+4)),actual=Array.from(out);
   expected.forEach((v,i)=>maximumAbsoluteError=Math.max(maximumAbsoluteError,Math.abs(v-actual[i])));
   samples.push({channel,cell:best,expected,actual});
  }
 } finally {r.setRenderTarget(old);r.toneMapping=tm;r.toneMappingExposure=exp;g.dispose();m.dispose();rt.dispose();}
 return {samples,maximumAbsoluteError};
}
