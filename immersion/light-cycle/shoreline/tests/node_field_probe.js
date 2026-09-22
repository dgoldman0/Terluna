/* Actual TSL render-target readback, independent of screenshot colour. */
async () => {
 const a=rendererLab,modules=await OM.loadThree(),T=modules.webgpu,N=modules.tsl,r=a.renderer;
 const old=r.getRenderTarget(),tm=r.toneMapping,exposure=r.toneMappingExposure;
 const rt=new T.RenderTarget(1,1,{type:T.FloatType,depthBuffer:false}),scene=new T.Scene(),camera=new T.OrthographicCamera(-1,1,1,-1,0,1);
 const f=a.geography.fieldUniforms,uv=N.uniform(new T.Vector2()),node=N.texture(f.uSurfaceState.value,uv),m=new T.NodeMaterial({toneMapped:false});m.fragmentNode=node;
 const g=new T.PlaneGeometry(2,2);scene.add(new T.Mesh(g,m));const records=[];
 try{
  r.toneMapping=T.NoToneMapping;r.toneMappingExposure=1;r.setRenderTarget(rt);
  for(const name of ['uSurfaceState','uPondLevels','uMaterialFields','uEnvironmentFields','uBathymetry']){
   const tex=f[name].value,{data,width,height}=tex.image;node.value=tex;
   for(const ch of name==='uSurfaceState'?[0,1,2,3]:[0]){
    let best=0;for(let i=0;i<width*height;i++)if(data[i*4+ch]>data[best*4+ch])best=i;
    uv.value.set((best%width+.5)/width,(Math.floor(best/width)+.5)/height);r.render(scene,camera);
    const gpu=Array.from(await r.readRenderTargetPixelsAsync(rt,0,0,1,1)),cpu=Array.from(data.slice(best*4,best*4+4));
    records.push({field:name,channel:ch,cell:best,cpu,gpu,maxAbsoluteError:Math.max(...cpu.map((v,i)=>Math.abs(v-gpu[i])))});
   }
  }
 }finally{r.setRenderTarget(old);r.toneMapping=tm;r.toneMappingExposure=exposure;g.dispose();m.dispose();rt.dispose();}
 return {records,maxAbsoluteError:Math.max(...records.map(x=>x.maxAbsoluteError)),actualBackend:a.report.backend};
}
