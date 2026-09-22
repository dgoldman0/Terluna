(function(root){
'use strict';const OM=root.OM,C=root.OpenMoonCore;
function makeTexture(T,kind,seed){const size=512,c=document.createElement('canvas');c.width=c.height=size;const ctx=c.getContext('2d'),img=ctx.createImageData(size,size),rnd=C.rng(seed);for(let y=0;y<size;y++)for(let x=0;x<size;x++){const i=(x+y*size)*4,n=C.fbm(x/40,y/(kind==='bark'?9:40)),grain=rnd();let v=.5+n*.46+(grain-.5)*.16;if(kind==='bark')v*=.7+.3*Math.pow(Math.abs(Math.sin(x*.12+C.noise(x*.06,y*.018)*5)),.6);for(let k=0;k<3;k++)img.data[i+k]=C.clamp(v)*255;img.data[i+3]=255;}ctx.putImageData(img,0,0);const t=new T.CanvasTexture(c);t.wrapS=t.wrapT=T.RepeatWrapping;t.colorSpace=T.SRGBColorSpace;t.anisotropy=8;t.repeat.set(kind==='bark'?2:40,kind==='bark'?3:40);return t;}
function mesh(T,geo,mat,group,pos,scale){const o=new T.Mesh(geo,mat);if(pos)o.position.set(...pos);if(scale)o.scale.set(...scale);o.castShadow=o.receiveShadow=true;group.add(o);return o;}
function cylinderBetween(T,a,b,r0,r1,mat,group){const p=new T.Vector3(...a),q=new T.Vector3(...b),v=q.clone().sub(p);const o=mesh(T,new T.CylinderGeometry(r1,r0,v.length(),7),mat,group);o.position.copy(p.add(q).multiplyScalar(.5));o.quaternion.setFromUnitVectors(new T.Vector3(0,1,0),v.normalize());return o;}
function prepareMaterial(T,mat,atm,state,{sway=0,sheltered=false,ground=false}={}){
 const uniforms={...atm.uniforms,uWet:{value:0},uCanopyWet:{value:0},uLeafWet:{value:0},uSway:{value:sway},uProtected:{value:sheltered?1:0},uLeafSurface:{value:sway?1:0}};
 mat.userData.omUniforms=uniforms;state.materials.push(mat);
 const coverDecl=ground?'varying float vOMCanopy;\n':'';
 const extra=coverDecl+`varying vec3 vOMWorld;uniform float uWet,uCanopyWet,uLeafWet,uSway,uProtected,uLeafSurface;\n`;
 mat.onBeforeCompile=shader=>{
  Object.assign(shader.uniforms,uniforms);
  shader.vertexShader=(ground?'attribute float omCanopy;varying float vOMCanopy;\n':'')+'varying vec3 vOMWorld;uniform float uTime,uWind,uSway;\n'+shader.vertexShader;
  shader.vertexShader=shader.vertexShader.replace('#include <begin_vertex>',`#include <begin_vertex>
   ${ground?'vOMCanopy=omCanopy;':''}
   vec4 omP=vec4(position,1.);\n#ifdef USE_INSTANCING\nomP=instanceMatrix*omP;\n#endif
   omP=modelMatrix*omP;
   transformed.x+=uSway*sin(omP.x*.36+omP.z*.16+uTime*(.65+uWind*.035))*min(1.5,uWind*.18);
   transformed.z+=uSway*.45*sin(omP.x*.18-omP.z*.28+uTime*.77)*min(1.5,uWind*.18);`);
  shader.vertexShader=shader.vertexShader.replace('#include <project_vertex>',`vec4 omW=vec4(transformed,1.);\n#ifdef USE_INSTANCING\nomW=instanceMatrix*omW;\n#endif
   vOMWorld=(modelMatrix*omW).xyz;\n#include <project_vertex>`);
  shader.fragmentShader=extra+OM.SKY_UNIFORMS+OM.GLSL_NOISE+OM.CLOUDS+OM.CLEAR_LOOKUP+shader.fragmentShader;
  shader.fragmentShader=shader.fragmentShader.replace('#include <roughnessmap_fragment>',`#include <roughnessmap_fragment>
   float omRoof=(1.-step(5.7,abs(vOMWorld.x-20.)))*(1.-step(4.2,abs(vOMWorld.z-33.)))*(1.-step(8.0,vOMWorld.y));
   float omWet=mix(${ground?'mix(uWet,uCanopyWet,vOMCanopy)':'mix(uWet,uLeafWet,uLeafSurface)'},0.,max(omRoof,uProtected));
   roughnessFactor=mix(roughnessFactor,max(.12,roughnessFactor*.34),omWet);
   diffuseColor.rgb*=mix(1.,.70,omWet);`);
  const lightChunk=T.ShaderChunk.lights_fragment_begin;
  const needle='getDirectionalLightInfo( directionalLight, directLight );';
  if(!lightChunk.includes(needle))throw Error('Pinned Three.js directional-light shader hook has changed.');
  shader.fragmentShader=shader.fragmentShader.replace('#include <lights_fragment_begin>',lightChunk.replace(needle,needle+'\n directLight.color *= omCloudShadow(vOMWorld);'));
  shader.fragmentShader=shader.fragmentShader.replace('#include <lights_fragment_end>',`#include <lights_fragment_end>
   reflectedLight.indirectDiffuse*=mix(1.,.24,max(omRoof,uProtected*.35))*${ground?'mix(1.,.58,vOMCanopy)':'1.'};
   reflectedLight.indirectSpecular*=mix(1.,.35,max(omRoof,uProtected*.35));`);
  shader.fragmentShader=shader.fragmentShader.replace('#include <tonemapping_fragment>',`vec3 omView=vOMWorld-cameraPosition;vec3 omFogTr=exp(-uLocalExtinction*length(omView));gl_FragColor.rgb=mix(omClear(normalize(omView)),gl_FragColor.rgb,omFogTr);\n#ifdef TONE_MAPPING\n gl_FragColor.rgb*=uWhiteBalance;\n#endif\n#include <tonemapping_fragment>`);
 };
 mat.customProgramCacheKey=()=>`open-moon-material-${sway}-${sheltered}-${ground}`;
 return mat;
}
function animateDepth(T,o,uniforms,sway){if(!sway)return;const d=new T.MeshDepthMaterial({depthPacking:T.RGBADepthPacking,side:T.DoubleSide});d.onBeforeCompile=s=>{Object.assign(s.uniforms,uniforms);s.vertexShader='uniform float uTime,uWind,uSway;\n'+s.vertexShader;s.vertexShader=s.vertexShader.replace('#include <begin_vertex>',`#include <begin_vertex>
 vec4 p=vec4(position,1.);\n#ifdef USE_INSTANCING\np=instanceMatrix*p;\n#endif
 p=modelMatrix*p;transformed.x+=uSway*sin(p.x*.36+p.z*.16+uTime*(.65+uWind*.035))*min(1.5,uWind*.18);transformed.z+=uSway*.45*sin(p.x*.18-p.z*.28+uTime*.77)*min(1.5,uWind*.18);`);};o.customDepthMaterial=d;}
function makeLeaf(T){const verts=[],uv=[],inds=[];const rows=5;for(let i=0;i<rows;i++){const t=i/(rows-1),w=Math.sin(Math.PI*t)*.105;verts.push(-w,t*.34,Math.sin(t*Math.PI)*.022,w,t*.34,Math.sin(t*Math.PI)*.022);uv.push(0,t,1,t);if(i<rows-1){const j=i*2;inds.push(j,j+1,j+2,j+1,j+3,j+2);}}const g=new T.BufferGeometry();g.setAttribute('position',new T.Float32BufferAttribute(verts,3));g.setAttribute('uv',new T.Float32BufferAttribute(uv,2));g.setIndex(inds);g.computeVertexNormals();return g;}
function batchInstances(T,source,scene,cell=24){
 const groups=new Map(),m=new T.Matrix4(),c=new T.Color();
 for(let i=0;i<source.count;i++){source.getMatrixAt(i,m);const key=Math.floor(m.elements[12]/cell)+':'+Math.floor(m.elements[14]/cell);if(!groups.has(key))groups.set(key,[]);groups.get(key).push(i);}
 for(const [key,ids]of groups){const group=new T.InstancedMesh(source.geometry,source.material,ids.length);ids.forEach((id,i)=>{source.getMatrixAt(id,m);group.setMatrixAt(i,m);if(source.instanceColor){source.getColorAt(id,c);group.setColorAt(i,c);}});group.castShadow=source.castShadow;group.receiveShadow=source.receiveShadow;group.customDepthMaterial=source.customDepthMaterial;group.name='spatial-batch-'+key;group.computeBoundingSphere();scene.add(group);}
}
function terrainGeometry(T,world){
 const angular=384,radii=[];
 for(const [a,b,step]of [[.5,64,.5],[66,256,2],[264,1024,8],[1056,4096,32],[4224,16384,128]])for(let r=a;r<=b;r+=step)radii.push(r);
 const pos=[0,C.groundHeight(0,0,world),0],norm=[0,1,0],index=[];
 for(let j=0;j<radii.length;j++){
  const r=radii[j];
  for(let i=0;i<angular;i++){
   const a=i*C.TAU/angular,x=Math.cos(a)*r,z=Math.sin(a)*r;
   pos.push(x,C.groundHeight(x,z,world),z);
   const h=Math.max(.15,r*.0008),g=C.surfaceGradient(x,z,h),R=C.worldRadius(world);
   const nx=-g.x+x/R,nz=-g.z+z/R,len=Math.hypot(nx,1,nz);norm.push(nx/len,1/len,nz/len);
   const b=1+j*angular+i,c=1+j*angular+(i+1)%angular;
   if(j===0)index.push(0,c,b);
   else {const a=b-angular,d=c-angular;index.push(a,d,b,b,d,c);}
  }
 }
 const geo=new T.BufferGeometry();geo.setAttribute('position',new T.Float32BufferAttribute(pos,3));geo.setAttribute('normal',new T.Float32BufferAttribute(norm,3));geo.setIndex(index);
 geo.setAttribute('omCanopy',new T.Float32BufferAttribute(new Float32Array(pos.length/3),1));
 geo.computeBoundingSphere();geo.userData={radialRings:radii.length,angular,outerRadius:16384,continuous:true};return geo;
}
function createScene(T,scene,atm,world='moon'){
 const state={materials:[],obstacles:[],trees:[],puddles:[],groundMeshes:[],animated:[],world},obj=new T.Object3D();
 const groundMat=OM.material(T,atm,state,'terrain');
 const ground=mesh(T,terrainGeometry(T,world),groundMat,scene);ground.castShadow=false;ground.name='continuous-terrain-and-bathymetry';state.groundMeshes.push(ground);
 const rockMat=OM.material(T,atm,state,'rock',{roughness:.88});
 const rockFamilies=[];
 for(let k=0;k<8;k++){
  const g=new T.SphereGeometry(1,64,40),p=g.attributes.position,shape=C.rng(88+k*371);
  const planes=[[1,0,0,.72+shape()*.26],[-1,0,0,.70+shape()*.25],[0,1,0,.66+shape()*.20],[0,-1,0,.52],[0,0,1,.74+shape()*.27],[0,0,-1,.78+shape()*.18]];
  for(let j=0;j<7;j++){const n=new T.Vector3(shape()*2-1,shape()*1.7-.6,shape()*2-1).normalize();planes.push([...n.toArray(),.76+shape()*.32]);}
  for(let i=0;i<p.count;i++){
   const x=p.getX(i),y=p.getY(i),z=p.getZ(i),cuts=[];
   for(const [nx,ny,nz,d]of planes){const dot=x*nx+y*ny+z*nz;if(dot>.001)cuts.push(d/dot);}
   const near=Math.min(...cuts);let blend=0;for(const t of cuts)blend+=Math.exp(-26*(t-near));
   const radius=near-Math.log(blend)/26;
   const erosion=(C.fbm(x*11+k*17,z*11+y*9)-.45)*.025;
   p.setXYZ(i,x*(radius+erosion),y*(radius+erosion),z*(radius+erosion));
  }g.computeVertexNormals();rockFamilies.push(g);
 }
 const rnd=C.rng(50724),rocks=[];
 for(const [cx,cz,n,span,scale]of [[-9,-5,18,5,1.3],[15,-7,20,6,1.2],[-30,-24,18,11,2.4],[39,-36,22,15,1.7],[-5,16,16,12,.5]]){
  for(let i=0;i<n;i++){const x=cx+(rnd()-.5)*span*2,z=cz+(rnd()-.5)*span*2;
   if(Math.abs(x)<3&&z<12&&z>-9)continue;
   const s=scale*(.25+Math.pow(rnd(),1.6)*1.5);rocks.push({x,z,s,y:C.groundHeight(x,z,world),rx:rnd()*.25,ry:rnd()*C.TAU,rz:rnd()*.3,k:Math.floor(rnd()*8)});
  }
 }
 for(let k=0;k<8;k++){
  const group=rocks.filter(r=>r.k===k),inst=new T.InstancedMesh(rockFamilies[k],rockMat,group.length);
  group.forEach((r,i)=>{obj.position.set(r.x,r.y+.20*r.s,r.z);obj.scale.set(r.s*(1+rnd()*.6),r.s,r.s*(.8+rnd()*.6));obj.rotation.set(r.rx,r.ry,r.rz);obj.updateMatrix();inst.setMatrixAt(i,obj.matrix);if(r.s>.6)state.obstacles.push({x:r.x,z:r.z,r:r.s*.65});});
  inst.castShadow=inst.receiveShadow=true;inst.computeBoundingSphere();scene.add(inst);
 }
 const pebbles=[],pr=C.rng(734),col=new T.Color();
 for(let i=0;i<16000;i++){
  const x=(pr()-.5)*100,z=-12+pr()*54,h=C.surfaceHeight(x,z);
  const patch=C.noise(x*.19+4,z*.19),band=Math.exp(-(((h-.60)/.8)**2));
  if(h<-.8||h>2.8||pr()>(.08+.70*band)*(.25+patch*.7))continue;
  const s=.014+pr()**4*.12;pebbles.push([x,z,s]);
 }
 const pebMat=new T.MeshStandardMaterial({color:0xffffff,roughness:.92});OM.prepareMaterial(T,pebMat,atm,state);
 const pebGeo=new T.IcosahedronGeometry(1,1),peb=new T.InstancedMesh(pebGeo,pebMat,pebbles.length);
 pebbles.forEach(([x,z,s],i)=>{obj.position.set(x,C.groundHeight(x,z,world)+s*.28,z);obj.scale.set(s*(1+pr()),s*.65,s);obj.rotation.set(pr(),pr()*C.TAU,pr());obj.updateMatrix();peb.setMatrixAt(i,obj.matrix);const v=.045+pr()*.19;col.setRGB(v,v*.95,v*.8);peb.setColorAt(i,col);});
 peb.castShadow=true;peb.receiveShadow=true;peb.computeBoundingSphere();batchInstances(T,peb,scene,12);
 const grassMat=OM.material(T,atm,state,'grass',{color:0x78884b,side:T.DoubleSide,roughness:.88});
 const gr=C.rng(9472),blades=[];
 for(let c=0;c<2400;c++){
  const cx=(gr()-.5)*135,cz=-4+gr()*82,h=C.surfaceHeight(cx,cz);
  if(h<.8||h>8||C.pathDistance(cx,cz)<1.7||C.roofMask(cx,cz,1))continue;
  const density=C.noise(cx*.11+4,cz*.11);
  if(gr()>C.smooth(.3,.7,density)*C.smooth(.9,2.5,h))continue;
  const clump=.10+gr()*.33;
  for(let j=0;j<26;j++){
   const theta=gr()*C.TAU,rr=Math.sqrt(gr())*clump,x=cx+Math.cos(theta)*rr,z=cz+Math.sin(theta)*rr;
   blades.push([x,z,.13+gr()*.31,gr()*C.TAU,.75+gr()*.7]);
  }
 }
 const bp=[],bi=[],bn=[];
 for(let j=0;j<=5;j++){const t=j/5,w=.0065*(1-t)+.0002;bp.push(-w,t,.22*t*t,w,t,.22*t*t);if(j<5){const i=j*2;bi.push(i,i+1,i+2,i+1,i+3,i+2);}}
 const bg=new T.BufferGeometry();bg.setAttribute('position',new T.Float32BufferAttribute(bp,3));bg.setIndex(bi);bg.computeVertexNormals();
 const grass=new T.InstancedMesh(bg,grassMat,blades.length);
 blades.forEach(([x,z,h,r,s],i)=>{obj.position.set(x,C.groundHeight(x,z,world),z);obj.scale.set(s,h,1);obj.rotation.set(0,r,0);obj.updateMatrix();grass.setMatrixAt(i,obj.matrix);const v=gr();col.setRGB(.65+v*.35,.74+v*.23,.42+v*.32);grass.setColorAt(i,col);});
 grass.receiveShadow=true;grass.computeBoundingSphere();animateDepth(T,grass,grassMat.userData.omUniforms,.008);batchInstances(T,grass,scene,18);
 const wood=OM.material(T,atm,state,'wood'),leafMat=OM.material(T,atm,state,'leaf',{color:0x63823e,side:T.DoubleSide,roughness:.8});
 const tr=C.rng(3242),segments=[],leaves=[];
 const sites=[[-10,-1,7.6],[25,3,6.7],[-24,23,10.5],[33,27,9.2]];
 for(let i=0;i<25;i++){const x=(tr()-.5)*150,z=24+tr()*100;if(C.pathDistance(x,z)<5||C.roofMask(x,z,3))continue;sites.push([x,z,7+tr()*6]);}
 function branch(a,dir,len,r,depth){
  const b=[a[0]+dir[0]*len,a[1]+dir[1]*len,a[2]+dir[2]*len];segments.push([a,b,r,r*.61]);
  if(depth>0){for(let j=0;j<3;j++){const angle=tr()*C.TAU,v=new T.Vector3(dir[0]*.4+Math.cos(angle)*.65,.30+tr()*.6,dir[2]*.4+Math.sin(angle)*.65).normalize();branch(b,v.toArray(),len*(.58+tr()*.14),r*.5,depth-1);}}
  else for(let j=0;j<48;j++){const u=tr()*C.TAU,c=tr()*2-1,rad=Math.cbrt(tr()),st=Math.sqrt(1-c*c);leaves.push([b[0]+Math.cos(u)*st*rad*.85,b[1]+c*rad*.48,b[2]+Math.sin(u)*st*rad*.85,.8+tr()*.7,tr()*C.TAU,tr()*C.TAU,tr()]);}
 }
 for(const [x,z,h]of sites){const y=C.groundHeight(x,z,world);if(y<.3)continue;state.trees.push({x,z,y,height:h});state.obstacles.push({x,z,r:.35});const a=[x,y,z],b=[x+.5,y+h*.45,z+.15];segments.push([a,b,h*.029,h*.018]);
  for(let j=0;j<5;j++){const theta=j*2.399+tr()*.4,v=new T.Vector3(Math.cos(theta)*.5,.50+tr()*.35,Math.sin(theta)*.5).normalize();branch([b[0],b[1]+j*.4,b[2]],v.toArray(),h*.21,h*.011,2);}
 }
 const tube=new T.CylinderGeometry(.64,1,1,10),trunks=new T.InstancedMesh(tube,wood,segments.length),up=new T.Vector3(0,1,0);
 segments.forEach(([a,b,r0,r1],i)=>{const v=new T.Vector3(b[0]-a[0],b[1]-a[1],b[2]-a[2]);obj.position.set((a[0]+b[0])*.5,(a[1]+b[1])*.5,(a[2]+b[2])*.5);obj.quaternion.setFromUnitVectors(up,v.clone().normalize());obj.scale.set(r0,v.length(),r0);obj.updateMatrix();trunks.setMatrixAt(i,obj.matrix);});trunks.castShadow=trunks.receiveShadow=true;trunks.computeBoundingSphere();batchInstances(T,trunks,scene,24);
 const leaf=new T.InstancedMesh(makeLeaf(T),leafMat,leaves.length);
 leaves.forEach(([x,y,z,s,rx,ry,v],i)=>{obj.position.set(x,y,z);obj.rotation.set(rx,ry,rx*.7);obj.scale.setScalar(s);obj.updateMatrix();leaf.setMatrixAt(i,obj.matrix);col.setRGB(.72+v*.23,.77+v*.21,.60+v*.22);leaf.setColorAt(i,col);});leaf.castShadow=leaf.receiveShadow=true;leaf.computeBoundingSphere();animateDepth(T,leaf,leafMat.userData.omUniforms,.012);batchInstances(T,leaf,scene,18);
 const S=C.SHELTER,deckY=C.groundHeight(S.x,S.z,world)+.22;state.pavilionY=deckY;
 const deckMat=OM.material(T,atm,state,'wood'),frameMat=new T.MeshStandardMaterial({color:0x313e39,metalness:.3,roughness:.5}),roofMat=new T.MeshStandardMaterial({color:0x616960,metalness:.35,roughness:.45});
 OM.prepareMaterial(T,frameMat,atm,state);OM.prepareMaterial(T,roofMat,atm,state);
 mesh(T,new T.BoxGeometry(11,.24,8),deckMat,scene,[S.x,deckY-.12,S.z]);
 for(const dx of[-4.8,4.8])for(const dz of[-3.3,3.3]){cylinderBetween(T,[S.x+dx,deckY,S.z+dz],[S.x+dx,deckY+4.1,S.z+dz],.115,.115,frameMat,scene);state.obstacles.push({x:S.x+dx,z:S.z+dz,r:.25});}
 const roof=mesh(T,new T.BoxGeometry(12,.16,9.2),roofMat,scene,[S.x,deckY+4.2,S.z]);roof.rotation.x=.045;
 const lamp=new T.PointLight(0xffd29b,0,22,2);lamp.position.set(S.x,deckY+3.5,S.z);scene.add(lamp);state.lamp=lamp;
 state.bulb=mesh(T,new T.CylinderGeometry(.14,.14,.06,16),new T.MeshStandardMaterial({color:0xffdfa7,emissive:0xffdfa7,emissiveIntensity:0}),scene,[S.x,deckY+3.92,S.z]);
 const ramp=new T.PlaneGeometry(4,7,1,12);ramp.rotateX(-Math.PI/2);const rp=ramp.attributes.position;
 for(let i=0;i<rp.count;i++){const x=rp.getX(i)+S.x,z=rp.getZ(i)+S.z-7.5;rp.setXYZ(i,x,C.mix(C.groundHeight(x,z,world),deckY,C.smooth(0,1,(z-S.z+11)/7)),z);}ramp.computeVertexNormals();mesh(T,ramp,deckMat,scene);
 state.height=(x,z)=>{let y=C.groundHeight(x,z,state.world);if(C.roofMask(x,z))y=Math.max(y,deckY);if(Math.abs(x-S.x)<2&&z>S.z-11&&z<S.z-4)y=C.mix(y,deckY,C.smooth(0,1,(z-S.z+11)/7));return y;};
 state.canMove=(x,z)=>{if(Math.abs(x)>130||z>140||C.surfaceHeight(x,z)<.12)return false;return !state.obstacles.some(o=>Math.hypot(x-o.x,z-o.z)<o.r+.28);};
 state.update=l=>{for(const m of state.materials){const u=m.userData.omUniforms;if(u){u.uWet.value=C.clamp(l.exposed/.9);u.uCanopyWet.value=C.clamp(l.canopy/.7);u.uLeafWet.value=C.clamp(l.leaf/.35);}}};
 state.counts={trees:state.trees.length,leaves:leaves.length,branchSegments:segments.length,grassBlades:blades.length,rocks:rocks.length,pebbles:pebbles.length,terrainVertices:ground.geometry.attributes.position.count};
 return state;
}
OM.createScene=createScene;OM.prepareMaterial=prepareMaterial;OM.terrainGeometry=terrainGeometry;
})(globalThis);
