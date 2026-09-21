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
  shader.fragmentShader=extra+OM.SKY_UNIFORMS+OM.GLSL_NOISE+OM.CLOUDS+shader.fragmentShader;
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
  shader.fragmentShader=shader.fragmentShader.replace('#include <tonemapping_fragment>',`vec3 omFog=(uDiffuse+uDirect*max(0.,uSun.y)*.14)/OM_PI;float omFogTr=exp(-3.912/max(80.,uVisibility)*length(cameraPosition-vOMWorld));gl_FragColor.rgb=mix(omFog,gl_FragColor.rgb,omFogTr);\n#include <tonemapping_fragment>`);
 };
 mat.customProgramCacheKey=()=>`open-moon-material-${sway}-${sheltered}-${ground}`;
 return mat;
}
function animateDepth(T,o,uniforms,sway){if(!sway)return;const d=new T.MeshDepthMaterial({depthPacking:T.RGBADepthPacking,side:T.DoubleSide});d.onBeforeCompile=s=>{Object.assign(s.uniforms,uniforms);s.vertexShader='uniform float uTime,uWind,uSway;\n'+s.vertexShader;s.vertexShader=s.vertexShader.replace('#include <begin_vertex>',`#include <begin_vertex>
 vec4 p=vec4(position,1.);\n#ifdef USE_INSTANCING\np=instanceMatrix*p;\n#endif
 p=modelMatrix*p;transformed.x+=uSway*sin(p.x*.36+p.z*.16+uTime*(.65+uWind*.035))*min(1.5,uWind*.18);transformed.z+=uSway*.45*sin(p.x*.18-p.z*.28+uTime*.77)*min(1.5,uWind*.18);`);};o.customDepthMaterial=d;}
function makeLeaf(T){const verts=[],uv=[],inds=[];const rows=5;for(let i=0;i<rows;i++){const t=i/(rows-1),w=Math.sin(Math.PI*t)*.105;verts.push(-w,t*.34,Math.sin(t*Math.PI)*.022,w,t*.34,Math.sin(t*Math.PI)*.022);uv.push(0,t,1,t);if(i<rows-1){const j=i*2;inds.push(j,j+1,j+2,j+1,j+3,j+2);}}const g=new T.BufferGeometry();g.setAttribute('position',new T.Float32BufferAttribute(verts,3));g.setAttribute('uv',new T.Float32BufferAttribute(uv,2));g.setIndex(inds);g.computeVertexNormals();return g;}
function createScene(T,scene,atm,world='moon'){
 const obj=new T.Object3D(),state={materials:[],obstacles:[],trees:[],puddles:[],groundMeshes:[],animated:[],pavilionY:0};
 const texture=makeTexture(T,'ground',84),bark=makeTexture(T,'bark',71);const bump=texture.clone();bump.colorSpace=T.NoColorSpace;const barkBump=bark.clone();barkBump.colorSpace=T.NoColorSpace;
 const terrainMat=prepareMaterial(T,new T.MeshStandardMaterial({color:0xffffff,vertexColors:true,map:texture,bumpMap:bump,bumpScale:.13,roughness:.95}),atm,state,{ground:true});
 function terrain(w,h,nx,nz,cx,cz,distant=false){const g=new T.PlaneGeometry(w,h,nx,nz);g.rotateX(-Math.PI/2);const p=g.attributes.position,colors=[];for(let i=0;i<p.count;i++){const x=p.getX(i)+cx,z=p.getZ(i)+cz;let y;if(!distant)y=C.groundHeight(x,z,world);else{const xx=x*.0014,zz=z*.0014;const ridges=45+170*Math.pow(C.fbm(xx*2,zz*2),2)+30*Math.sin(x*.0014);y=ridges*(.45+.55*C.smooth(-.3,.5,Math.sin(x*.0025+.4)))-(x*x+z*z)/(2*(world==='earth'?6371000:1737400));}p.setXYZ(i,x,y,z);const d=z-C.shore(x),path=C.pathDistance(x,z);const col=new T.Color();if(distant){const v=C.fbm(x*.023,z*.023);col.setRGB(.08+v*.11,.12+v*.12,.10+v*.12);}else if(d<7){const v=.82+.18*C.noise(x*2,z*2);col.setRGB(.61*v,.55*v,.41*v);}else{const f=C.fbm(x*.09,z*.09),trail=1-C.smooth(1.1,2.5,path);col.setRGB(C.mix(.105+f*.1,.30,trail),C.mix(.15+f*.13,.23,trail),C.mix(.055+f*.055,.145,trail));}colors.push(col.r,col.g,col.b);}g.setAttribute('color',new T.Float32BufferAttribute(colors,3));g.computeVertexNormals();const o=mesh(T,g,terrainMat,scene);o.castShadow=false;state.groundMeshes.push(o);return o;}
 terrain(360,300,220,180,0,85);terrain(5000,1600,180,50,300,-2200,true);
 // Rock shapes are shared, while proportions, erosion noise and colours vary.
 const rockMat=prepareMaterial(T,new T.MeshStandardMaterial({color:0x7c8274,map:texture,bumpMap:bump,bumpScale:.3,roughness:.93}),atm,state);
 const rockGeo=new T.IcosahedronGeometry(1,2);for(let i=0;i<rockGeo.attributes.position.count;i++){const p=rockGeo.attributes.position,x=p.getX(i),y=p.getY(i),z=p.getZ(i),f=.8+.25*C.noise(x*4+7,z*4+y*6);p.setXYZ(i,x*f,y*f,z*f);}rockGeo.computeVertexNormals();
 const rr=C.rng(1953),rocks=new T.InstancedMesh(rockGeo,rockMat,90);for(let i=0;i<90;i++){const x=(rr()-.5)*240,z=C.shore(x)+rr()*75,y=C.groundHeight(x,z,world),s=.4+rr()**3*4;obj.position.set(x,y+s*.3,z);obj.scale.set(s*(1+rr()),s*.65,s);obj.rotation.set(rr()*.4,rr()*C.TAU,rr()*.2);obj.updateMatrix();rocks.setMatrixAt(i,obj.matrix);if(s>1)state.obstacles.push({x,z,r:s*.8});}rocks.castShadow=rocks.receiveShadow=true;rocks.computeBoundingSphere();scene.add(rocks);
 const woodMat=prepareMaterial(T,new T.MeshStandardMaterial({color:0x78664d,map:bark,bumpMap:barkBump,bumpScale:.08,roughness:.93}),atm,state);
 const leafMat=prepareMaterial(T,new T.MeshStandardMaterial({color:0x789447,roughness:.77,side:T.DoubleSide}),atm,state,{sway:.095});
 // Branch cylinders and genuine leaf meshes: moving sideways changes every occlusion.
 const treeRng=C.rng(9246),segments=[],leaves=[];let attempts=0;
 while(state.trees.length<54&&attempts++<1500){let x=(treeRng()-.5)*190,z=17+treeRng()*112;if(C.pathDistance(x,z)<3||C.roofMask(x,z,8)||Math.hypot(x,z-9)<14)continue;const y=C.groundHeight(x,z,world),height=7+treeRng()*7.5,r=.13+height*.018;state.trees.push({x,z,y,height});state.obstacles.push({x,z,r:r+.16});let prev=[x,y,z];const leanX=(treeRng()-.5)*1.5,leanZ=(treeRng()-.5)*1.5;
 for(let j=1;j<=5;j++){const t=j/5,next=[x+leanX*t*t,y+height*t,z+leanZ*t*t];segments.push([prev,next,r*(1-(j-1)*.15),r*(1-j*.15)]);prev=next;}
 for(let b=0;b<11;b++){const t=.45+treeRng()*.47,az=b*2.4+treeRng(),len=(1-t)*height*.7+1.8;let p=[x+leanX*t*t,y+height*t,z+leanZ*t*t];for(let k=1;k<=3;k++){const q=[p[0]+Math.cos(az)*len/3,p[1]+len*(.12+.16*k)/3,p[2]+Math.sin(az)*len/3];segments.push([p,q,r*.22*(1-k*.2),r*.2*(1-k*.25)]);p=q;}
 for(let l=0;l<125;l++){const theta=treeRng()*C.TAU,ct=treeRng()*2-1,rad=Math.cbrt(treeRng()),sx=Math.sqrt(1-ct*ct);leaves.push([p[0]+Math.cos(theta)*sx*rad*1.9,p[1]+ct*rad*1.4,p[2]+Math.sin(theta)*sx*rad*1.9,.8+treeRng()*1.0,treeRng()*C.TAU,treeRng()*C.TAU,treeRng()]);}
 }}
 const trunks=new T.InstancedMesh(new T.CylinderGeometry(1,1,1,7),woodMat,segments.length);const up=new T.Vector3(0,1,0);for(let i=0;i<segments.length;i++){const [a,b,r0,r1]=segments[i],v=new T.Vector3(b[0]-a[0],b[1]-a[1],b[2]-a[2]);obj.position.set((a[0]+b[0])*.5,(a[1]+b[1])*.5,(a[2]+b[2])*.5);obj.quaternion.setFromUnitVectors(up,v.clone().normalize());obj.scale.set((r0+r1)*.5,v.length(),(r0+r1)*.5);obj.updateMatrix();trunks.setMatrixAt(i,obj.matrix);}trunks.castShadow=trunks.receiveShadow=true;trunks.computeBoundingSphere();scene.add(trunks);
 const foliage=new T.InstancedMesh(makeLeaf(T),leafMat,leaves.length),c=new T.Color();for(let i=0;i<leaves.length;i++){const [x,y,z,s,rx,ry,v]=leaves[i];obj.position.set(x,y,z);obj.rotation.set(rx,ry,rx*.7);obj.scale.setScalar(s);obj.updateMatrix();foliage.setMatrixAt(i,obj.matrix);c.setRGB(.68+v*.28,.76+v*.23,.5+v*.24);foliage.setColorAt(i,c);}foliage.castShadow=foliage.receiveShadow=true;foliage.computeBoundingSphere();animateDepth(T,foliage,leafMat.userData.omUniforms,.095);scene.add(foliage);
 // Low plants: crossed curved blades, with a matching displaced shadow shader.
 const grassGeo=new T.BufferGeometry();grassGeo.setAttribute('position',new T.Float32BufferAttribute([-.035,0,0,.035,0,0,-.023,.35,.02,.023,.35,.02,0,.8,.10],3));grassGeo.setIndex([0,1,2,1,3,2,2,3,4]);grassGeo.computeVertexNormals();
 const grassMat=prepareMaterial(T,new T.MeshStandardMaterial({color:0x728647,roughness:.85,side:T.DoubleSide}),atm,state,{sway:.035});const gr=C.rng(146),placements=[];for(let i=0;i<55000;i++){const x=(gr()-.5)*180,z=gr()*125;if(z-C.shore(x)<7||C.pathDistance(x,z)<1.8||C.roofMask(x,z,1))continue;placements.push([x,z,.22+gr()*.65,gr()*C.TAU]);if(placements.length>=26000)break;}
 const grass=new T.InstancedMesh(grassGeo,grassMat,placements.length);for(let i=0;i<placements.length;i++){const [x,z,s,rot]=placements[i];obj.position.set(x,C.groundHeight(x,z,world),z);obj.rotation.set(0,rot,0);obj.scale.set(s,s*1.3,s);obj.updateMatrix();grass.setMatrixAt(i,obj.matrix);c.setRGB(.7+gr()*.3,.7+gr()*.3,.6+gr()*.3);grass.setColorAt(i,c);}grass.receiveShadow=true;grass.castShadow=false;grass.computeBoundingSphere();scene.add(grass);
 // Rain shelter. Walking under this roof changes rain exposure and direct shadows.
 const S=C.SHELTER,deckY=C.groundHeight(S.x,S.z,world)+.22;state.pavilionY=deckY;
 const deckMat=prepareMaterial(T,new T.MeshStandardMaterial({color:0x6d6052,roughness:.88}),atm,state,{sheltered:true});
 const frameMat=prepareMaterial(T,new T.MeshStandardMaterial({color:0x263b3b,metalness:.5,roughness:.38}),atm,state);
 const roofMat=prepareMaterial(T,new T.MeshStandardMaterial({color:0x7d8884,metalness:.4,roughness:.44}),atm,state);
 mesh(T,new T.BoxGeometry(11,.24,8),deckMat,scene,[S.x,deckY-.12,S.z]);for(let i=0;i<28;i++)mesh(T,new T.BoxGeometry(.008,.008,8),frameMat,scene,[S.x-5.3+i*.39,deckY+.002,S.z]);
 for(const dx of[-4.8,4.8])for(const dz of[-3.3,3.3]){cylinderBetween(T,[S.x+dx,deckY,S.z+dz],[S.x+dx,deckY+4.1,S.z+dz],.115,.115,frameMat,scene);state.obstacles.push({x:S.x+dx,z:S.z+dz,r:.25});}
 const roof=mesh(T,new T.BoxGeometry(12,.16,9.2),roofMat,scene,[S.x,deckY+4.2,S.z]);roof.rotation.x=.045;
 for(const z of[-3.4,3.4])mesh(T,new T.BoxGeometry(10.3,.16,.12),frameMat,scene,[S.x,deckY+3.95,S.z+z]);
 mesh(T,new T.BoxGeometry(6.7,.25,.65),deckMat,scene,[S.x,deckY+.58,S.z+2.3]);for(const x of[-2.8,2.8])mesh(T,new T.BoxGeometry(.14,.65,.6),frameMat,scene,[S.x+x,deckY+.28,S.z+2.3]);
 const lamp=new T.PointLight(0xffd29b,0,22,2);lamp.position.set(S.x,deckY+3.5,S.z);scene.add(lamp);state.lamp=lamp;const bulb=mesh(T,new T.CylinderGeometry(.14,.14,.06,16),new T.MeshStandardMaterial({color:0xffdfa7,emissive:0xffdfa7,emissiveIntensity:0}),scene,[S.x,deckY+3.92,S.z]);state.bulb=bulb;
 // Ramped access supports collision-safe entry without a jump requirement.
 const rampGeo=new T.PlaneGeometry(4,7,1,12);rampGeo.rotateX(-Math.PI/2);for(let i=0;i<rampGeo.attributes.position.count;i++){const p=rampGeo.attributes.position,x=p.getX(i)+S.x,z=p.getZ(i)+S.z-7.5,t=C.clamp((z-(S.z-11))/(7));p.setXYZ(i,x,C.mix(C.groundHeight(x,z,world),deckY,C.smooth(0,1,t)),z);}rampGeo.computeVertexNormals();mesh(T,rampGeo,deckMat,scene);state.rampGeo=rampGeo;
 // Shallow wet-film patches with irregular real geometry.
 const puddleMat=prepareMaterial(T,new T.MeshPhysicalMaterial({color:0x33443c,metalness:0,roughness:.12,clearcoat:1,clearcoatRoughness:.06,transparent:true,opacity:0,depthWrite:false}),atm,state);
 const pr=C.rng(512);for(let i=0;i<17;i++){const x=(pr()-.5)*60,z=12+pr()*37;if(C.roofMask(x,z))continue;const g=new T.CircleGeometry(1,30);g.rotateX(-Math.PI/2);const p=g.attributes.position;for(let j=0;j<p.count;j++){const px=p.getX(j)*(1.1+pr()*.55),pz=p.getZ(j)*(.45+pr()*.45);p.setXYZ(j,px,C.groundHeight(x+px,z+pz,world)-C.groundHeight(x,z,world)+.018,pz);}g.computeVertexNormals();const o=mesh(T,g,puddleMat,scene,[x,C.groundHeight(x,z,world),z]);o.castShadow=false;state.puddles.push(o);}state.puddleMaterial=puddleMat;
 // A discreet marker of human scale at the overlook.
 const ox=-56,oz=47,oy=C.groundHeight(ox,oz,world);for(const x of[-5,5])cylinderBetween(T,[ox+x,oy,oz-3],[ox+x,oy+1.05,oz-3],.045,.045,frameMat,scene);cylinderBetween(T,[ox-5,oy+1.05,oz-3],[ox+5,oy+1.05,oz-3],.035,.035,frameMat,scene);
 for(const gm of state.groundMeshes){const p=gm.geometry.attributes.position,a=new Float32Array(p.count);for(let i=0;i<p.count;i++){let cover=0;for(const tr of state.trees){const d=Math.hypot(p.getX(i)-tr.x,p.getZ(i)-tr.z);cover=Math.max(cover,1-C.smooth(2,5.5,d));}a[i]=cover;}gm.geometry.setAttribute('omCanopy',new T.BufferAttribute(a,1));}
 state.counts={trees:state.trees.length,leaves:leaves.length,branchSegments:segments.length,grassBlades:placements.length,rocks:90};
 state.height=(x,z)=>{let y=C.groundHeight(x,z,world);if(C.roofMask(x,z))y=Math.max(y,deckY);if(Math.abs(x-S.x)<2&&z>S.z-11&&z<S.z-4){const t=C.smooth(0,1,(z-S.z+11)/7);y=C.mix(y,deckY,t);}return y;};
 state.canMove=(x,z)=>{if(Math.abs(x)>130||z>140||z<C.shore(x)+1.8)return false;for(const o of state.obstacles)if(Math.hypot(x-o.x,z-o.z)<o.r+.32)return false;return true;};
 state.update=(ledger)=>{for(const m of state.materials){const u=m.userData.omUniforms;if(u){u.uWet.value=C.clamp(ledger.exposed/.9);u.uCanopyWet.value=C.clamp(ledger.canopy/.7);u.uLeafWet.value=C.clamp(ledger.leaf/.35);}}state.puddleMaterial.opacity=C.smooth(.10,.6,ledger.exposed)*.63;};
 return state;
}
OM.createScene=createScene;OM.prepareMaterial=prepareMaterial;
})(globalThis);
