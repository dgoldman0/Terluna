/* Stable habitat-driven plant/stone identities and connected plant architecture.
 * Communities are authored suitability rules for an introduced biosphere.
 * Structural variety and age classes are explicit; growth is not time-simulated.
 */
(function(root){
'use strict';const OM=root.OM=root.OM||{},C=root.OpenMoonCore,L=root.OpenMoonLandscape;
const norm=v=>{const n=Math.hypot(...v);return v.map(x=>x/n);};
const plus=(a,b)=>a.map((x,i)=>x+b[i]);
const mul=(a,t)=>a.map(x=>x*t);
function plan(){
 L.initialize();L.setCanopies([]);const trees=[],rocks=[],gravel=[],tufts=[],understory=[];
 for(let j=0;j<18;j++)for(let i=-15;i<=17;i++){
  const id=`tree:${i}:${j}`,r=C.rng(Math.floor(L.hash(i,j,873)*4294967295)),x=(i+(r()-.5)*.8)*8,z=7+(j+(r()-.5)*.8)*8,f=L.sample(x,z);
  if(C.pathDistance(x,z)<3.7||C.roofMask(x,z,4)||f.elevation<1.2||f.slope>.52||f.soilDepth<.16)continue;
  const probability=.45*f.community[1]+.25*f.community[2]+.025*f.community[0];
  if(r()>probability)continue;
  const family=f.moisture>.45?2:f.exposure>.72&&f.elevation<5?0:1;
  const mature=r()>.24,h=mature?(family===0?5.5:7.7)+r()*4.5:2.1+r()*3.1;
  trees.push({id,x,z,y:C.groundHeight(x,z),height:h,family,ageClass:mature?'mature':'juvenile',crownRadius:h*(family===1?.34:.39),canopyOpacity:family===2?.85:.80,seed:Math.floor(r()*4294967295),habitat:{soilDepth:f.soilDepth,moisture:f.moisture,exposure:f.exposure,suitability:probability}});
 }
 // Two managed edge trees are part of the authored planting scenario. Their
 // improved rooting zones occur in the shared substrate/material field.
 for(const [id,x,z,h,seed]of [['managed:west-edge',-10,-1,7.2,381712],['managed:east-edge',18,0,7.2,292013]]){
  const f=L.sample(x,z);if(f.elevation>.70&&f.soilDepth>.30)trees.push({id,x,z,y:C.groundHeight(x,z),height:h,family:0,ageClass:'mature',crownRadius:h*.39,canopyOpacity:.85,seed,managed:true,habitat:{soilDepth:f.soilDepth,moisture:f.moisture,exposure:f.exposure,suitability:f.community[1]}});
 }
 // A stable identifier priority plus crown-scale spacing, independent of render LOD.
 trees.sort((a,b)=>a.seed-b.seed);const accepted=[];
 for(const t of trees)if(!accepted.some(a=>Math.hypot(a.x-t.x,a.z-t.z)<(a.crownRadius+t.crownRadius)*.55))accepted.push(t);
 L.setCanopies(accepted);
 for(let j=-12;j<34;j++)for(let i=-33;i<34;i++){
  const r=C.rng(Math.floor(L.hash(i,j,218)*4294967295)),x=(i+r())*4,z=(j+r())*4,f=L.sample(x,z);
  if(C.pathDistance(x,z)<1.8||C.roofMask(x,z,1)||f.elevation<-1.0||f.elevation>24)continue;
  if(r()<f.weights[2]*.55+f.weights[1]*.12){const s=.32+r()**1.3*(1.4+f.substrate*2.4);if(C.pathDistance(x,z)<2.0+s*1.2)continue;rocks.push({id:`rock:${i}:${j}`,x,z,s,angle:r()*C.TAU,family:Math.floor(r()*6),stretch:.9+r()*.7});}
 }
 for(let j=-45;j<92;j++)for(let i=-130;i<130;i++){
  const r=C.rng(Math.floor(L.hash(i,j,94)*4294967295)),x=(i+r())*.5,z=(j+r())*.5,f=L.sample(x,z);
  if(f.elevation<-.55||f.elevation>3||C.roofMask(x,z)||r()>f.weights[1]*.8+f.drainage*.18)continue;
  gravel.push({id:`pebble:${i}:${j}`,x,z,s:.018+r()**3*.09,angle:r()*C.TAU,tone:r()});
 }
 for(let j=-2;j<97;j++)for(let i=-83;i<84;i++){
  const r=C.rng(Math.floor(L.hash(i,j,512)*4294967295)),x=(i+r())*1.5,z=(j+r())*1.5;
  if(C.pathDistance(x,z)<1.55||C.roofMask(x,z,1))continue;const f=L.sample(x,z);
  const community=f.community[2]>.2?2:f.canopy>.35?1:0;
  const density=(f.community[0]*.65+f.community[2]*.60+f.community[3]*.18)*(.45+.55*L.noise(x*.085,z*.085,91));
  if(r()<density)tufts.push({id:`tuft:${i}:${j}`,x,z,h:(community===2?.38:.19)+r()*.30,angle:r()*C.TAU,scale:.6+r()*.9,community,tone:r()});
  if(f.canopy>.28&&f.soilDepth>.2&&f.elevation>1.8&&r()<(.035+.12*f.moisture)*f.canopy)understory.push({id:`understory:${i}:${j}`,x,z,h:.22+r()*.4,angle:r()*C.TAU,scale:.5+r()*.75});
 }
 return {trees:accepted,rocks,gravel,tufts,understory,seed:873,version:'community-2'};
}
function skeleton(tree){
 const r=C.rng(tree.seed),segments=[],leaves=[],H=tree.height,base=[tree.x,tree.y,tree.z],family=tree.family,lean=(family===0?.085:.032)*H,az=r()*C.TAU,leanX=Math.cos(az)*lean,leanZ=Math.sin(az)*lean;
 const trunk=t=>[base[0]+leanX*t*t,base[1]+H*t,base[2]+leanZ*t*t];
 const trunkRadius=t=>Math.max(.013,H*.023*(1-t)**.72);
 const branches=tree.ageClass==='juvenile'?5:9;
 const trunkNodes=[...new Set([...Array.from({length:10},(_,i)=>i/9),...Array.from({length:branches},(_,i)=>.32+.56*i/branches)])].sort((a,b)=>a-b);
 for(let i=1;i<trunkNodes.length;i++){const a=trunkNodes[i-1],b=trunkNodes[i];segments.push({a:trunk(a),b:trunk(b),r0:trunkRadius(a),r1:trunkRadius(b)});}
 function limb(a,dir,len,radius,depth){
  const points=[a];let d=dir;
  for(let s=1;s<=3;s++){d=norm([d[0]+(r()-.5)*.16,d[1]+.15,d[2]+(r()-.5)*.16]);const b=plus(points.at(-1),mul(d,len/3));segments.push({a:points.at(-1),b,r0:radius*(1-(s-1)*.22),r1:radius*(1-s*.22)});points.push(b);}
  if(depth>0){for(let i=0;i<3;i++){const theta=az+r()*C.TAU+i*2.399,nd=norm([d[0]*.43+Math.cos(theta)*.66,.15+r()*.48+d[1]*.2,d[2]*.43+Math.sin(theta)*.66]);limb(points[i===0?2:3],nd,len*(.52+r()*.16),radius*.43,depth-1);}}
  else{
   const a=points[1],b=points[3],axis=norm(b.map((x,i)=>x-a[i]));
   for(let j=0;j<11;j++){
    const t=.1+j*.087,part=t<.5?1:2,u=t<.5?t*2:(t-.5)*2,centre=points[part].map((x,i)=>C.mix(x,points[part+1][i],u));
    for(let side=0;side<3;side++){
     const theta=j*2.399+side*2.094+r()*.25,d=norm([Math.cos(theta)*.72+axis[0]*.25,.10+r()*.65,Math.sin(theta)*.72+axis[2]*.25]);
     leaves.push({position:centre,dir:d,roll:r()*.7-.35,size:.7+r()*.65,family,tone:r()});
    }
   }
  }
 }
 for(let i=0;i<branches;i++){
  const t=.32+.56*i/branches,theta=i*2.399+az+(r()-.5)*.5;
  const direction=norm([Math.cos(theta)*.9,(family===1?.24:.08)+r()*.28,Math.sin(theta)*.9]);
  const length=H*(.20+(1-t)*.18)*(family===0?1.12:1);
  limb(trunk(t),direction,length,trunkRadius(t)*.45,2);
 }
 return {segments,leaves};
}
function batch(T,geometry,material,records,transform,scene,name,cell=24,shadow=true){
 const groups=new Map(),o=new T.Object3D(),c=new T.Color();
 for(const r of records){const x=r.x??r.position?.[0]??r.a?.[0]??0,z=r.z??r.position?.[2]??r.a?.[2]??0,key=Math.floor(x/cell)+':'+Math.floor(z/cell);if(!groups.has(key))groups.set(key,[]);groups.get(key).push(r);}
 const output=[];
 for(const [key,items]of groups){const mesh=new T.InstancedMesh(geometry,material,items.length);items.forEach((r,i)=>{o.position.set(0,0,0);o.rotation.set(0,0,0);o.scale.set(1,1,1);c.setRGB(1,1,1);transform(r,o,c);o.updateMatrix();mesh.setMatrixAt(i,o.matrix);mesh.setColorAt(i,c);});mesh.castShadow=shadow;mesh.receiveShadow=true;mesh.name=name+':'+key;mesh.userData.ecology=true;mesh.computeBoundingSphere();scene.add(mesh);output.push(mesh);}
 return output;
}
function leafGeometry(T,family){
 const positions=[],uv=[],indices=[],rows=5,length=family===2?.31:.34,width=family===0?.085:family===1?.12:.07;
 for(let j=0;j<rows;j++){const t=j/(rows-1),w=width*Math.sin(Math.PI*t)**.8;positions.push(-w,.07+t*length,Math.sin(Math.PI*t)*.015,w,.07+t*length,Math.sin(Math.PI*t)*.015);uv.push(0,t,1,t);if(j<rows-1){const a=j*2;indices.push(a,a+1,a+2,a+1,a+3,a+2);}}
 const stem=positions.length/3;positions.push(-.0018,0,0,.0018,0,0,-.0018,.071,0,.0018,.071,0);uv.push(.5,0,.5,0,.5,.05,.5,.05);indices.push(stem,stem+1,stem+2,stem+1,stem+3,stem+2);
 const geo=new T.BufferGeometry();geo.setAttribute('position',new T.Float32BufferAttribute(positions,3));geo.setAttribute('uv',new T.Float32BufferAttribute(uv,2));geo.setIndex(indices);geo.computeVertexNormals();return geo;
}
function tuftGeometry(T,family){
 const r=C.rng(456+family),p=[],uv=[],idx=[],blades=family===2?11:15;
 for(let b=0;b<blades;b++){
  const a=b*2.399+r()*.7,bend=.20+r()*.50,len=.50+r()*.50,width=family===2?.015:.009,x=(r()-.5)*.2,z=(r()-.5)*.2,start=p.length/3;
  for(let j=0;j<=4;j++){const t=j/4,w=width*(1-t)+.0001,reach=bend*t*t;
   for(const sign of [-1,1]){p.push(x+Math.cos(a)*reach+Math.sin(a)*w*sign,t*len,z+Math.sin(a)*reach-Math.cos(a)*w*sign);uv.push((sign+1)/2,t);}
   if(j<4){const i=start+j*2;idx.push(i,i+1,i+2,i+1,i+3,i+2);}
  }
 }
 const geo=new T.BufferGeometry();geo.setAttribute('position',new T.Float32BufferAttribute(p,3));geo.setAttribute('uv',new T.Float32BufferAttribute(uv,2));geo.setIndex(idx);geo.computeVertexNormals();return geo;
}
function fernGeometry(T){
 const p=[],uv=[],idx=[];
 for(let frond=0;frond<7;frond++){
  const angle=frond*2.399,dx=Math.cos(angle),dz=Math.sin(angle);
  for(let j=1;j<10;j++){
   const t=j/10,cx=dx*t*.65,cz=dz*t*.65,cy=Math.sin(t*Math.PI*.86)*.75,len=Math.sin(t*Math.PI)*.19;
   for(const side of [-1,1]){const a=p.length/3;p.push(cx,cy,cz,cx-dz*len*side-dx*.035,cy-.035,cz+dx*len*side-dz*.035,cx-dz*len*side+dx*.035,cy-.045,cz+dx*len*side+dz*.035);uv.push(.5,0,0,.8,1,1);idx.push(a,a+1,a+2);}
  }
 }
 const geo=new T.BufferGeometry();geo.setAttribute('position',new T.Float32BufferAttribute(p,3));geo.setAttribute('uv',new T.Float32BufferAttribute(uv,2));geo.setIndex(idx);geo.computeVertexNormals();return geo;
}
function populate(T,scene,atm,state,layout){
 const wood=OM.material(T,atm,state,'wood'),rockMat=OM.material(T,atm,state,'rock'),up=new T.Vector3(0,1,0),v=new T.Vector3(),leafMat=OM.material(T,atm,state,'leaf',{color:0x647e37,side:T.DoubleSide,roughness:.78});
 const allSegments=[],allLeaves=[];
 for(const tree of layout.trees){const model=skeleton(tree);allSegments.push(...model.segments);allLeaves.push(...model.leaves);state.trees.push(tree);state.obstacles.push({x:tree.x,z:tree.z,r:tree.height*.026+.04});}
 const ratios=[.40,.60,.80,.95];
 for(let k=0;k<ratios.length;k++){
  const segments=allSegments.filter(s=>{if(s.r0<.003)return false;const ratio=s.r1/s.r0;let b=0;for(let j=1;j<ratios.length;j++)if(Math.abs(ratios[j]-ratio)<Math.abs(ratios[b]-ratio))b=j;return b===k;});
  const geo=new T.CylinderGeometry(ratios[k],1,1,8,2);
  batch(T,geo,wood,segments,(r,o)=>{v.set(r.b[0]-r.a[0],r.b[1]-r.a[1],r.b[2]-r.a[2]);const len=v.length();o.position.set(...r.a.map((x,i)=>(x+r.b[i])*.5));o.quaternion.setFromUnitVectors(up,v.normalize());o.scale.set(r.r0,len,r.r0);},scene,'branches');
 }
 for(let family=0;family<3;family++){
  const leaves=allLeaves.filter(l=>l.family===family),meshes=batch(T,leafGeometry(T,family),leafMat,leaves,(r,o,c)=>{o.position.set(...r.position);v.set(...r.dir);o.quaternion.setFromUnitVectors(up,v.normalize());o.rotateY(r.roll);o.scale.setScalar(r.size);c.setRGB(.76+r.tone*.28,.83+r.tone*.22,.58+r.tone*.31);},scene,'leaves-'+family);
  for(const m of meshes)OM.animateDepth(T,m,leafMat.userData.omUniforms,.015);
 }
 for(let family=0;family<6;family++){
  const geo=new T.SphereGeometry(1,28,18),p=geo.attributes.position,r=C.rng(912+family),planes=[];
  for(let j=0;j<12;j++){const n=norm([r()*2-1,r()*2-1,r()*2-1]);planes.push([...n,.68+r()*.3]);}
  for(let i=0;i<p.count;i++){const x=p.getX(i),y=p.getY(i),z=p.getZ(i);let radius=1;for(const [nx,ny,nz,d]of planes){const dot=x*nx+y*ny+z*nz;if(dot>0)radius=Math.min(radius,d/dot);}radius+=(L.noise(x*5+family,z*6+y*4,8)-.5)*.035;p.setXYZ(i,x*radius,y*radius*.72,z*radius);}
  geo.computeVertexNormals();const rocks=layout.rocks.filter(r=>r.family===family);
  batch(T,geo,rockMat,rocks,(r,o)=>{o.position.set(r.x,C.groundHeight(r.x,r.z,state.world)+r.s*.16,r.z);o.rotation.set(.1*Math.sin(r.angle),r.angle,.13*Math.cos(r.angle));o.scale.set(r.s*r.stretch,r.s,r.s*.9);if(r.s>.65)state.obstacles.push({x:r.x,z:r.z,r:r.s*.58});},scene,'substrate-rocks');
 }
 const pebbleMat=OM.material(T,atm,state,'rock');
 batch(T,new T.IcosahedronGeometry(1,1),pebbleMat,layout.gravel,(r,o,c)=>{o.position.set(r.x,C.groundHeight(r.x,r.z,state.world)+r.s*.22,r.z);o.scale.set(r.s*1.4,r.s*.65,r.s);o.rotation.set(.3,r.angle,.1);c.setRGB(.75+r.tone*.50,.72+r.tone*.48,.67+r.tone*.42);},scene,'deposited-gravel',18,true);
 for(let family=0;family<3;family++){
  const mat=OM.material(T,atm,state,'grass',{color:family===2?0x688442:family===1?0x60753c:0x85824a,side:T.DoubleSide});
  const groups=batch(T,tuftGeometry(T,family),mat,layout.tufts.filter(t=>t.community===family),(r,o,c)=>{o.position.set(r.x,C.groundHeight(r.x,r.z,state.world),r.z);o.scale.set(r.scale,r.h,r.scale);o.rotation.y=r.angle;c.setRGB(.80+r.tone*.3,.84+r.tone*.24,.65+r.tone*.31);},scene,'community-grass-'+family,18,false);
  for(const m of groups)OM.animateDepth(T,m,mat.userData.omUniforms,.025);
 }
 const fernMat=OM.material(T,atm,state,'leaf',{color:0x5a7639,side:T.DoubleSide});
 const ferns=batch(T,fernGeometry(T),fernMat,layout.understory,(r,o)=>{o.position.set(r.x,C.groundHeight(r.x,r.z,state.world),r.z);o.scale.set(r.scale,r.h*2,r.scale);o.rotation.y=r.angle;},scene,'woodland-understory',18,true);
 for(const m of ferns)OM.animateDepth(T,m,fernMat.userData.omUniforms,.015);
 state.ecology=layout;return {trees:layout.trees.length,leaves:allLeaves.length,branchSegments:allSegments.length,grassClumps:layout.tufts.length,understory:layout.understory.length,rocks:layout.rocks.length,pebbles:layout.gravel.length};
}
OM.Ecology={plan,skeleton,populate,leafGeometry,tuftGeometry};if(typeof module!=='undefined'&&module.exports)module.exports=OM.Ecology;
})(globalThis);
