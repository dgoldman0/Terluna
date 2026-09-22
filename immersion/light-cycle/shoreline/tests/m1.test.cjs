'use strict';
const {test, after}=require('node:test');
const assert=require('node:assert/strict');
const fs=require('node:fs');
const path=require('node:path');
const crypto=require('node:crypto');
const C=require('../src/core.js');
const T=require('../vendor/three.cjs');
globalThis.OM={}; require('../src/scene.js');
const root=path.join(__dirname,'..');
const close=(a,b,tol=1e-7)=>assert.ok(Math.abs(a-b)<=tol, `${a} != ${b}, tolerance ${tol}`);
let terrain;
function geometry(){return terrain ||= OM.terrainGeometry(T,'moon');}
after(()=>terrain?.dispose());
test('M1: supplied Three.js has exact pinned Git blob identity',()=>{
 const bytes=fs.readFileSync(path.join(root,'vendor/three.cjs'));
 const hash=crypto.createHash('sha1').update(`blob ${bytes.length}\0`).update(bytes).digest('hex');
 assert.equal(hash,'ca4833532c363b72477b2e8a6f47cc0e2fc7b09a');assert.equal(T.REVISION,'180');
});
test('M1: terrain field stays continuous across the coast guide',()=>{
 const e=1e-5;for(let x=-500;x<=500;x+=7){const z=C.shore(x);close(C.surfaceHeight(x,z-e),C.surfaceHeight(x,z+e),1e-4);}
});
test('M1: land and seabed depths use the same sea datum',()=>{
 for(let x=-2000;x<=2000;x+=71)for(let z=-4000;z<=400;z+=89){close(C.waterDepth(x,z),Math.max(0,-C.surfaceHeight(x,z)));close(C.groundHeight(x,z)+C.curvatureSag(x,z),C.surfaceHeight(x,z));}
});
test('M1: submerged points have positive depth and regional relief emerges',()=>{
 assert.ok(C.waterDepth(0,-180)>0);assert.ok(C.surfaceHeight(1280,-2800)>0);assert.equal(C.waterDepth(1280,-2800),0);
});
test('M1: spherical sag is zero at origin, symmetric and radius-dependent',()=>{
 assert.equal(C.curvatureSag(0,0),0);close(C.curvatureSag(1500,600),C.curvatureSag(-1500,-600));assert.ok(C.curvatureSag(1500,600)>C.curvatureSag(1500,600,'earth'));
 const r=15000,R=C.worldRadius();close(C.curvatureSag(r,0),R-Math.sqrt(R*R-r*r),1e-8);
});
test('M1: 196609 terrain vertices are finite and retain unit upward normals',()=>{
 const g=geometry(),p=g.attributes.position,n=g.attributes.normal;
 assert.equal(p.count,196609);assert.equal(n.count,p.count);
 for(let i=0;i<p.count;i++){assert.ok([p.getX(i),p.getY(i),p.getZ(i)].every(Number.isFinite));close(Math.hypot(n.getX(i),n.getY(i),n.getZ(i)),1,2e-6);assert.ok(n.getY(i)>0);}
});
test('M1: terrain vertices sample the authoritative field at every radial scale',()=>{
 const p=geometry().attributes.position;
 for(let i=0;i<p.count;i+=127){close(p.getY(i),C.groundHeight(p.getX(i),p.getZ(i)),.025);}
});
test('M1: all terrain triangle indices are valid and have consistent upward winding',()=>{
 const g=geometry(),p=g.attributes.position,idx=g.index.array;
 for(let i=0;i<idx.length;i+=3){const [a,b,c]=idx.slice(i,i+3);assert.ok(a<p.count&&b<p.count&&c<p.count);
 const area=(p.getZ(b)-p.getZ(a))*(p.getX(c)-p.getX(a))-(p.getX(b)-p.getX(a))*(p.getZ(c)-p.getZ(a));assert.ok(area>0,`degenerate or reversed triangle ${i/3}`);}
});
test('M1: the mesh has a single 384-edge outer boundary and no internal cracks',()=>{
 const g=geometry(),idx=g.index.array,edges=new Map(),N=g.attributes.position.count;
 for(let i=0;i<idx.length;i+=3)for(let j=0;j<3;j++){let a=idx[i+j],b=idx[i+(j+1)%3];if(a>b)[a,b]=[b,a];const k=a*N+b;edges.set(k,(edges.get(k)||0)+1);}
 let boundary=0;for(const [k,n]of edges){assert.ok(n===1||n===2);if(n===1){boundary++;const a=Math.floor(k/N),b=k%N;assert.ok(a>=N-384&&b>=N-384);}}
 assert.equal(boundary,384);
});
test('M1: clear-noon forcing has zero cloud water, coverage and rainfall',()=>{
 for(const t of[0,1800,8100,14400]){const w=C.weatherAt(t,'clear');assert.equal(w.coverage,0);assert.equal(w.lwc,0);assert.equal(w.tau,0);assert.equal(w.rain,0);assert.equal(w.wind,1.5);}
});
test('M1: all twelve wave bands remain finite over the scene and month',()=>{
 assert.equal(C.WATER_BANDS.length,12);
 for(const t of[0,1000,1e5,C.PERIOD])for(const x of[-16384,-8,0,330,16384]){const w=C.renderWaveAt(x,-x/2,t);assert.ok(Object.values(w).every(Number.isFinite));assert.ok(Math.abs(w.y)<.3);}
});
test('M1: twelve-band wave normals match height finite differences',()=>{
 const e=1e-5;for(const t of[0,41]){const x=3,z=-6,w=C.renderWaveAt(x,z,t);close(w.nx,-(C.renderWaveAt(x+e,z,t).y-C.renderWaveAt(x-e,z,t).y)/(2*e),1e-7);close(w.nz,-(C.renderWaveAt(x,z+e,t).y-C.renderWaveAt(x,z-e,t).y)/(2*e),1e-7);}
});
test('M1: wave timing is consistent with square-root gravity scaling',()=>{
 const a=C.renderWaveAt(4,7,50,1.5,1.62),b=C.renderWaveAt(4,7,50*Math.sqrt(1.62/9.80665),1.5,9.80665);
 close(a.y,b.y);close(a.nx,b.nx);close(a.nz,b.nz);
});
test('M1: overflowing GLSL hyperbolic tangent is replaced by its bounded equivalent',()=>{
 const src=fs.readFileSync(path.join(root,'src/water.js'),'utf8');assert.ok(src.includes('exp(-2.*k*depth)'));assert.ok(!src.includes('tanh('));
 for(const k of C.WATER_BANDS){const e=Math.exp(-2*k*12);close(Math.sqrt(1.62*k*(1-e)/(1+e)),C.waveOmega(k,12,1.62));}
});
test('M1: HTML embeds verified vendor and all eight runtime modules',()=>{
 const html=fs.readFileSync(path.join(root,'Open_Moon_Shoreline.html'),'utf8');
 const encoded=html.match(/<script id="three-vendor"[^>]*>([^<]+)<\/script>/)[1];assert.deepEqual(Buffer.from(encoded,'base64'),fs.readFileSync(path.join(root,'vendor/three.cjs')));
 for(const name of['core','atmosphere','materials','scene','water','audio','app','loader'])assert.ok(html.includes(fs.readFileSync(path.join(root,`src/${name}.js`),'utf8').trim()));
 assert.ok(!/__SKY_DATA__|__THREE_VENDOR__|__MATERIALS__/.test(html));
});

test('M1: clear environment cache ignores observer drift but follows optical state',()=>{
 require('../src/atmosphere.js');const a=new OM.Atmosphere(T,{worlds:{}});const first=a.clearEnvironmentKey();
 assert.ok(first);a.uniforms.uObserver.value.set(20,30);a.uniforms.uEyeHeight.value=12;a.uniforms.uTime.value=88;assert.equal(a.clearEnvironmentKey(),first);
 a.uniforms.uSun.value.set(.1,.99,0);assert.notEqual(a.clearEnvironmentKey(),first);
 a.uniforms.uCover.value=.5;a.uniforms.uTau.value=12;assert.equal(a.clearEnvironmentKey(),null);
});
