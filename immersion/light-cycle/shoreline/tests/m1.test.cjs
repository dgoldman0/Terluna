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
// Moving-grid topology and geometry invariants are in landscape.test.cjs.
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
test('M1: HTML embeds verified vendor and all twelve runtime modules',()=>{
 const html=fs.readFileSync(path.join(root,'Open_Moon_Shoreline.html'),'utf8');
 const encoded=html.match(/<script id="three-vendor"[^>]*>([^<]+)<\/script>/)[1];assert.deepEqual(Buffer.from(encoded,'base64'),fs.readFileSync(path.join(root,'vendor/three.cjs')));
 for(const name of['landscape','core','atmosphere','materials','terrain','surface-water','ecology','scene','water','audio','app','loader'])assert.ok(html.includes(fs.readFileSync(path.join(root,`src/${name}.js`),'utf8').trim()));
 assert.ok(!/__SKY_DATA__|__THREE_VENDOR__|__MATERIALS__/.test(html));
});

test('M1: clear environment cache ignores observer drift but follows optical state',()=>{
 require('../src/atmosphere.js');const a=new OM.Atmosphere(T,{worlds:{}});const first=a.clearEnvironmentKey();
 assert.ok(first);a.uniforms.uObserver.value.set(20,30);a.uniforms.uEyeHeight.value=12;a.uniforms.uTime.value=88;assert.equal(a.clearEnvironmentKey(),first);
 a.uniforms.uSun.value.set(.1,.99,0);assert.notEqual(a.clearEnvironmentKey(),first);
 a.uniforms.uCover.value=.5;a.uniforms.uTau.value=12;assert.equal(a.clearEnvironmentKey(),null);
});
