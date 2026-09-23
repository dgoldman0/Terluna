'use strict';
const test=require('node:test'),assert=require('node:assert/strict');
const W=require('../weather-column.js');
const near=(a,b,e=1e-7)=>assert.ok(Math.abs(a-b)<e,`${a} vs ${b} exceeds ${e}`);
const models={};function get(k='fair',w='moon'){return models[k+w]??=W.create(k,w);}
test('Water and ice saturation share the thermodynamic triple point',()=>{near(W.es(W.K.Tt),W.K.et);near(W.es(W.K.Tt,'ice'),W.K.et);});
test('Saturation matches the independent MetPy/Ambaum example within 0.15% (different e0 and rounded Cp)',()=>near(W.es(298.15)/100,31.623456,.047));
test('Ice saturation lies below liquid saturation at subfreezing temperatures',()=>{for(const T of[240,250,260])assert.ok(W.es(T,'ice')<W.es(T));});
test('Parcel LCL agrees with the independently documented Romps/MetPy case',()=>{const q=W.lcl(306.15,94300,W.es(301.15)/W.es(306.15));near(q.p/100,877.033549,.20);near(q.T,299.9091908,.03);});
test('Saturated surface air has zero lifting-condensation displacement',()=>{const a=W.lcl(288,121590,1);near(a.x,0,1e-12);near(a.T,288,1e-9);});
test('Saturation pressure rises monotonically through the working interval',()=>{let previous=0;for(let t=180;t<=330;t++){const p=W.es(t);assert.ok(p>previous);previous=p;}});
test('Invalid humidity, temperature, pressure, world and step fail early',()=>{for(const h of[0,-.1,1.1,NaN])assert.throws(()=>W.lcl(294,121590,h));assert.throws(()=>W.es(80));assert.throws(()=>W.rs(100,300));assert.throws(()=>W.create('unknown'));assert.throws(()=>W.create('fair','mars'));assert.throws(()=>W.create('fair','moon',{dx:0}));});
test('Liquid pseudoadiabat matches the independent MetPy pressure-profile example',()=>{
 const levels=[925,850,700,500,300,200],expected=[5,.99635104,-8.88958079,-28.38862857,-60.12003999,-83.34321585];let T=278.15,x=0;
 for(let i=0;i<levels.length;i++){const target=Math.log(925/levels[i]);while(x<target-1e-12){const h=Math.min(.0005,target-x);T=W.rkMoist(T,x,h,92500);x+=h;}near(T-273.15,expected[i],.14);}
});
test('Spherical hydrostatics reduce to the thin-shell limit near the ground',()=>{const z=W.heightFromPotential(1000,1737400,1.62);near(z,1000/1.62,.23);});
test('Spherical gravity inversion matches the exact potential',()=>{for(const z of[20,10000,60000]){const R=1737400,g=1.62,phi=g*R*z/(R+z);near(W.heightFromPotential(phi,R,g),z,1e-8);}});
test('All soundings have ordered pressure/height and positive thermodynamic states',()=>{for(const w of['moon','earth'])for(const k of Object.keys(W.PRESETS)){const c=get(k,w);let prevZ=-1,prevP=Infinity;for(const r of c.rows){assert.ok(r.z>prevZ&&r.p<prevP);prevZ=r.z;prevP=r.p;assert.ok(r.T>150&&r.rho>0&&r.rh>=0&&r.rh<=1.000001);assert.ok(Object.values(r).every(Number.isFinite));}}});
test('Pressure drop equals integrated density times spherical gravity',()=>{for(const w of['moon','earth'])for(const k of Object.keys(W.PRESETS))assert.ok(get(k,w).summary.hydrostaticWeightRelativeError<2e-7);});
test('Relative humidity uses partial pressure, not the mixing-ratio approximation',()=>{for(const c of Object.values(models))for(const r of c.rows.filter((r,i)=>i%67===0))near(r.rh,(r.p*r.r/(W.EPS+r.r))/W.es(r.T),1e-10);});
test('A low-gravity column has the expected dry lapse rate and large scale height',()=>{const m=get().summary,e=get('fair','earth').summary;near(m.dryLapse_K_per_km,1.6135458167,1e-8);assert.ok(m.surfaceScaleHeight_m>50000);assert.ok(m.lcl_m/e.lcl_m>6);});
test('Heating-free reference parcel reproduces the roughly 8.6 km lunar LCL diagnostic',()=>{const spec=JSON.parse(JSON.stringify(W.PRESETS.fair));spec.parcelHeating=0;const s=W.create('fair','moon',{preset:spec}).summary;assert.ok(s.lcl_m>8100&&s.lcl_m<8900);});
test('Fair-weather cloud is bounded by condensation and the prescribed inversion',()=>{const c=get('fair');assert.ok(c.summary.parcelReachesCloud);assert.ok(c.summary.cloudBase_m>9000&&c.summary.cloudTop_m<14000);for(const r of c.rows)if(r.z<c.summary.cloudBase_m||r.z>c.summary.cloudTop_m)assert.equal(r.extinction,0);});
test('Deep convection has a high top and both liquid and ice optical regions',()=>{const s=get('convection').summary;assert.ok(s.cloudTop_m>50000&&s.cloudTop_m<70000);assert.ok(s.iceWaterPath_kg_m2>0&&s.liquidWaterPath_kg_m2>0);});
test('Coastal fog remains a low saturated layer beneath an inversion',()=>{const c=get('fog');assert.equal(c.summary.cloudBase_m,0);assert.ok(c.summary.cloudTop_m>50&&c.summary.cloudTop_m<300);assert.ok(c.rows[0].extinction>0);assert.equal(c.summary.iceWaterPath_kg_m2,0);});
test('A stable column can have an LCL without reachable convective cloud',()=>{const p=JSON.parse(JSON.stringify(W.PRESETS.fair));p.parcelHeating=0;p.launchSpeed=0;p.temperature=p.temperature.map(([x])=>[x,294]);const c=W.create('fair','moon',{preset:p});assert.ok(c.summary.lcl_m>0);assert.equal(c.summary.cloudBase_m,null);assert.equal(c.summary.inCloudOpticalDepth,0);});
test('Halving vertical grid spacing converges cloud boundaries and optical depth',()=>{for(const k of['fair','convection','fog']){const a=get(k),b=W.create(k,'moon',{dx:.0005});near(a.summary.lcl_m,b.summary.lcl_m,.1);near(a.summary.cloudTop_m,b.summary.cloudTop_m,55);assert.ok(Math.abs(a.summary.inCloudOpticalDepth/b.summary.inCloudOpticalDepth-1)<.025);}});
test('Texture channels sample diagnosed extinction, ice fraction and wind shear',()=>{const c=get('convection'),t=W.textureData(c);for(let i=0;i<t.n;i++){const r=c.sample(t.base+(t.top-t.base)*i/(t.n-1));near(t.data[i*4],r.extinction,1e-7);near(t.data[i*4+2],r.windX,2e-6);near(t.data[i*4+3],r.windZ,2e-6);}});
test('Diagnostic recreation is deterministic and does not mutate prescribed inputs',()=>{const before=JSON.stringify(W.PRESETS);assert.deepEqual(W.create('fair').rows,W.create('fair').rows);assert.equal(JSON.stringify(W.PRESETS),before);});
