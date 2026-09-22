/* SI, deterministic column diagnostics; no renderer or simulation-clock dependency.
 * Pressure-coordinate soundings are selected experiments, not a lunar forecast.
 * Hydrostatic heights include spherical g. Parcel temperatures follow dry / liquid
 * pseudoadiabats. Diagnosed retained condensate and mixed-phase partition are
 * explicit optical closures, not a prognostic precipitation budget.
 * Equations and source attribution: CLOUD_METHODS.md.
 */
(function(root){
'use strict';
const K=Object.freeze({Rd:287.05,Rv:461.5,cpd:1004,cpv:1850,cpl:4218,cpi:2106,Tt:273.16,et:611.657,Lv:2500800,Ls:2834500});
const EPS=K.Rd/K.Rv,clamp=(v,a=0,b=1)=>Math.max(a,Math.min(b,v)),mix=(a,b,t)=>a+(b-a)*t;
function finite(v,label){if(!Number.isFinite(v))throw new TypeError(label+' must be finite');return v;}
function es(T,phase='liquid'){
 if(!Number.isFinite(T)||T<150||T>350)throw new RangeError('Saturation temperature must be 150–350 K');
 if(!['liquid','ice'].includes(phase))throw new RangeError('Unknown condensate phase');
 const cp=phase==='ice'?K.cpi:K.cpl,L0=phase==='ice'?K.Ls:K.Lv,dcp=K.cpv-cp,L=L0+dcp*(T-K.Tt);
 return K.et*Math.pow(K.Tt/T,(cp-K.cpv)/K.Rv)*Math.exp(L0/(K.Rv*K.Tt)-L/(K.Rv*T));
}
function rs(p,T,phase='liquid'){finite(p,'Pressure');const e=es(T,phase);if(p<=e)throw new RangeError('Vapour pressure must be smaller than total pressure');return EPS*e/(p-e);}
function virtualT(T,r){return T*(1+r/EPS)/(1+r);}
function moistSlope(p,T){const r=rs(p,T),L=K.Lv;return -(K.Rd*T+L*r)/(K.cpd+L*L*r*EPS/(K.Rd*T*T));}
function rkMoist(T,x,dx,ps){const f=(x,t)=>moistSlope(ps*Math.exp(-x),t),a=f(x,T),b=f(x+dx/2,T+a*dx/2),c=f(x+dx/2,T+b*dx/2),d=f(x+dx,T+c*dx);return T+dx*(a+2*b+2*c+d)/6;}
function lcl(T,p,rh){
 if(!Number.isFinite(rh)||rh<=0||rh>1)throw new RangeError('Parcel relative humidity must be in (0,1]');
 const e=rh*es(T);if(p<=e)throw new RangeError('Invalid parcel pressure');
 const r=EPS*e/(p-e),kappa=(K.Rd+r*K.Rv)/(K.cpd+r*K.cpv);
 let lo=0,hi=Math.min(2,Math.log(T/150)/kappa-.000001);if(hi<=0||rs(p*Math.exp(-hi),T*Math.exp(-kappa*hi))>r)throw new RangeError('LCL falls outside the 150 K parcel limit');
 for(let i=0;i<60;i++){const x=(lo+hi)/2;if(rs(p*Math.exp(-x),T*Math.exp(-kappa*x))>r)lo=x;else hi=x;}
 const x=(lo+hi)/2;return {x,p:p*Math.exp(-x),T:T*Math.exp(-kappa*x),mixingRatio:r,kappa};
}
function interp(knots,x,col=1){if(x<=knots[0][0])return knots[0][col];let i=1;while(i<knots.length-1&&knots[i][0]<x)i++;return mix(knots[i-1][col],knots[i][col],clamp((x-knots[i-1][0])/(knots[i][0]-knots[i-1][0])));}
const PRESETS={
 fair:{label:'Elevated fair-weather cloud',surfaceT:294,rh:.48,parcelHeating:.8,launchSpeed:1.5,retention:.12,
  temperature:[[0,294],[.23,275.2],[.29,284],[.7,269],[1.5,230],[2.2,216],[4,210]],
  humidity:[[0,.48],[.2,.48],[.3,.25],[.8,.35],[2,.15],[4,.02]],wind:[[0,1.5,.3],[.25,6,1],[.7,13,6],[2,20,12],[4,20,12]],
  morphology:{type:1,aspect:1.8,coverage:.42,iceRadius_um:35,liquidRadius_um:12}},
 convection:{label:'Deep convective column',surfaceT:299,rh:.78,parcelHeating:1,launchSpeed:2,retention:.055,
  temperature:[[0,299],[.18,285],[.55,265],[1.10,236],[1.45,243],[2.2,221],[4,210]],
  humidity:[[0,.78],[.18,.72],[.7,.62],[1.5,.30],[4,.02]],wind:[[0,3,.6],[.3,7,2],[.8,17,6],[1.5,30,18],[4,35,20]],
  morphology:{type:2,aspect:.4,coverage:.27,iceRadius_um:45,liquidRadius_um:14}},
 fog:{label:'Coastal inversion fog',surfaceT:288,rh:1,parcelHeating:0,launchSpeed:0,retention:0,
  temperature:[[0,288],[.006,289.4],[.02,292],[.12,290],[.7,265],[1.5,230],[2.2,216],[4,210]],
  humidity:[[0,1],[.006,1],[.025,.8],[.2,.45],[1,.35],[4,.02]],wind:[[0,1.2,.2],[.03,2,.6],[.2,5,1],[1,15,7],[4,20,12]],
  fogPreCooling_K:.8,morphology:{type:0,aspect:5,coverage:1,iceRadius_um:35,liquidRadius_um:8}}
};
// Deep freeze the experiment inputs. Callers customize through create(...,{preset}).
function freeze(o){Object.values(o).forEach(v=>{if(v&&typeof v==='object')freeze(v);});return Object.freeze(o);}freeze(PRESETS);
function planet(world){if(!['moon','moon_no_ozone','earth'].includes(world))throw new RangeError('Unknown world');return world==='earth'?{radius:6371000,g0:9.80665,ps:101325}:{radius:1737400,g0:1.62,ps:121590};}
function heightFromPotential(phi,R,g){if(phi<0||phi>=g*R)throw new RangeError('Hydrostatic potential outside finite spherical column');return R*phi/(g*R-phi);}
function create(key='fair',world='moon',options={}){
 const spec=options.preset||PRESETS[key];if(!spec)throw new RangeError('Unknown column scenario');
 const P={...planet(world),...options.planet},dx=options.dx??.001;
 if(!Number.isFinite(dx)||dx<=0||dx>.02)throw new RangeError('Log-pressure step must be (0,.02]');
 for(const k of ['radius','g0','ps'])if(!Number.isFinite(P[k])||P[k]<=0)throw new RangeError('Invalid planet '+k);
 const T0=spec.surfaceT+spec.parcelHeating;
 // Parcel heating keeps vapour mixing ratio; its relative humidity falls.
 const parcelRH=spec.rh*es(spec.surfaceT)/es(T0),lift=lcl(T0,P.ps,parcelRH),rows=[];
 let phi=0,prevTv=0,parcel=T0,prevX=0,kinetic=spec.launchSpeed**2/2,reachable=true,lfc=null,el=null,cloudSeen=false,cin=0,cape=0;
 const fogTotal=spec.fogPreCooling_K?rs(P.ps,spec.surfaceT+spec.fogPreCooling_K):0;
 const grid=new Set([0,2.2]);if(lift.x>1e-10)grid.add(lift.x);for(let x=dx;x<2.2;x+=dx)grid.add(x);for(let x=.00005;x<.03;x+=.00005)grid.add(x);
 const xs=[...grid].filter(x=>x<=2.2).sort((a,b)=>a-b).filter((x,i,a)=>i===0||x-a[i-1]>1e-10);
 for(let i=0;i<xs.length;i++){
  const x=xs[i],p=P.ps*Math.exp(-x),T=interp(spec.temperature,x),rh=clamp(interp(spec.humidity,x));
  const sat=rs(p,T),vaporPressure=rh*es(T),r=key==='fog'&&x<.03?Math.min(sat,fogTotal*Math.exp(-x/.18)):EPS*vaporPressure/(p-vaporPressure);
  const Tv=virtualT(T,r),h=i?x-prevX:0;
  if(i)phi+=K.Rd*(prevTv+Tv)*h/2;
  const z=heightFromPotential(phi,P.radius,P.g0),g=P.g0*(P.radius/(P.radius+z))**2;
  if(x<=lift.x)parcel=T0*Math.exp(-lift.kappa*x);
  else if(prevX<lift.x)parcel=rkMoist(lift.T,lift.x,x-lift.x,P.ps);
  else parcel=rkMoist(parcel,prevX,h,P.ps);
  const pr=x>=lift.x?rs(p,parcel):lift.mixingRatio,B=g*(virtualT(parcel,pr)-Tv)/Tv;
  if(i){const last=rows.at(-1),work=(last.buoyancy+B)/2*(z-last.z);if(reachable){kinetic+=work;if(kinetic<0)reachable=false;}if(x>=lift.x)cape+=Math.max(0,work);if(lfc===null)cin+=Math.min(0,work);}
  if(x>=lift.x&&B>0&&lfc===null)lfc=z;
  if(cloudSeen&&B<=0&&el===null)el=z;
  const allowed=x>lift.x&&reachable&&B>0&&el===null&&key!=='fog';
  if(allowed)cloudSeen=true;
  // A selected retention fraction of the undiluted parcel's total condensation
  // diagnoses optical condensate; it is independent of ground rainfall forcing.
  let qc=allowed?spec.retention*Math.max(0,lift.mixingRatio-pr)*Math.exp(-(x-lift.x)/1.1):0;
  if(key==='fog')qc=x<.03?Math.max(0,fogTotal*Math.exp(-x/.18)-sat):0;
  const ice=key==='fog'?0:clamp((273.15-parcel)/40),rho=p/(K.Rd*Tv),condensate=rho*qc/(1+r);
  const liquid=condensate*(1-ice),frozen=condensate*ice;
  const ext=1.5*(liquid/(1000*spec.morphology.liquidRadius_um*1e-6)+frozen/(917*spec.morphology.iceRadius_um*1e-6));
  rows.push({x,z,p,T,rh:(p*r/(EPS+r))/es(T),r,Tv,rho,g,parcelT:parcel,parcelR:pr,buoyancy:B,liquid,ice:frozen,extinction:ext,windX:interp(spec.wind,x,1),windZ:interp(spec.wind,x,2)});
  prevX=x;prevTv=Tv;
 }
 const sample=z=>sampleRows(rows,z);
 const wet=rows.filter(r=>r.extinction>1e-9),base=wet.length?(key==='fog'?0:sampleX(rows,lift.x).z):null,top=wet.length?rows[Math.min(rows.length-1,rows.indexOf(wet.at(-1))+1)].z:null;
 let liquidPath=0,icePath=0,tau=0,weight=0;
 for(let i=1;i<rows.length;i++){const a=rows[i-1],b=rows[i],dz=b.z-a.z;liquidPath+=(a.liquid+b.liquid)*dz/2;icePath+=(a.ice+b.ice)*dz/2;tau+=(a.extinction+b.extinction)*dz/2;weight+=(a.rho*a.g+b.rho*b.g)*dz/2;}
 const pressureDrop=P.ps-rows.at(-1).p;
 const column={schema:'open-moon-column/1',key,world,planet:P,inputs:JSON.parse(JSON.stringify(spec)),rows,sample,
  summary:{name:spec.label,surfacePressure_Pa:P.ps,surfaceT_K:spec.surfaceT,surfaceRH:spec.rh,surfaceScaleHeight_m:K.Rd*rows[0].Tv/P.g0,dryLapse_K_per_km:1000*P.g0/K.cpd,lcl_m:sampleX(rows,lift.x).z,lclPressure_Pa:lift.p,lclTemperature_K:lift.T,
  cloudBase_m:base,cloudTop_m:top,lfc_m:lfc,equilibriumLevel_m:el,parcelReachesCloud:cloudSeen,positiveBuoyancyIntegral_J_per_kg:cape,preLFCNegativeWork_J_per_kg:cin,liquidWaterPath_kg_m2:liquidPath,iceWaterPath_kg_m2:icePath,inCloudOpticalDepth:tau,hydrostaticWeightRelativeError:Math.abs(weight-pressureDrop)/pressureDrop,pressureFloor_Pa:rows.at(-1).p,
  assumptions:['Prescribed pressure-coordinate sounding; no GCM or climate-equilibrium solution','Undiluted liquid pseudoadiabat for parcel buoyancy; mixed-phase optics is a diagnostic closure','Condensate retention and horizontal cloud morphology are selected inputs','Atmospheric study supplies no predicted rain to the surface-water model']}};
 return column;
}
function sampleRows(rows,z){
 finite(z,'Altitude');if(z<=rows[0].z)return {...rows[0]};if(z>=rows.at(-1).z)return {...rows.at(-1)};
 let lo=0,hi=rows.length-1;while(hi-lo>1){const m=(lo+hi)>>1;if(rows[m].z<=z)lo=m;else hi=m;}
 const a=rows[lo],b=rows[hi],t=(z-a.z)/(b.z-a.z),out={};for(const k of Object.keys(a))out[k]=mix(a[k],b[k],t);out.z=z;return out;
}
function sampleX(rows,x){let i=1;while(i<rows.length-1&&rows[i].x<x)i++;return {z:mix(rows[i-1].z,rows[i].z,(x-rows[i-1].x)/(rows[i].x-rows[i-1].x))};}
function textureData(column,n=256){
 const {cloudBase_m:base,cloudTop_m:top}=column.summary;if(base===null||top===null)return {data:new Float32Array(n*4),base:0,top:1,n};
 const data=new Float32Array(n*4);
 for(let i=0;i<n;i++){const z=mix(base,top,i/(n-1)),r=column.sample(z);data.set([r.extinction,r.ice/Math.max(1e-15,r.ice+r.liquid),r.windX,r.windZ],i*4);}
 return {data,base,top,n};
}
const API={K,EPS,PRESETS,planet,es,rs,virtualT,moistSlope,rkMoist,lcl,create,sampleRows,heightFromPotential,textureData};
root.OpenMoonWeatherColumn=API;if(typeof module!=='undefined')module.exports=API;
})(globalThis);
