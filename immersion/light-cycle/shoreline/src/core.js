/* Pure scenario / geometry / water-accounting functions. Metres, seconds, mm. */
(function(root){
'use strict';
const Landscape=root.OpenMoonLandscape||(typeof require==='function'?require('./landscape.js'):null);
if(!Landscape)throw Error('Load landscape.js before core.js');
const TAU=2*Math.PI, DAY=86400, PERIOD=29.53059*DAY;
const clamp=(x,a=0,b=1)=>Math.max(a,Math.min(b,x));
const mix=(a,b,t)=>a+(b-a)*t;
const smooth=(a,b,x)=>{const t=clamp((x-a)/(b-a));return t*t*(3-2*t);};
function rng(seed=18374){let s=seed>>>0;return()=>{s+=0x6D2B79F5;let t=Math.imul(s^s>>>15,1|s);t^=t+Math.imul(t^t>>>7,61|t);return((t^t>>>14)>>>0)/4294967296;};}
function hash(x,z){let n=Math.imul(x|0,374761393)+Math.imul(z|0,668265263);n=Math.imul(n^n>>>13,1274126177);return((n^n>>>16)>>>0)/4294967295;}
function noise(x,z){let i=Math.floor(x),j=Math.floor(z),u=x-i,v=z-j;u=u*u*(3-2*u);v=v*v*(3-2*v);return mix(mix(hash(i,j),hash(i+1,j),u),mix(hash(i,j+1),hash(i+1,j+1),u),v);}
function fbm(x,z){let a=.55,s=0;for(let i=0;i<5;i++){s+=a*noise(x,z);x=x*2.03+17.4;z=z*2.03+8.1;a*=.48;}return s;}
function shore(x){return Landscape.shorelineAtX(x);}
function surfaceHeight(x,z){return Landscape.height(x,z);}
function worldRadius(world='moon'){return world==='earth'?6371000:1737400;}
function curvatureSag(x,z,world='moon'){
 const R=worldRadius(world),r2=x*x+z*z;
 return r2/(R+Math.sqrt(Math.max(0,R*R-r2)));
}
function groundHeight(x,z,world='moon'){return surfaceHeight(x,z)-curvatureSag(x,z,world);}
function surfaceGradient(x,z,h=.2){return Landscape.gradient(x,z,h);}
function waterDepth(x,z){return Math.max(0,-surfaceHeight(x,z));}
const SHELTER={x:20,z:33,width:11,depth:8,roofHeight:4.3};
function roofMask(x,z,pad=0){return Math.abs(x-SHELTER.x)<SHELTER.width/2+pad&&Math.abs(z-SHELTER.z)<SHELTER.depth/2+pad?1:0;}
const WAYPOINTS=[{id:'shore',name:'Shoreline',x:0,z:9,yaw:0,pitch:-.04},{id:'shelter',name:'Rain shelter',x:20,z:33,yaw:.05,pitch:-.06},{id:'path',name:'Woodland path',x:-21,z:39,yaw:-.7,pitch:-.06},{id:'overlook',name:'High overlook',x:-56,z:47,yaw:.04,pitch:-.12}];
function pathDistance(x,z){let best=1e9;const a=[[0,-18],[0,9],[7,20],[20,33],[-6,43],[-30,43],[-56,47]];for(let i=1;i<a.length;i++){const [ax,az]=a[i-1],[bx,bz]=a[i],dx=bx-ax,dz=bz-az,t=clamp(((x-ax)*dx+(z-az)*dz)/(dx*dx+dz*dz));best=Math.min(best,Math.hypot(x-ax-t*dx,z-az-t*dz));}return best;}
function sunAt(phase){const h=((phase%1)+1)%1*TAU;return {x:Math.sin(h),y:Math.cos(h),z:0,elevation:Math.asin(Math.cos(h))*180/Math.PI,side:Math.sin(h)>=0?1:-1};}
function phaseForElevation(deg,rising=false){const h=Math.acos(clamp(Math.sin(deg*Math.PI/180),-1,1));return rising?1-h/TAU:h/TAU;}
function calibratedFov(screenHeightCm,distanceCm){if(!(screenHeightCm>0&&distanceCm>0))throw Error('Screen height and distance must be positive.');return 2*Math.atan(screenHeightCm/(2*distanceCm))*180/Math.PI;}
function cloudTau(lwc_g_m3,thickness_m,effectiveRadius_um){if(lwc_g_m3<0||thickness_m<=0||effectiveRadius_um<=0)throw Error('Invalid cloud input');return 3*(lwc_g_m3*.001)*thickness_m/(2*1000*effectiveRadius_um*1e-6);}
/* Authored four-hour episode. These are forcing inputs, not weather forecasts. */
const EPISODE=[
 [0, .10,.015,1100,850,1.5,.48,294,0,35000],
 [1800,.25,.040,1100,900,2.4,.59,293,0,30000],
 [3600,.52,.080,1050,1000,3.8,.72,292,0,22000],
 [5100,.88,.200,900,1100,5.4,.87,290,1.5,8000],
 [6300,.97,.290,800,1200,6.7,.96,289,7,3800],
 [8100,.99,.340,750,1300,6.2,.98,288.5,10,2500],
 [9600,.92,.190,900,1150,4.6,.94,289,3.5,5000],
 [10800,.70,.120,1100,1000,3.2,.87,290,0.3,9500],
 [12600,.37,.050,1350,950,2.6,.75,292,0,22000],
 [14400,.12,.020,1500,900,1.7,.64,293,0,35000]
];
function weatherAt(t,kind='episode'){
 let a,b;if(kind==='clear')a=b=[0,0,0,1100,850,1.5,.48,294,0,180000];else if(kind==='fog')a=b=[0,.62,.08,800,800,1.2,.99,288,0,110];else {t=clamp(t,0,14400);let j=1;while(j<EPISODE.length-1&&EPISODE[j][0]<t)j++;a=EPISODE[j-1];b=EPISODE[j];}
 const q=a===b?0:smooth(a[0],b[0],t),v=a.map((x,i)=>mix(x,b[i],q));
 return {coverage:v[1],lwc:v[2],cloudBase:v[3],thickness:v[4],wind:v[5],humidity:v[6],temperature:v[7],rain:v[8],visibility:v[9],radius:12,tau:cloudTau(v[2],v[4],12)};
}
function windDistance(t,kind='episode',step=10){let distance=0;for(let a=0;a<t;a+=step){const h=Math.min(step,t-a);distance+=weatherAt(a+h*.5,kind).wind*h;}return distance;}
/* A bounded linear reservoir. Exact solution for constant input and rates.
   dS/dt=I-(ke+kd)S; overflow removes inflow exceeding capacity. */
function reservoirStep(s,input,ke,kd,capacity,dt){
 if(![s,input,ke,kd,capacity,dt].every(Number.isFinite)||Math.min(s,input,ke,kd,capacity,dt)<0||s>capacity+1e-9)throw Error('Invalid reservoir input');
 const k=ke+kd,decay=Math.exp(-k*dt);let next,integral,overflow=0;
 if(k===0){next=Math.min(capacity,s+input*dt);integral=(s+next)*.5*dt;overflow=Math.max(0,s+input*dt-capacity);}
 else {const eq=input/k;let hit=Infinity;if(eq>capacity&&s<capacity)hit=-Math.log((capacity-eq)/(s-eq))/k;else if(s===capacity&&eq>=capacity)hit=0;
  const u=Math.min(dt,Math.max(0,hit));next=eq+(s-eq)*Math.exp(-k*u);integral=eq*u+(s-eq)*(-Math.expm1(-k*u))/k;
  if(u<dt){next=capacity;integral+=capacity*(dt-u);overflow=(input-k*capacity)*(dt-u);}}
 const evaporation=ke*integral,drainage=kd*integral;next=clamp(next,0,capacity);
 return {storage:next,evaporation,drainage,overflow:Math.max(0,overflow),input:input*dt,residual:s+input*dt-next-evaporation-drainage-overflow};
}
function emptyLedger(){return {exposed:0,canopy:0,leaf:0,sheltered:0,totalInput:0,evaporated:0,drained:0,residual:0};}
function stepLedger(l,w,dt){
 const evaporation=.00006*(1-w.humidity)*(1+w.wind*.2)*Math.exp((w.temperature-290)/25);
 const input=w.rain/3600;
 const leaf=reservoirStep(l.leaf,.45*input,evaporation*1.5,.00038,.35,dt);
 const canopy=reservoirStep(l.canopy,.55*input+(leaf.drainage+leaf.overflow)/Math.max(dt,1e-12),evaporation*.45,.00006,1.4,dt);
 const exposed=reservoirStep(l.exposed,input,evaporation,.00009,2.0,dt);
 const sheltered=reservoirStep(l.sheltered,0,evaporation*.4,.00006,1.4,dt);
 const next={exposed:exposed.storage,canopy:canopy.storage,leaf:leaf.storage,sheltered:sheltered.storage,
 totalInput:l.totalInput+2*input*dt,
 evaporated:l.evaporated+leaf.evaporation+canopy.evaporation+exposed.evaporation+sheltered.evaporation,
 drained:l.drained+canopy.drainage+canopy.overflow+exposed.drainage+exposed.overflow+sheltered.drainage+sheltered.overflow,residual:0};
 next.residual=next.totalInput-next.evaporated-next.drained-next.exposed-next.canopy-next.leaf-next.sheltered;
 return next;
}
function ledgerAt(t,kind='episode',step=10){let l=emptyLedger();for(let s=0;s<t;s+=step){const dt=Math.min(step,t-s);l=stepLedger(l,weatherAt(s+dt*.5,kind),dt);}return l;}
function dropletTerminalSpeed(radius=.0007,gravity=1.62,rhoAir=1.45){const mu=1.8e-5;let lo=0,hi=100;for(let i=0;i<70;i++){const v=(lo+hi)*.5,Re=2*radius*rhoAir*v/mu,Cd=Re<1000?24/Math.max(Re,1e-10)*(1+.15*Math.pow(Re,.687)):.44;const drag=.5*Cd*rhoAir*Math.PI*radius*radius*v*v,weight=(4/3)*Math.PI*radius**3*(1000-rhoAir)*gravity;if(drag<weight)lo=v;else hi=v;}return .5*(lo+hi);}
function waveOmega(k,depth=12,gravity=1.62){return Math.sqrt(gravity*k*Math.tanh(k*depth));}
function waveAt(x,z,t,wind=2,gravity=1.62){const dirs=[[.94,.342],[.36,.933],[-.45,.89],[.78,-.62]],ks=[.095,.22,.48,1.25],amps=[.15,.085,.04,.012];let y=0,nx=0,nz=0;for(let i=0;i<4;i++){const a=amps[i]*(.4+wind*.15),k=ks[i],q=k*(x*dirs[i][0]+z*dirs[i][1])-waveOmega(k,12,gravity)*t+i*1.5;y+=a*Math.sin(q);nx-=a*k*dirs[i][0]*Math.cos(q);nz-=a*k*dirs[i][1]*Math.cos(q);}return {y,nx,nz};}
const WATER_BANDS=[.095,.14,.22,.34,.48,.76,1.25,2.,3.5,5.8,9.5,15.];
function renderWaveAt(x,z,t,wind=1.5,gravity=1.62){
 let y=0,nx=0,nz=0;
 WATER_BANDS.forEach((k,i)=>{const theta=.95+Math.sin(i*2.399)*.71,dx=Math.cos(theta),dz=Math.sin(theta),a=.095*Math.pow(.095/k,1.12)*(.4+wind*.18),q=k*(x*dx+z*dz)-waveOmega(k,12,gravity)*t+i*2.721;y+=a*Math.sin(q);nx-=a*k*dx*Math.cos(q);nz-=a*k*dz*Math.cos(q);});
 return {y,nx,nz};
}
const API={WATER_BANDS,renderWaveAt,TAU,DAY,PERIOD,clamp,mix,smooth,rng,noise,fbm,shore,groundHeight,surfaceHeight,worldRadius,curvatureSag,surfaceGradient,waterDepth,pathDistance,SHELTER,roofMask,WAYPOINTS,sunAt,phaseForElevation,calibratedFov,cloudTau,EPISODE,weatherAt,windDistance,reservoirStep,emptyLedger,stepLedger,ledgerAt,dropletTerminalSpeed,waveOmega,waveAt};
root.OpenMoonCore=API;if(typeof module!=='undefined'&&module.exports)module.exports=API;
})(globalThis);
