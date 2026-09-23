/* Solver-neutral frozen Cartesian cloud field. Vertex samples, x-fastest.
 * The field's voxelization error is independent of either solver's error.
 * Shared optical assumptions: gray extinction, selectable single-scattering
 * albedo, and the revision-04 liquid/ice two-HG phase mixture. */
(function(root){
'use strict';
function decode(s){const b=typeof Buffer!=='undefined'?new Uint8Array(Buffer.from(s,'base64')):Uint8Array.from(atob(s),x=>x.charCodeAt(0));if(b.byteLength%4)throw new RangeError('Invalid float32 byte length');const out=new Float32Array(b.length/4),v=new DataView(b.buffer,b.byteOffset,b.byteLength);for(let i=0;i<out.length;i++)out[i]=v.getFloat32(i*4,true);return out;}
function create(spec){
 if(spec?.schema!=='open-moon-frozen-cloud/1')throw new TypeError('Unsupported frozen cloud schema');
 spec=JSON.parse(JSON.stringify(spec));
 const freeze=o=>{if(o&&typeof o==='object'){Object.values(o).forEach(freeze);Object.freeze(o);}};freeze(spec);
 const dims=spec.dimensions,lo=spec.bounds_min_m,hi=spec.bounds_max_m;
 if(!Array.isArray(dims)||dims.length!==3||dims.some(n=>!Number.isInteger(n)||n<2)||!Array.isArray(lo)||!Array.isArray(hi)||lo.length!==3||hi.length!==3||lo.some((v,i)=>!Number.isFinite(v)||!Number.isFinite(hi[i])||hi[i]<=v))throw new RangeError('Invalid frozen grid');
 const ext=decode(spec.extinction_float32le),ice=decode(spec.ice_fraction_float32le),n=dims.reduce((a,b)=>a*b,1);
 if(ext.length!==n||ice.length!==n||ext.some(x=>!Number.isFinite(x)||x<0)||ice.some(x=>!Number.isFinite(x)||x<0||x>1))throw new RangeError('Invalid optical voxel values');
 const albedo=spec.single_scattering_albedo??1;if(!Number.isFinite(albedo)||albedo<0||albedo>1)throw new RangeError('Invalid single scattering albedo');
 let maximum=0;for(const x of ext)maximum=Math.max(maximum,x);
 function sample(p){
  if(!Array.isArray(p)||p.length!==3||p.some(v=>!Number.isFinite(v)))throw new RangeError('Invalid sampling point');
  if(p.some((v,i)=>v<lo[i]||v>hi[i]))return {extinction:0,iceFraction:0};
  const q=p.map((v,i)=>(v-lo[i])/(hi[i]-lo[i])*(dims[i]-1)),a=q.map((v,i)=>Math.min(dims[i]-2,Math.floor(v))),f=q.map((v,i)=>v-a[i]);let e=0,g=0;
  for(let z=0;z<2;z++)for(let y=0;y<2;y++)for(let x=0;x<2;x++){const w=(x?f[0]:1-f[0])*(y?f[1]:1-f[1])*(z?f[2]:1-f[2]),k=((a[2]+z)*dims[1]+a[1]+y)*dims[0]+a[0]+x;e+=w*ext[k];g+=w*ice[k];}
  return {extinction:e,iceFraction:g};
 }
 return Object.freeze({spec,get extinction(){return new Float32Array(ext);},get iceFraction(){return new Float32Array(ice);},sample,majorant:maximum*(1+1e-6),singleScatteringAlbedo:albedo});
}
const GLSL=`
uniform highp sampler3D a1FrozenExt;
uniform highp sampler3D a1FrozenIce;
uniform vec3 a1FrozenMin,a1FrozenMax,a1FrozenSize;
vec2 a1FrozenSample(vec3 p){
 if(any(lessThan(p,a1FrozenMin))||any(greaterThan(p,a1FrozenMax)))return vec2(0.);
 vec3 uv=((p-a1FrozenMin)/(a1FrozenMax-a1FrozenMin)*(a1FrozenSize-1.)+.5)/a1FrozenSize;
 return vec2(texture(a1FrozenExt,uv).r,texture(a1FrozenIce,uv).r);
}
vec2 a1Box(vec3 p,vec3 d){
 float entry=0.,exit=1e30;
 for(int k=0;k<3;k++){
  if(abs(d[k])<1e-12){if(p[k]<a1FrozenMin[k]||p[k]>a1FrozenMax[k])return vec2(1.,0.);}
  else {float a=(a1FrozenMin[k]-p[k])/d[k],b=(a1FrozenMax[k]-p[k])/d[k];entry=max(entry,min(a,b));exit=min(exit,max(a,b));}
 }
 return vec2(entry,exit);
}
float a1Phase(float mu,float ice){float g=mix(.78,.82,ice),back=-.25;return .85*(1.-g*g)/(12.566370614359*pow(1.+g*g-2.*g*mu,1.5))+.15*(1.-back*back)/(12.566370614359*pow(1.+back*back-2.*back*mu,1.5));}
`;
const API={create,decode,GLSL};root.OpenMoonFrozenCloud=API;if(typeof module!=='undefined')module.exports=API;
})(globalThis);
