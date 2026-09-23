/* Independent adaptive path integration of the shared profile. No sky atlas,
 * exposure, cloud closure, or inherited exponential density enters this module. */
(function(root) {
'use strict';
const A = root.OpenMoonAtmosphericProfile || (typeof require === 'function' ? require('./atmospheric-profile.js') : null);
function rayInterval(profile, h, mu, length = Infinity) {
  if (!Number.isFinite(h) || h < 0 || h > profile.state.upper.top_m || !Number.isFinite(mu) || Math.abs(mu) > 1 || !(length >= 0)) throw new RangeError('Invalid optical ray');
  const R = profile.state.planet.radius, r = R + h, top = R + profile.state.upper.top_m;
  const tangent = r * Math.sqrt(Math.max(0, 1 - mu * mu));
  const end = -r * mu + Math.sqrt(Math.max(0, r * r * mu * mu + (top - r) * (top + r)));
  let ground = Infinity;
  if (mu < 0 && tangent < R) ground = -r * mu - Math.sqrt(Math.max(0, (R - tangent) * (R + tangent)));
  return {R, r, mu, end: Math.min(end, ground, length), blocked: ground < length && ground <= end, ground};
}
function integrate(profile, h, mu, options = {}) {
  const q = rayInterval(profile, h, mu, options.length_m ?? Infinity);
  if (q.blocked && options.length_m === undefined) return {air_m: Infinity, ozone_DU: Infinity, blocked: true, evaluations: 0, converged: true};
  const tol = options.relativeTolerance ?? 1e-7, maxDepth = options.maxDepth ?? 22;
  if (!(tol > 0 && tol < 1) || !Number.isInteger(maxDepth) || maxDepth < 1) throw new RangeError('Invalid integration tolerance');
  const ozone = profile.state.optics.ozone;
  let evaluations = 0, converged = true;
  const at = s => {
    evaluations++;
    const z = Math.max(0, Math.min(profile.state.upper.top_m, ((h + 2 * q.R) * h + 2 * q.r * mu * s + s * s) / (Math.sqrt(q.r * q.r + 2 * q.r * mu * s + s * s) + q.R)));
    return [profile.sample(z).airRelative, A.ozoneShape(z, ozone) * ozone.column_DU];
  };
  const simp = (a,b,fa,fm,fb) => fa.map((v,k) => (b-a)*(v+4*fm[k]+fb[k])/6);
  function recurse(a,b,fa,fm,fb,whole,eps,depth) {
    const mid = (a+b)/2, fl = at((a+mid)/2), fr = at((mid+b)/2);
    const left = simp(a,mid,fa,fl,fm), right = simp(mid,b,fm,fr,fb), sum = left.map((v,k) => v+right[k]);
    if (sum.every((v,k) => Math.abs(v-whole[k]) <= 15*eps[k])) return sum.map((v,k) => v+(v-whole[k])/15);
    if (!depth) { converged = false; return sum; }
    const l = recurse(a,mid,fa,fl,fm,left,eps.map(v=>v/2),depth-1), r = recurse(mid,b,fm,fr,fb,right,eps.map(v=>v/2),depth-1);
    return l.map((v,k)=>v+r[k]);
  }
  if (q.end === 0) return {air_m:0,ozone_DU:0,blocked:q.blocked,evaluations:0,converged:true};
  // Split at tangent, ozone corners and selected height shells so a thin layer
  // cannot be skipped by an initially empty Simpson stencil.
  const splits = [0,q.end]; if (-q.r*mu > 0 && -q.r*mu < q.end) splits.push(-q.r*mu);
  for (const z of [100,1000,5000,10000,25000,40000,80000,160000,320000,ozone.bottom_m,ozone.peak_m,ozone.top_m]) {
    const disc = q.r*q.r*mu*mu + (q.R+z-q.r)*(q.R+z+q.r);
    if (disc >= 0) for (const s of [-q.r*mu-Math.sqrt(disc),-q.r*mu+Math.sqrt(disc)]) if (s>0&&s<q.end) splits.push(s);
  }
  const pts = [...new Set(splits)].sort((a,b)=>a-b), result=[0,0];
  for (let i=1;i<pts.length;i++) { const a=pts[i-1],b=pts[i],fa=at(a),fm=at((a+b)/2),fb=at(b),whole=simp(a,b,fa,fm,fb), eps=[Math.max(1e-6,tol*Math.abs(whole[0])),Math.max(1e-9,tol*Math.abs(whole[1]))];
    const v=recurse(a,b,fa,fm,fb,whole,eps,maxDepth);v.forEach((x,k)=>result[k]+=x); }
  return {air_m:result[0],ozone_DU:result[1],blocked:q.blocked,evaluations,converged};
}
function transmission(profile,h,mu,wavelengths_nm,ozoneCrossSections_m2,options={}) {
  if (!Array.isArray(wavelengths_nm) || wavelengths_nm.length !== ozoneCrossSections_m2.length || wavelengths_nm.some(v=>!Number.isFinite(v)||v<=0) || ozoneCrossSections_m2.some(v=>!Number.isFinite(v)||v<0)) throw new RangeError('Invalid spectral coefficients');
  const c=integrate(profile,h,mu,options);
  if(!c.converged)throw new Error('Optical integration did not meet its tolerance');
  return {columns:c, values:wavelengths_nm.map((lambda,i)=>c.blocked?0:Math.exp(-1.24062e-6*(lambda/1000)**-4*c.air_m-ozoneCrossSections_m2[i]*A.DU*c.ozone_DU))};
}
const API={rayInterval,integrate,transmission};root.OpenMoonProfileOptics=API;if(typeof module!=='undefined')module.exports=API;
})(globalThis);
