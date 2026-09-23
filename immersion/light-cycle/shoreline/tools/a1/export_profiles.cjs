'use strict';
const fs=require('node:fs'),path=require('node:path'),crypto=require('node:crypto');
const W=require('../../src/weather-column.js'),A=require('../../src/atmospheric-profile.js'),O=require('../../src/profile-optics.js');
const root=path.resolve(__dirname,'../..'),out=path.resolve(process.argv[2]||path.join(root,'data/a1/generated'));
(async()=>{
 fs.mkdirSync(out,{recursive:true});const manifest={schema:'open-moon-a1-profile-manifest/1',profiles:[],sources:{}};
 for(const f of['src/weather-column.js','src/atmospheric-profile.js','src/profile-optics.js'])manifest.sources[f]=crypto.createHash('sha256').update(fs.readFileSync(path.join(root,f))).digest('hex');
 for(const world of['moon','earth','moon_no_ozone'])for(const key of['fair','convection','fog']){
  const c=W.create(key,world),p=A.fromColumn(c),id=world+'-'+key,filename=id+'.profile.json';fs.writeFileSync(path.join(out,filename),p.serialize());
  const refs=[];for(const[h,mu]of[[0,1],[2,0],[10000,.3],[60000,-.1]]){if(h>=p.state.upper.top_m)continue;const x=O.integrate(p,h,mu);refs.push({h_m:h,mu,...x,air_m:Number.isFinite(x.air_m)?x.air_m:null,ozone_DU:Number.isFinite(x.ozone_DU)?x.ozone_DU:null});}
  const surface=p.sample(0),vertical=refs[0];
  const alternatives=[];for(const opt of[{top_m:world==='earth'?180000:900000},{upperTemperature_K:200},{upperTemperature_K:250},{upperTemperature_K:300}]){const q=A.fromColumn(c,opt),v=O.integrate(q,0,1);alternatives.push({options:opt,profile_sha256:await q.fingerprint(),verticalAir_m:v.air_m,relativeVerticalChange:(v.air_m/vertical.air_m)-1});}
  manifest.profiles.push({id,file:filename,profile_sha256:await p.fingerprint(),world,regime:key,summary:c.summary,surface,vertical,independent_js_paths:refs,upperSensitivity:alternatives});
 }
 fs.writeFileSync(path.join(out,'profiles.json'),JSON.stringify(manifest,null,2)+'\n');console.log(JSON.stringify(manifest.profiles.map(p=>({id:p.id,air_m:p.vertical.air_m,surfaceAir:p.surface.airRelative,upper:p.upperSensitivity})),null,2));
})().catch(e=>{console.error(e);process.exitCode=1;});
