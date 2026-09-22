/* Verified ESM dependency graph. Original release bytes are hashed before the
 * two declared import edges are linked to local Blob URLs. No eval/CJS shim.
 */
(function(root){
'use strict';const OM=root.OM=root.OM||{};
const REQUIRED=['three.core.js','three.module.js'];
const EDGES={'three.core.js':{},'three.module.js':{'./three.core.js':'three.core.js'},'three.webgpu.js':{'./three.core.js':'three.core.js'},'three.tsl.js':{'three/webgpu':'three.webgpu.js'}};
const hex=bytes=>Array.from(new Uint8Array(bytes),x=>x.toString(16).padStart(2,'0')).join('');
const decode=s=>Uint8Array.from(atob(s),c=>c.charCodeAt(0));
let loaded=null;
async function verifyModule(bytes,expected,name='module'){if(!expected||bytes.length!==expected.bytes||hex(await crypto.subtle.digest('SHA-256',bytes))!==expected.sha256)throw Error('Three.js integrity check failed: '+name);return true;}
async function loadThree(){
 if(loaded)return loaded;
 const pack=JSON.parse(document.getElementById('three-vendor').textContent);
 if(pack.manifest.version!=='0.186.0'||pack.manifest.revision!=='186')throw Error('Expected the pinned Three.js r186 ESM package.');
 if(!root.crypto?.subtle)throw Error('Web Crypto is required to verify the embedded dependency. Use a current browser or localhost.');
 const sources={};
 for(const [name,b64]of Object.entries(pack.modules)){
  if(!(name in EDGES))throw Error('Undeclared dependency module: '+name);
  const bytes=decode(b64),expected=pack.manifest.files[name];
  await verifyModule(bytes,expected,name);
  sources[name]=new TextDecoder('utf-8',{fatal:true}).decode(bytes);
 }
 for(const name of REQUIRED)if(!sources[name])throw Error('Missing dependency module: '+name);
 const urls={};
 try{
  for(const name of Object.keys(EDGES)){
   if(!sources[name])continue;
   let source=sources[name];
   for(const [specifier,dependency]of Object.entries(EDGES[name])){
    if(!urls[dependency])throw Error('Incomplete dependency graph: '+name);
    const literals=["'"+specifier+"'",'"'+specifier+'"'];
    if(!literals.some(l=>source.includes(l)))throw Error('Pinned import layout changed: '+name);
    for(const literal of literals)source=source.split(literal).join(JSON.stringify(urls[dependency]));
   }
   urls[name]=URL.createObjectURL(new Blob([source],{type:'text/javascript'}));
  }
  const webgl=await import(urls['three.module.js']);
  if(webgl.REVISION!=='186')throw Error('Loaded Three.js revision mismatch.');
  loaded={webgl,webgpu:urls['three.webgpu.js']?await import(urls['three.webgpu.js']):null,tsl:urls['three.tsl.js']?await import(urls['three.tsl.js']):null,manifest:pack.manifest};
  return loaded;
 }finally{for(const url of Object.values(urls))URL.revokeObjectURL(url);}
}
OM.loadThree=loadThree;
OM.downloadBlob=function(blob,name){const url=URL.createObjectURL(blob),a=document.createElement('a');a.href=url;a.download=name;a.click();setTimeout(()=>URL.revokeObjectURL(url),30000);};
OM.saveOffline=function(){const doc=document.documentElement.cloneNode(true);doc.querySelector('#boot')?.removeAttribute('hidden');doc.querySelector('#runtime-error')?.setAttribute('hidden','');doc.querySelector('#help')?.removeAttribute('open');doc.querySelector('#panel')?.classList.remove('open');doc.querySelector('body')?.classList.remove('minimal');OM.downloadBlob(new Blob(['<!doctype html>\n'+doc.outerHTML],{type:'text/html'}),'Open_Moon_Shoreline_r186_Offline.html');};
async function start(){
 const status=document.getElementById('boot-status');
 try{status.textContent='Verifying Three.js r186 modules…';const modules=await loadThree();root.THREE=modules.webgl;OM.dependency=modules.manifest;await (OM.bootLab?OM.bootLab(modules):OM.boot(modules.webgl));}
 catch(e){console.error(e);root.openMoonShorelineError=e.message;status.textContent=e.message;document.getElementById('boot-actions')?.removeAttribute('hidden');}
}
if(typeof document==='undefined'){if(typeof module!=='undefined')module.exports={REQUIRED,EDGES,verifyModule};return;}
const retry=document.getElementById('retry');if(retry)retry.onclick=()=>location.reload();
const about=document.getElementById('boot-about');if(about)about.onclick=()=>document.getElementById('help').showModal();
start();
})(globalThis);
