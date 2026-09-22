/* CPU geometry-update microbenchmark. Both implementations use r186 and the
 * same unchanged landscape. This records no GPU time or display frame rate. */
'use strict';
const path=require('node:path'),fs=require('node:fs'),crypto=require('node:crypto');
const root=path.join(__dirname,'..'),file=path.resolve(process.argv[2]||path.join(root,'src/terrain.js'));
const T=require('../vendor/three-r186/three.module.js');globalThis.OM={};
const L=require('../src/landscape.js');require('../src/core.js');const {TerrainSystem}=require(file);
L.initialize();const scene=new T.Scene(),mat=new T.MeshBasicMaterial(),start=performance.now(),system=new TerrainSystem(T,scene,mat);const initialization_ms=performance.now()-start;
for(let i=0;i<24;i++)system.update(-10+i*.51,9+i*.05);
const times=[];for(let i=0;i<96;i++){const t=performance.now();system.update(i*.51,10+i*.07);times.push(performance.now()-t);}
const sorted=[...times].sort((a,b)=>a-b),bytes=fs.readFileSync(file),hash=crypto.createHash('sha256').update(bytes).digest('hex');
console.log(JSON.stringify({node:process.version,source_sha256:hash,initialization_ms,samples:times.length,median_ms:sorted[48],p95_ms:sorted[91],max_ms:sorted.at(-1),times_ms:times,scope:'CPU geometry and attribute updates only; r186 shared between implementations; no GPU benchmark'},null,2));system.dispose();mat.dispose();
