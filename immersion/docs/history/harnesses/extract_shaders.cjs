/* Extract our GLSL without a Three.js renderer. Constructors only collect text.
   This does not emulate, execute, or validate the Three.js rendering pipeline. */
const fs=require('node:fs'),path=require('node:path');
const ROOT=path.resolve(__dirname,'..');
const collected=[];
class Value {constructor(){this.texture={};this.image={data:new Float32Array(33*41*4)};this.position=this;this.up=this;this.matrixWorldInverse=this;this.projectionMatrix=this;this.matrix=this;} set(){return this;} copy(){return this;} multiply(){return this;} clone(){return new Value();} rotateX(){return this;} setAttribute(){return this;} setIndex(){return this;} add(){return this;} dispose(){} }
class ShaderMaterial extends Value {constructor(p){super();Object.assign(this,p);collected.push(p);}}
const T=new Proxy({ShaderMaterial,ShaderChunk:{lights_fragment_begin:'getDirectionalLightInfo( directionalLight, directLight );'}},{get(o,k){return o[k]||Value;}});
require('../src/core.js');require('../src/atmosphere.js');require('../src/water.js');require('../src/scene.js');
(async()=>{
 const atm=new OM.Atmosphere(T,{worlds:{}});await atm.init(new Value());
 OM.createWater(T,new Value(),atm,new Value());OM.createRain(T,new Value(),atm);
 const result={scope:'Custom GLSL extraction with constructor-only stand-ins. Integrated Three.js rendering remains untested.',programs:[]};
 const names=['sky','water','rain'];for(let i=0;i<collected.length;i++)result.programs.push({name:names[i],vertex:collected[i].vertexShader,fragment:collected[i].fragmentShader});
 for(const ground of [false,true]){
  const mat={userData:{}},state={materials:[]};OM.prepareMaterial(T,mat,atm,state,{sway:ground?0:.095,ground});
  const shader={uniforms:{},vertexShader:'void main(){\n#include <begin_vertex>\n#include <project_vertex>\n}',fragmentShader:'void main(){ vec4 diffuseColor=vec4(.5);\n#include <roughnessmap_fragment>\n#include <lights_fragment_begin>\n#include <lights_fragment_end>\ngl_FragColor=diffuseColor;\n#include <tonemapping_fragment>\n}'};
  mat.onBeforeCompile(shader);result.programs.push({name:ground?'terrain-material-additions':'foliage-material-additions',vertex:shader.vertexShader,fragment:shader.fragmentShader});
 }
 const out=path.join(ROOT,'validation/custom_shaders.json');fs.writeFileSync(out,JSON.stringify(result,null,2));console.log(out);
})().catch(e=>{console.error(e);process.exitCode=1;});
