/* World-space materials. Colour is linear; detail amplitudes are in metres.
 * Spatial derivatives suppress detail smaller than the pixel footprint.
 */
(function(root){
'use strict';
const OM=root.OM;
OM.SURFACE_GLSL=`
float omNoise(vec3 p);
float omHash2(vec2 p){uint n=uint(int(p.x))*374761393u+uint(int(p.y))*668265263u;n=(n^(n>>13u))*1274126177u;return float(n^(n>>16u))/4294967295.;}
float omN2(vec2 p){vec2 i=floor(p),f=fract(p);f=f*f*(3.-2.*f);return mix(mix(omHash2(i),omHash2(i+vec2(1,0)),f.x),mix(omHash2(i+vec2(0,1)),omHash2(i+vec2(1,1)),f.x),f.y);}
float omF2(vec2 p){return .56*omN2(p)+.27*omN2(p*2.07+19.)+.12*omN2(p*4.13+7.)+.05*omN2(p*8.3);}
float omBand(vec2 p,float frequency){float footprint=max(length(dFdx(p)),length(dFdy(p)))*frequency;return 1.-smoothstep(.35,1.4,footprint);}
float omGrain(vec2 p,float frequency){return (omN2(p*frequency)-.5)*omBand(p,frequency);}
vec3 omPerturb(vec3 p,vec3 n,float h){vec3 dx=dFdx(p),dy=dFdy(p);vec3 r1=cross(dy,n),r2=cross(n,dx);float det=dot(dx,r1);return normalize(abs(det)*n-sign(det)*(dFdx(h)*r1+dFdy(h)*r2));}
`;
OM.material=function(T,atm,state,kind,options={}){
 const mat=new T.MeshStandardMaterial({color:0xffffff,roughness:.9,...options});
 OM.prepareMaterial(T,mat,atm,state,{ground:kind==='terrain',sway:kind==='leaf'?.012:kind==='grass'?.008:0});
 const previous=mat.onBeforeCompile;
 mat.onBeforeCompile=shader=>{
  previous(shader);
  shader.fragmentShader=OM.SURFACE_GLSL+shader.fragmentShader;
  const terrain=`
   vec2 p=vOMWorld.xz;float elevation=vOMWorld.y+dot(p,p)/(2.*uR);
   float slope=1.-abs(normalize(cross(dFdx(vOMWorld),dFdy(vOMWorld))).y);
   float broad=omF2(p*.12),fine=omGrain(p,22.);
   float grassy=smoothstep(2.4,4.7,elevation+(broad-.5)*2.)*(1.-smoothstep(.035,.18,slope));
   float stone=smoothstep(.025,.15,slope)*smoothstep(.3,3.,elevation);
   vec3 sand=mix(vec3(.29,.265,.218),vec3(.43,.40,.32),broad);
   sand*=1.+fine*.30+omGrain(p,130.)*.15;
   float ridge=pow(.5+.5*sin(p.x*8.4+p.y*3.7+omF2(p*.55)*7.),5.)*omBand(p,1.6);
   sand*=1.-ridge*.055;
   vec3 loam=mix(vec3(.067,.081,.029),vec3(.14,.18,.059),broad);
   loam*=.82+.25*omN2(p*3.7);
   vec3 basalt=mix(vec3(.074,.075,.068),vec3(.16,.156,.133),omF2(p*.23));
   diffuseColor.rgb*=mix(mix(sand,loam,grassy),basalt,stone);
   diffuseColor.rgb*=mix(.58,1.,smoothstep(-.02,.55,elevation));
   float surfaceMicro=(.0025*omGrain(p,35.)+.0009*omGrain(p,115.)+.014*ridge)*(1.-stone)+stone*.018*omGrain(p,3.);
  `;
  const rock=`
   vec3 p=vOMWorld;vec2 q=p.xz+p.y*vec2(.53,.28);
   float broad=omF2(q*1.4),grain=omGrain(q,80.);
   float seam=pow(1.-abs(sin(p.y*9.+p.x*1.7+omF2(q*2.)*5.)),18.)*omBand(q,3.);
   vec3 stone=mix(vec3(.055,.060,.056),vec3(.17,.162,.143),broad);
   stone*=1.-.32*seam;stone*=1.+grain*.25;
   float fleck=smoothstep(.67,.84,omN2(q*140.))*omBand(q,140.);
   stone+=fleck*vec3(.07,.063,.045);
   float lichen=smoothstep(.63,.76,omF2(q*3.))*smoothstep(.4,.9,abs(normalize(cross(dFdx(vOMWorld),dFdy(vOMWorld))).y));
   stone=mix(stone,vec3(.21,.23,.102),lichen*.22);
   float elevation=p.y+dot(p.xz,p.xz)/(2.*uR);
   stone*=mix(.60,1.,smoothstep(-.02,.18,elevation));
   diffuseColor.rgb*=stone;float surfaceMicro=.003*omGrain(q,25.)-.003*seam+.0008*grain;
  `;
  const wood=`
   vec3 p=vOMWorld;float stripe=omN2(vec2(p.x*14.+p.z*8.,p.y*.75));
   diffuseColor.rgb*=mix(vec3(.045,.035,.023),vec3(.20,.15,.085),stripe);
   float surfaceMicro=(stripe-.5)*.018;
  `;
  const leaf=`
   float face=gl_FrontFacing?1.:.75;
   diffuseColor.rgb*=face;
   float surfaceMicro=0.;
  `;
  const detail=kind==='terrain'?terrain:kind==='rock'?rock:kind==='wood'?wood:leaf;
  shader.fragmentShader=shader.fragmentShader.replace('#include <color_fragment>','#include <color_fragment>\n'+detail);
  shader.fragmentShader=shader.fragmentShader.replace('#include <normal_fragment_maps>',`#include <normal_fragment_maps>\nnormal=omPerturb(-vViewPosition,normal,surfaceMicro);`);
  if(kind==='leaf'||kind==='grass'){
   shader.fragmentShader=shader.fragmentShader.replace('#include <lights_fragment_end>',`#include <lights_fragment_end>
    reflectedLight.indirectDiffuse+=diffuseColor.rgb*uDiffuse/OM_PI*.06;
   `);
  }
 };
 mat.customProgramCacheKey=()=>`open-moon-m1-${kind}`;
 return mat;
};
})(globalThis);
