(function(root){'use strict';const OM=root.OM,C=root.OpenMoonCore;
function createWater(T,scene,atm,renderer){
 const reflectRT=new T.WebGLRenderTarget(512,512,{type:T.HalfFloatType,depthBuffer:true});reflectRT.texture.colorSpace=T.LinearSRGBColorSpace;
 const u={...atm.uniforms,uReflection:{value:reflectRT.texture},uMirrorMatrix:{value:new T.Matrix4()},uGravity:{value:1.62},uCamera:{value:new T.Vector3()},uRain:{value:0}};
 const wave=`
uniform float uTime,uWind,uGravity,uR;
float wo(float k){float e=exp(2.*k*12.);return sqrt(uGravity*k*(e-1.)/(e+1.));}
vec3 waterWave(vec2 p){vec3 v=vec3(0.);vec2 dirs[4];dirs[0]=vec2(.94,.342);dirs[1]=vec2(.36,.933);dirs[2]=vec2(-.45,.89);dirs[3]=vec2(.78,-.62);float ks[4];ks[0]=.095;ks[1]=.22;ks[2]=.48;ks[3]=1.25;float aa[4];aa[0]=.15;aa[1]=.085;aa[2]=.04;aa[3]=.012;for(int i=0;i<4;i++){float a=aa[i]*(.4+uWind*.15),k=ks[i],q=k*dot(p,dirs[i])-wo(k)*uTime+float(i)*1.5;v.x+=a*sin(q);v.yz-=a*k*dirs[i]*cos(q);}return v;}
`;
 const vs=wave+`uniform mat4 uMirrorMatrix;varying vec3 vW,vN;varying vec4 vMirror;void main(){vec3 p=position;vec3 w=waterWave(p.xz);p.y=w.x-dot(p.xz,p.xz)/(2.*uR);vN=normalize(vec3(w.y+p.x/uR,1.,w.z+p.z/uR));vW=p;vMirror=uMirrorMatrix*vec4(p,1.);gl_Position=projectionMatrix*viewMatrix*vec4(p,1.);}`;
 const fs=OM.SKY_UNIFORMS+OM.GLSL_NOISE+OM.CLOUDS+`uniform sampler2D uReflection;uniform vec3 uCamera;uniform float uRain;varying vec3 vW,vN;varying vec4 vMirror;
void main(){vec3 n=normalize(vN),v=normalize(uCamera-vW);float detail=omNoise(vec3(vW.xz*.75,uTime*.4));n=normalize(n+vec3((detail-.5)*.02,0.,(omNoise(vec3(vW.zx*.81,uTime*.33))-.5)*.02)*(1.+uRain*.04));float fres=.0204+.9796*pow(1.-max(0.,dot(n,v)),5.);vec2 uv=vMirror.xy/vMirror.w;uv+=n.xz*.008;vec3 reflected=texture2D(uReflection,clamp(uv,.001,.999)).rgb;
 float depth=max(0.,(-9.+4.*sin(vW.x*.035)+2.1*sin(vW.x*.083+.8)-vW.z)*.075);vec3 body=mix(vec3(.105,.16,.135),vec3(.013,.049,.048),1.-exp(-depth*.2))*(uDiffuse+uDirect*max(0.,uSun.y)*omCloudShadow(vW))/OM_PI;
 vec3 colour=mix(body,reflected,fres);vec3 h=normalize(v+uSun);float rough=.04+uWind*.006,a2=pow(rough,4.),nh=max(0.,dot(n,h)),den=nh*nh*(a2-1.)+1.;float D=a2/(OM_PI*den*den);float spec=min(6.,.0204*D/(4.*max(.05,dot(n,v))));colour+=uDirect*max(0.,dot(n,uSun))*spec*omCloudShadow(vW);
 float foam=(1.-smoothstep(.02,.75,depth))*smoothstep(.4,.88,sin(vW.z*2.6-uTime*1.1+sin(vW.x*.5))*.5+.5)*(.3+uWind*.045);colour=mix(colour,(uDiffuse+uDirect*max(uSun.y,0.)*.7)/OM_PI,foam);
 float dist=length(uCamera-vW),fog=1.-exp(-3.912/max(80.,uVisibility)*dist);colour=mix(colour,(uDiffuse+uDirect*max(0.,uSun.y)*.14)/OM_PI,fog);
 gl_FragColor=vec4(max(vec3(0.),colour),1.);\n#include <tonemapping_fragment>\n#include <colorspace_fragment>\n}`;
 // Fine local waves, progressively coarser water beyond the walking area.
 const verts=[],indices=[];
 function grid(x0,x1,z0,z1,nx,nz){const first=verts.length/3;for(let j=0;j<=nz;j++)for(let i=0;i<=nx;i++)verts.push(x0+(x1-x0)*i/nx,0,z0+(z1-z0)*j/nz);for(let j=0;j<nz;j++)for(let i=0;i<nx;i++){const a=first+j*(nx+1)+i,b=a+1,c=a+nx+1,d=c+1;indices.push(a,c,b,b,c,d);}}
 grid(-300,300,-300,300,256,256);
 grid(-13000,13000,-13000,-300,80,50);grid(-13000,13000,300,13000,80,50);
 grid(-13000,-300,-300,300,50,64);grid(300,13000,-300,300,50,64);
 const geo=new T.BufferGeometry();geo.setAttribute('position',new T.Float32BufferAttribute(verts,3));geo.setIndex(indices);
 const mat=new T.ShaderMaterial({uniforms:u,vertexShader:vs,fragmentShader:fs,side:T.FrontSide,toneMapped:true});const water=new T.Mesh(geo,mat);water.frustumCulled=false;scene.add(water);
 const mirrorCam=new T.PerspectiveCamera(),bias=new T.Matrix4().set(.5,0,0,.5,0,.5,0,.5,0,0,.5,.5,0,0,0,1),target=new T.Vector3();let last=-Infinity;
 function reflection(camera,now,force=false){u.uCamera.value.copy(camera.position);if(!force&&now-last<.075)return;last=now;mirrorCam.copy(camera);mirrorCam.position.y=-camera.position.y;camera.getWorldDirection(target);target.y=-target.y;mirrorCam.up.set(0,-1,0);mirrorCam.lookAt(mirrorCam.position.clone().add(target));mirrorCam.updateMatrixWorld();u.uMirrorMatrix.value.copy(bias).multiply(mirrorCam.projectionMatrix).multiply(mirrorCam.matrixWorldInverse);
  const rt=renderer.getRenderTarget(),clip=renderer.clippingPlanes,shadows=renderer.shadowMap.autoUpdate,tm=renderer.toneMapping,exp=renderer.toneMappingExposure;water.visible=false;renderer.shadowMap.autoUpdate=false;renderer.clippingPlanes=[new T.Plane(new T.Vector3(0,1,0),-.07)];renderer.toneMapping=T.NoToneMapping;renderer.toneMappingExposure=1;renderer.setRenderTarget(reflectRT);renderer.render(scene,mirrorCam);renderer.setRenderTarget(rt);renderer.clippingPlanes=clip;renderer.shadowMap.autoUpdate=shadows;renderer.toneMapping=tm;renderer.toneMappingExposure=exp;water.visible=true;
 }
 return {mesh:water,uniforms:u,reflection};
}
function createRain(T,scene,atm){const random=C.rng(719),positions=[],seeds=[];for(let i=0;i<5000;i++){const x=random(),y=random(),z=random();for(let j=0;j<2;j++){positions.push(x,y,z);seeds.push(j);}}
 const g=new T.BufferGeometry();g.setAttribute('position',new T.Float32BufferAttribute(positions,3));g.setAttribute('endPoint',new T.Float32BufferAttribute(seeds,1));const u={...atm.uniforms,uCamera:{value:new T.Vector3()},uFall:{value:1.56},uRain:{value:0},uRoofY:{value:8.24}};
 const vs=`uniform float uTime,uWind,uFall,uRain,uRoofY,uDrift;uniform vec3 uCamera;attribute float endPoint;varying float vAlpha;void main(){vec3 p=position;float x=fract(p.x+uDrift*.018),z=fract(p.z+uDrift*.005);float y=mod(p.y*24.-uTime*uFall,24.);vec3 w=vec3((x-.5)*50.+uCamera.x,y+uCamera.y-10.,(z-.5)*50.+uCamera.z);w-=endPoint*vec3(uWind*.018,-uFall/35.,uWind*.005);float protectedArea=(1.-step(6.,abs(w.x-20.)))*(1.-step(4.6,abs(w.z-33.)))*(1.-step(uRoofY,w.y));vAlpha=step(position.x,clamp(uRain/10.,0.,1.))*(1.-protectedArea)*smoothstep(-1.,1.,w.y)*.30;gl_Position=projectionMatrix*viewMatrix*vec4(w,1.);}`;
 const fs=`uniform vec3 uDiffuse,uDirect;uniform vec3 uSun;varying float vAlpha;void main(){vec3 c=(uDiffuse+uDirect*max(0.,uSun.y)*.1)*.18;gl_FragColor=vec4(c,vAlpha);\n#include <tonemapping_fragment>\n#include <colorspace_fragment>\n}`;
 const m=new T.ShaderMaterial({uniforms:u,vertexShader:vs,fragmentShader:fs,transparent:true,depthWrite:false});const rain=new T.LineSegments(g,m);rain.frustumCulled=false;scene.add(rain);return {mesh:rain,uniforms:u};
}
OM.createWater=createWater;OM.createRain=createRain;
})(globalThis);
