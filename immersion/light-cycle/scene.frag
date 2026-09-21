#version 300 es
precision highp float;precision highp sampler2DArray;
in vec2 uv;out vec4 fragColor;
uniform sampler2DArray sky;uniform sampler2D metadata,albedoMap,normalMap,geometryMap,shadowMap;
uniform float layer0,layer1,blend,exposureValue,sunElevation,aspect,fov,yaw,pitch,radius,sunSide,weatherTime,cloudCover,cloudTau,fogK,fogHeight,rainAmount,wetness,cameraY,quality;
uniform vec2 texSize;uniform vec3 rayBeta;
const float PI=3.141592653589793;
vec3 meta(int col,float layer){return texelFetch(metadata,ivec2(col,int(layer)),0).rgb;}
vec3 met(int col){return mix(meta(col,layer0),meta(col,layer1),blend);}
float hash(vec3 p){p=fract(p*.1031);p+=dot(p,p.yzx+33.33);return fract((p.x+p.y)*p.z);}
float noise3(vec3 p){vec3 i=floor(p),f=fract(p);f=f*f*(3.-2.*f);return mix(mix(mix(hash(i),hash(i+vec3(1,0,0)),f.x),mix(hash(i+vec3(0,1,0)),hash(i+vec3(1,1,0)),f.x),f.y),mix(mix(hash(i+vec3(0,0,1)),hash(i+vec3(1,0,1)),f.x),mix(hash(i+vec3(0,1,1)),hash(i+vec3(1,1,1)),f.x),f.y),f.z);}
float fbm(vec3 p){float s=0.,a=.53;for(int i=0;i<4;i++){s+=a*noise3(p);p=p*2.03+vec3(7.1,3.4,5.3);a*=.48;}return s;}
vec3 sun(){return vec3(sunSide*cos(sunElevation),sin(sunElevation),0.);}
vec3 radiance(vec3 d){float el=asin(clamp(d.y,-1.,1.));float v=.5+.5*sign(el)*sqrt(abs(el)/(PI*.5));float az=acos(clamp(d.x*sunSide/max(length(d.xz),1e-6),-1.,1.));vec2 tc=(vec2(az/PI,v)*(texSize-1.)+.5)/texSize;return mix(texture(sky,vec3(tc,layer0)).rgb*meta(0,layer0).x,texture(sky,vec3(tc,layer1)).rgb*meta(0,layer1).x,blend);}
vec3 skyE(vec3 n){float x=n.x*sunSide,y=n.z*sunSide,z=n.y;return max(vec3(0.),.282095*met(4)+.488603*y*met(5)+.488603*z*met(6)+.488603*x*met(7)+1.092548*x*y*met(8)+1.092548*y*z*met(9)+.315392*(3.*z*z-1.)*met(10)+1.092548*x*z*met(11)+.546274*(x*x-y*y)*met(12));}
float cloudT(){return 1./(1.+.75*.15*cloudTau);}
float beamT(){return exp(-cloudTau/max(sin(sunElevation),.035));}
vec3 diffuseWeather(){return met(3)*(1.-cloudCover+cloudCover*cloudT())+met(2)*cloudCover*max(0.,cloudT()-beamT());}
float cloudDensity(vec3 p){float h=(p.y+dot(p.xz,p.xz)/(2.*radius)-1500.)/2100.;if(h<0.||h>1.)return 0.;vec3 q=vec3((p.x+weatherTime*5.)*.00052,h*1.5,(p.z+weatherTime*1.2)*.00052);float n=fbm(q);float threshold=.80-cloudCover*.55;float shape=smoothstep(0.,.13,h)*(1.-smoothstep(.62,1.,h));return clamp((n-threshold)/.25,0.,1.)*shape;}
float cloudShadow(vec3 p){if(cloudCover<.001||cloudTau<.001||sunElevation<0.)return 1.;vec3 s=sun();float t=(2500.-p.y)/max(.025,s.y);vec3 q=p+s*t;float d=cloudDensity(q);return exp(-cloudTau*d/max(s.y,.08));}
float hg(float v,float g){return (1.-g*g)/(4.*PI*pow(max(.005,1.+g*g-2.*g*v),1.5));}
vec4 cloudLayer(vec3 d,bool cheap){if(cloudCover<.001||cloudTau<.001||d.y<-.075)return vec4(0.,0.,0.,1.);
float r=radius+cameraY,b=r*d.y;float start=-b+sqrt(max(0.,b*b+2.*radius*(1500.-cameraY)+1500.*1500.-cameraY*cameraY));float end=-b+sqrt(max(0.,b*b+2.*radius*(3600.-cameraY)+3600.*3600.-cameraY*cameraY));
end=min(end,200000.);if(end<=start)return vec4(0.,0.,0.,1.);
int count=cheap?16:(quality<.5?16:quality>1.5?64:32);float ds=(end-start)/float(count);vec3 s=sun(),C=vec3(0.);float T=1.;vec3 ambient=max(met(14),vec3(0.))/PI;vec3 sunlight=max(met(13),vec3(0.));float phase=.65*hg(dot(d,s),.62)+.35*hg(dot(d,s),-.15);
for(int i=0;i<64;i++){if(i>=count||T<.007)break;vec3 p=vec3(0,cameraY,0)+d*(start+(float(i)+.5)*ds);float den=cloudDensity(p);if(den<.001)continue;float h=(p.y+dot(p.xz,p.xz)/(2.*radius)-1500.)/2100.;float lightDepth=0.;if(sunlight.x+sunlight.y+sunlight.z>1e-7){lightDepth=(cloudDensity(p+s*150.)*180.+cloudDensity(p+s*500.)*450.+cloudDensity(p+s*1200.)*900.)*cloudTau/2100.;}
float dt=den*cloudTau*ds/2100.;float a=1.-exp(-dt);vec3 source=ambient*(.35+.4*h)+sunlight*phase*exp(-lightDepth);C+=T*a*source;T*=1.-a;}
return vec4(C,T);}
vec3 fogColour(){vec3 c=(diffuseWeather()+met(2)*.22)/PI;float L=dot(c,vec3(.2126,.7152,.0722));return mix(c,vec3(L),.25);}
vec3 skyColour(vec3 d,bool cheap){vec3 light=radiance(d);float sd=acos(clamp(dot(d,sun()),-1.,1.));float pixel=max(fwidth(sd),.00002);float cover=1.-smoothstep(.004654211-pixel,.004654211+pixel,sd);float dep=acos(radius/(radius+1.7));float f=clamp((sunElevation+dep)/.004654211,-1.,1.);float visible=(acos(-f)+f*sqrt(max(0.,1.-f*f)))/PI;
if(d.y>-.002&&visible>1e-6){vec3 disk=max(met(1),vec3(0.))/(PI*.004654211*.004654211*visible);light+=disk*cover;float halo=exp(-sd*sd/.00007)*.00008;light+=disk*halo;}
vec4 cc=cloudLayer(d,cheap);light=light*cc.a+cc.rgb;if(fogK>0.){float optical=fogK*fogHeight/max(.025,d.y);float tr=exp(-optical);light=mix(fogColour(),light,tr);}return max(light,vec3(0.));}
vec3 displayMap(vec3 c){float Y=max(dot(c,vec3(.2126,.7152,.0722)),0.);float mn=min(c.x,min(c.y,c.z));if(mn<0.)c=mix(vec3(Y),c,clamp(Y/(Y-mn+1e-20),0.,1.));c=max(c,vec3(0.))*exposureValue;float L=dot(c,vec3(.2126,.7152,.0722));c/=1.+L;float mx=max(c.x,max(c.y,c.z));if(mx>1.){float l=dot(c,vec3(.2126,.7152,.0722));c=mix(vec3(l),c,clamp((1.-l)/(mx-l+1e-8),0.,1.));}return mix(12.92*c,1.055*pow(max(c,vec3(0.)),vec3(1./2.4))-.055,greaterThan(c,vec3(.0031308)));}
void main(){vec2 v=uv*2.-1.;vec3 forward=vec3(cos(yaw)*cos(pitch),sin(pitch),sin(yaw)*cos(pitch));vec3 right=vec3(-sin(yaw),0,cos(yaw));vec3 up=vec3(-cos(yaw)*sin(pitch),cos(pitch),-sin(yaw)*sin(pitch));vec3 d=normalize(forward+v.x*aspect*tan(fov*.5)*right+v.y*tan(fov*.5)*up);
vec2 tc=vec2(atan(d.z,d.x)/(2.*PI)+.5,acos(clamp(d.y,-1.,1.))/PI);vec3 geo=texture(geometryMap,tc).rgb;int kind=int(floor(geo.b*255.+.5));vec3 light;
if(kind==0)light=skyColour(d,false);else{
 float encoded=floor(geo.r*255.+.5)*256.+floor(geo.g*255.+.5);float distance=exp2(encoded/65535.*16.)-1.;vec3 p=vec3(0,cameraY,0)+d*distance;
 vec3 n=normalize(texture(normalMap,tc).rgb*2.-1.),col=pow(texture(albedoMap,tc).rgb,vec3(2.));vec3 ss=texture(shadowMap,tc).rgb;float horizon=(sunSide>0.?ss.r:ss.g)*PI*.5;float visibility=smoothstep(horizon-.015,horizon+.015,sunElevation)*cloudShadow(p);
 float directFog=exp(-fogK*fogHeight/max(.025,sin(sunElevation)));vec3 D=max(met(1),vec3(0.))*visibility*directFog;
 vec3 E=skyE(n)*(1.-cloudCover+cloudCover*cloudT())+met(2)*cloudCover*max(0.,cloudT()-beamT())*(.5+.5*n.y);E=max(E,vec3(0.))*ss.b;
 if(kind==2){
  float g=radius<3000000.?1.62:9.80665;float t=weatherTime;float phase0=p.x*.6+p.z*.14-sqrt(g*.616)*t,phase1=p.x*1.35-p.z*.75-sqrt(g*1.544)*t,phase2=p.x*3.4+p.z*2.2-sqrt(g*4.05)*t;
  float footprint=max(length(dFdx(p)),length(dFdy(p)));vec3 wf=exp(-vec3(.616* .616,1.544*1.544,4.05*4.05)*footprint*footprint*.2);n=normalize(vec3(.07*cos(phase0)*wf.x+.035*cos(phase1)*wf.y+.018*cos(phase2)*wf.z,1.,.035*cos(phase0)*wf.x-.04*cos(phase1)*wf.y+.018*sin(phase2)*wf.z));
  vec3 refl=reflect(d,n);refl.y=max(.001,refl.y);vec3 reflected=skyColour(normalize(refl),true);float fres=.02+.98*pow(1.-max(0.,dot(-d,n)),5.);
  float shallow=exp(-distance/160.);vec3 waterCol=mix(vec3(.01,.043,.05),vec3(.028,.09,.075),shallow);
  light=reflected*fres+diffuseWeather()*waterCol*(1.-fres)/PI;
  vec3 H=normalize(sun()-d);float nv=max(.04,dot(n,-d)),nh=max(0.,dot(n,H));float a=.06;float spec=a*a/(PI*pow(max(.003,nh*nh*(a*a-1.)+1.),2.));light+=D*.018*spec/(4.*nv);
 }else{
  col*=1.-wetness*.32;light=col*(E+D*max(dot(n,sun()),0.))/PI;
  if(kind==3)light+=col*max(met(1),vec3(0.))*.055*visibility*pow(max(dot(d,sun()),0.),4.)/PI;
  if(wetness>.01){vec3 refl=reflect(d,n);refl.y=abs(refl.y);float fr=.035+.2*pow(1.-max(dot(-d,n),0.),5.);light+=wetness*fr*radiance(refl);}
 }
 vec3 molecular=exp(-rayBeta*distance);light=light*molecular+max(radiance(d),vec3(0.))*(1.-molecular);
 if(fogK>0.)light=mix(fogColour(),light,exp(-fogK*distance));
}
// Screen-space rain is an illustrative exposure effect, driven by the weather clock.
if(rainAmount>.001){float rain=0.;for(int i=0;i<3;i++){float size=45.+float(i)*23.;vec2 q=vec2(uv.x*aspect,uv.y)*size;float column=floor(q.x+q.y*.17);float seed=hash(vec3(column,float(i),7.));q.y+=weatherTime*(10.+seed*12.);float x=fract(q.x+q.y*.17)-.5;float y=fract(q.y+seed*25.);rain+=(1.-smoothstep(.005,.035,abs(x)))*smoothstep(.0,.03,y)*(1.-smoothstep(.15,.27,y))*step(.3,seed);}light+=max(diffuseWeather(),vec3(0.))/PI*rainAmount*rain*.055;}
vec3 outc=displayMap(light);float vignette=1.-.035*dot(v,v);outc*=vignette;float grain=(hash(vec3(gl_FragCoord.xy,1.))-.5)/255.;fragColor=vec4(clamp(outc+grain,0.,1.),1.);}
