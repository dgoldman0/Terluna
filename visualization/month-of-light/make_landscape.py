"""Deterministic illustrative cove. Pretrace fixed-viewpoint geometry for both renderers.
Geometry and biology are design choices, not a lunar terrain survey or ecosystem model.
All geometry is local, in metres. No climate or cloud statistics are inferred.
"""
from pathlib import Path
import numpy as np, math, json, time
from numba import njit, prange
from PIL import Image
ROOT=Path(__file__).resolve().parent

@njit(cache=True)
def fract(x):return x-math.floor(x)
@njit(cache=True)
def hash2(x,z):return fract(math.sin(x*127.1+z*311.7)*43758.5453123)
@njit(cache=True)
def noise(x,z):
 i=math.floor(x);j=math.floor(z);u=x-i;v=z-j;u=u*u*(3-2*u);v=v*v*(3-2*v)
 return ((1-u)*hash2(i,j)+u*hash2(i+1,j))*(1-v)+((1-u)*hash2(i,j+1)+u*hash2(i+1,j+1))*v
@njit(cache=True)
def fbm(x,z):
 a=.52;s=0.
 for k in range(5):s+=a*noise(x,z);x=x*2.03+12.2;z=z*2.03+9.5;a*=.49
 return s
@njit(cache=True)
def terrain(x,z):
 # A wooded headland around an open bay, with an archipelago beyond it.
 coast=70+20*math.sin(z*.018)+abs(z)*.17
 land=1/(1+math.exp(max(-50.,min(50.,(x-coast)*.045))))
 local=-13+(14+3*(fbm(x*.022,z*.022)-.5))*land
 # Coastline hills rise off to either side, remaining out of the central water view.
 upland=max(0.,1-math.exp(-max(abs(z)-150.,0.)/600.))
 local+=land*upland*(50+100*fbm(x*.002,z*.002))
 # Several ridged peaks. Their positions are deliberately illustrative.
 mountains=0.
 for cx,cz,sx,sz,h in ((3600.,-1950.,1100.,1050.,530.),(6200.,-700.,1350.,1800.,790.),(7800.,2300.,1900.,1400.,1080.),(3400.,3100.,1200.,1600.,420.),(11200.,-2600.,1900.,2200.,1300.)):
  q=((x-cx)/sx)**2+((z-cz)/sz)**2
  if q<13:mountains+=h*math.exp(-q)*(.48+.78*fbm(x*.0023+cx,z*.0023+cz))
 return local+mountains

# Primitive: centre xyz, halfsize xyz, albedo rgb, kind, axis xyz, half length.
# Kind 3 foliage, 4 bark, 5 stone. Shapes are ellipsoids; trunks are oriented cigars.
def scene():
 rng=np.random.default_rng(7124);P=[];T=[]
 def ell(c,r,col,k,axis=(0,1,0)):
  P.append([*c,*r,*col,k,*axis])
 def branch(a,b,r,col):
  a=np.array(a);b=np.array(b);d=b-a;le=np.linalg.norm(d);ell((a+b)/2,(r,le/2+r,r),col,4,d/le)
 def tree(x,z,h,w):
  y=terrain(x,z);T.append([x,y,z,h,w]);col=np.array([.11,.069,.037])
  branch((x,y,z),(x+.35,y+h*.70,z+.1),h*.022,col)
  for k in range(10):
   ang=k*2.39996+rng.uniform(-.2,.2);rr=w*rng.uniform(.40,.95);cy=y+h*rng.uniform(.64,.94)
   center=np.array([x+math.cos(ang)*rr,cy,z+math.sin(ang)*rr]);r=w*rng.uniform(.28,.46)
   branch((x,y+h*.35,z),center,h*.0075,col*rng.uniform(.85,1.15))
   ell(center,(r,r*rng.uniform(.8,1.3),r),np.array([.055,.125,.017])*rng.uniform(.7,1.4),3)
   # Sub-clusters produce a irregular crown outline at the viewer's scale.
   for j in range(13):
    v=rng.normal(size=3);v/=np.linalg.norm(v);c=center+v*r*.87;s=r*rng.uniform(.20,.36)
    ell(c,(s,s*rng.uniform(.7,1.4),s),np.array([.060,.14,.018])*rng.uniform(.72,1.3),3)
 near=[(10,-16,13,4.3),(-7,12,15,4.4),(34,33,11,3.6),(51,-32,12,3.9),(2,-42,17,5.1),(-28,-28,18,5.5),(25,62,15,4.8),(-28,44,14,4.2)]
 for x,z,h,w in near:tree(x,z,h,w)
 for _ in range(58):
  x=rng.uniform(-240,85);z=rng.uniform(-230,230)
  if x*x+z*z<2500 or terrain(x,z)<-.1:continue
  h=rng.uniform(7,18);tree(x,z,h,h*rng.uniform(.22,.32))
 for _ in range(70):
  x=rng.uniform(15,85);z=rng.uniform(-120,120);y=terrain(x,z)
  if y< -4:continue
  r=rng.uniform(.2,1.8);c=[x,y+r*.10,z];col=np.array([.20,.22,.20])*rng.uniform(.7,1.4)
  ell(c,(r,r*.65,r*.75),col,5)
 # Low shrubs near the path; their size is illustrative.
 for _ in range(60):
  x=rng.uniform(-12,40);z=rng.uniform(-30,30)
  if abs(z)<2+max(0,x)*.018:continue
  y=terrain(x,z);s=rng.uniform(.15,.45)
  ell((x,y+s*.6,z),(s,s*.7,s),np.array([.06,.13,.02])*rng.uniform(.7,1.3),3)
 return np.asarray(P,np.float64),np.asarray(T,np.float64)

def bvh(P):
 # Rotation-independent bounding spheres keep box intersection conservative.
 lo=P[:,:3]-P[:,3:6].max(axis=1)[:,None];hi=P[:,:3]+P[:,3:6].max(axis=1)[:,None]
 nodes=[];order=[]
 def rec(ids):
  n=len(nodes);mn=lo[ids].min(0);mx=hi[ids].max(0);nodes.append([*mn,*mx,-1,-1,-1,-1])
  if len(ids)<=7:
   start=len(order);order.extend(ids);nodes[n][8:]=[start,len(ids)]
  else:
   axis=np.argmax(np.ptp(P[ids,:3],axis=0));ids=ids[np.argsort(P[ids,axis])];m=len(ids)//2
   a=rec(ids[:m]);b=rec(ids[m:]);nodes[n][6:8]=[a,b]
  return n
 rec(np.arange(len(P)));return np.array(nodes),np.array(order,np.int32)

@njit(cache=True)
def boxhit(ro,rd,b,limit):
 low=0.;high=limit
 for i in range(3):
  inv=1/(rd[i]+1e-30);a=(b[i]-ro[i])*inv;c=(b[i+3]-ro[i])*inv
  low=max(low,min(a,c));high=min(high,max(a,c))
  if high<low:return False
 return True

@njit(cache=True)
def intersect(ro,rd,P,B,order,limit):
 stack=np.empty(96,np.int32);stack[0]=0;n=1;best=limit;bid=-1
 while n:
  n-=1;idx=stack[n];bb=B[idx]
  if not boxhit(ro,rd,bb,best):continue
  if bb[8]<0:
   stack[n]=int(bb[6]);stack[n+1]=int(bb[7]);n+=2;continue
  for k in range(int(bb[8]),int(bb[8]+bb[9])):
   ip=order[k];p=P[ip];ax=p[10:13];oc=ro-p[:3]
   # Implicit oriented spheroid; the x/z radii match for foliage and trunk.
   if p[9]==4:
    oy=np.dot(oc,ax);dy=np.dot(rd,ax);rr=p[3];ry=p[4]
    A=(np.dot(rd,rd)-dy*dy)/(rr*rr)+dy*dy/(ry*ry)
    Q=(np.dot(oc,rd)-oy*dy)/(rr*rr)+oy*dy/(ry*ry)
    C=(np.dot(oc,oc)-oy*oy)/(rr*rr)+oy*oy/(ry*ry)-1
   else:
    oo=oc/p[3:6];dd=rd/p[3:6];A=np.dot(dd,dd);Q=np.dot(oo,dd);C=np.dot(oo,oo)-1
   disc=Q*Q-A*C
   if disc<=0:continue
   t=(-Q-math.sqrt(disc))/A
   if t>.01 and t<best:best=t;bid=ip
 return best,bid

@njit(cache=True)
def groundhit(ro,rd):
 # Conservative-enough practical height-field tracing, verified visually for this cove.
 t=.1;old=.1
 for k in range(2600):
  x=ro[0]+rd[0]*t;z=ro[2]+rd[2]*t;y=ro[1]+rd[1]*t
  dh=y-terrain(x,z)
  if dh<.03:
   lo=old;hi=t
   for j in range(7):
    mid=(lo+hi)/2
    if ro[1]+rd[1]*mid>terrain(ro[0]+rd[0]*mid,ro[2]+rd[2]*mid):lo=mid
    else:hi=mid
   return (lo+hi)/2
  old=t;t+=max(.15,min(dh*.25,max(5.,t*.022)))
  if t>19000 or (rd[1]>0 and ro[1]+rd[1]*t>1600):break
 return 1e8

@njit(cache=True)
def horizon(p,T,sign):
 # Visual soft-shadow proxy. Canopy silhouettes are treated as opaque obstacles.
 h=0.
 for j in range(len(T)):
  dx=(T[j,0]-p[0])*sign
  if dx<.2:continue
  dz=abs(T[j,2]-p[2]);w=T[j,4]
  if dz>w:continue
  top=T[j,1]+T[j,3]-.1*T[j,3]*(dz/w)**2
  h=max(h,math.atan2(top-p[1],dx))
 for dx in (12.,30.,65.,150.,350.,700.,1300.,2400.,4000.,7000.):
  h=max(h,math.atan2(terrain(p[0]+dx*sign,p[2])-p[1],dx))
 return h

@njit(parallel=True,cache=True)
def render(W,H,P,B,O,T):
 alb=np.zeros((H,W,3),np.uint8);norm=np.zeros_like(alb);info=np.zeros_like(alb);shad=np.zeros_like(alb)
 ro=np.array([0.,terrain(0.,0.)+1.7,0.])
 for j in prange(H):
  elev=math.pi/2-(j+.5)/H*math.pi;ce=math.cos(elev);se=math.sin(elev)
  for i in range(W):
   az=(i+.5)/W*2*math.pi-math.pi;rd=np.array([ce*math.cos(az),se,ce*math.sin(az)])
   if elev>1.1:continue
   t=groundhit(ro,rd) if se<.32 else 1e8;kind=1
   if se<0:
    wt=(-3.-ro[1])/se
    if wt>0 and wt<t and wt<14000:t=wt;kind=2
   tt,idx=intersect(ro,rd,P,B,O,t)
   if idx>=0:t=tt;kind=int(P[idx,9])
   if t>19000:continue
   p=ro+rd*t
   ao=1.;col=np.zeros(3)
   if idx>=0:
    pp=P[idx];v=p-pp[:3];ax=pp[10:13]
    if kind==4:
     vv=np.dot(v,ax);n=(v-ax*vv)/(pp[3]**2)+ax*vv/(pp[4]**2)
    else:n=v/(pp[3:6]**2)
    n=n/np.sqrt(np.dot(n,n));col=pp[6:9].copy()
    if kind==3:
     grain=noise(p[0]*35+p[1]*10,p[2]*35+p[1]*8)
     col*=.58+.70*grain;ao=.56+.34*max(n[1],0.)
    elif kind==4:
     col*=.6+.6*noise(p[0]*50+p[2]*50,p[1]*1.6);ao=.8
    else:
     col*=.55+.70*fbm(p[0]*4+p[1],p[2]*4)
     bump=np.array([noise(p[0]*8,p[1]*8),noise(p[1]*7,p[2]*7),noise(p[2]*8,p[0]*8)])-.5
     n+=bump*.23;n/=math.sqrt(np.dot(n,n))
   elif kind==2:
    n=np.array([0.,1.,0.]);col=np.array([.015,.08,.075])
   else:
    ep=max(.04,min(2.,t*.0004));n=np.array([terrain(p[0]-ep,p[2])-terrain(p[0]+ep,p[2]),2*ep,terrain(p[0],p[2]-ep)-terrain(p[0],p[2]+ep)])
    n/=math.sqrt(np.dot(n,n));v=fbm(p[0]*.10,p[2]*.10)
    sand=max(0.,min(1.,(-p[1]-.20)*.38))
    stone=max(0.,min(1.,(1-n[1])/.55))
    grass=np.array([.08,.145,.028])*(.68+.7*v)
    earth=np.array([.30,.25,.16])*(.75+.4*v)
    rock=np.array([.18,.195,.18])*(.65+.5*v)
    col=grass*(1-sand)+earth*sand
    col=col*(1-stone)+rock*stone
    # A narrow footpath is an intentionally familiar scale cue.
    pathz=1.2*math.sin(p[0]*.075)
    if p[0]>-9 and p[0]<59 and abs(p[2]-pathz)<.68+p[0]*.005:
     blend=max(0.,min(1.,(.75+p[0]*.005-abs(p[2]-pathz))/.18));col=col*(1-blend)+np.array([.21,.175,.115])*blend
    grain=noise(p[0]*27,p[2]*27);col*=.80+.38*grain
    if t<100:
     n[0]+=(noise(p[0]*8,p[2]*8)-.5)*.12;n[2]+=(noise(p[0]*9+77,p[2]*9)-.5)*.12;n/=math.sqrt(np.dot(n,n))
   alb[j,i]=np.minimum(255,np.maximum(0,np.sqrt(col)*255)).astype(np.uint8)
   norm[j,i]=np.minimum(255,np.maximum(0,(n*.5+.5)*255)).astype(np.uint8)
   enc=int(min(65535.,math.log2(1+t)/16.*65535))
   info[j,i,0]=enc//256;info[j,i,1]=enc%256;info[j,i,2]=kind
   # Save the horizon as a shadow aid; cap under-canopy self-occlusion for leaves.
   hw=horizon(p+np.array([0.,.3,0.]),T,1);he=horizon(p+np.array([0.,.3,0.]),T,-1)
   if kind==3 or kind==4:hw*=.4;he*=.4
   shad[j,i,0]=int(min(255.,max(0.,hw/(math.pi/2)*255)))
   shad[j,i,1]=int(min(255.,max(0.,he/(math.pi/2)*255)))
   shad[j,i,2]=int(ao*255)
 return alb,norm,info,shad

if __name__=='__main__':
 import argparse
 p=argparse.ArgumentParser();p.add_argument('--width',type=int,default=3072);a=p.parse_args()
 dest=ROOT/'data'/'landscape';dest.mkdir(exist_ok=True)
 P,T=scene();B,O=bvh(P);print('Primitives',len(P),'trees',len(T),'BVH',len(B),flush=True)
 t=time.time();imgs=render(a.width,a.width//2,P,B,O,T)
 for name,img in zip(('albedo','normal','geometry','shadow'),imgs):Image.fromarray(img).save(dest/(name+'.png'))
 (dest/'scene.json').write_text(json.dumps(dict(width=a.width,height=a.width//2,primitives=len(P),trees=len(T),camera=[0,float(terrain(0.,0.)+1.7),0],water_level_m=-3,scene='Illustrative wooded cove; local tangent geometry',seconds=time.time()-t),indent=2))
 print('Done',round(time.time()-t,1),flush=True)
