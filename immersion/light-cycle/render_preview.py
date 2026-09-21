"""Render the exact GPU shader offscreen with Mesa EGL. This is not browser validation."""
import ctypes as C,json,math,argparse
from pathlib import Path
import numpy as np
from PIL import Image
import egl_scene as E
ROOT=Path(__file__).parent;gl=E.gl;u=C.c_uint;i=C.c_int;f=C.c_float;p=C.c_void_p
P=E.program;gl('glUseProgram',None,[u])(P)
vert=np.array([-1,-1,1,-1,-1,1,-1,1,1,-1,1,1],np.float32);buf=u();gl('glGenBuffers',None,[i,C.POINTER(u)])(1,C.byref(buf));gl('glBindBuffer',None,[u,u])(0x8892,buf);gl('glBufferData',None,[u,C.c_ssize_t,p,u])(0x8892,vert.nbytes,vert.ctypes.data,0x88E4);loc=gl('glGetAttribLocation',i,[u,C.c_char_p])(P,b'a');gl('glEnableVertexAttribArray',None,[u])(loc);gl('glVertexAttribPointer',None,[u,i,u,u,i,p])(loc,2,0x1406,0,0,None)
uniformloc=gl('glGetUniformLocation',i,[u,C.c_char_p]);uf=gl('glUniform1f',None,[i,f]);ui=gl('glUniform1i',None,[i,i])
def texture(unit,target,arr,internal,format_,type_,linear=True):
 t=u();gl('glGenTextures',None,[i,C.POINTER(u)])(1,C.byref(t));gl('glActiveTexture',None,[u])(0x84C0+unit);gl('glBindTexture',None,[u,u])(target,t)
 for key in [0x2801,0x2800]:gl('glTexParameteri',None,[u,u,i])(target,key,0x2601 if linear else 0x2600)
 for key in [0x2802,0x2803]:gl('glTexParameteri',None,[u,u,i])(target,key,0x812F)
 gl('glPixelStorei',None,[u,i])(0x0CF5,1)
 if target==0x8C1A:
  gl('glTexImage3D',None,[u,i,i,i,i,i,i,u,u,p])(target,0,internal,arr.shape[2],arr.shape[1],arr.shape[0],0,format_,type_,arr.ctypes.data)
 else:gl('glTexImage2D',None,[u,i,i,i,i,i,u,u,p])(target,0,internal,arr.shape[1],arr.shape[0],0,format_,type_,arr.ctypes.data)
 return t

def run(world='moon',angle=90,weather='scattered',yaw=0,pitch=.04,ev=0,out='preview_noon.png',width=1200,height=750):
 d=np.load(ROOT/'data'/f'{world}_atlas.npz');meta=json.loads(str(d['metadata']));atm=meta['atmosphere'];rgb=d['rgb'];layers,nh,nw,_=rgb.shape;scale=np.maximum(np.max(np.abs(rgb),axis=(1,2,3))/1000,1e-20)
 sky=np.ones((layers,nh,nw,4),np.float16);sky[...,:3]=(rgb/scale[:,None,None,None]).astype(np.float16)
 tx=texture(0,0x8C1A,sky,0x881A,0x1908,0x140B);ui(uniformloc(P,b'sky'),0)
 md=np.zeros((layers,15,4),np.float32);md[:,0,0]=scale
 for col,name in [(1,'direct'),(2,'direct_horizontal'),(3,'diffuse'),(13,'cloud_direct'),(14,'cloud_diffuse')]:md[:,col,:3]=d[name]
 md[:,4:13,:3]=d['sh'];tm=texture(1,0x0DE1,md,0x8814,0x1908,0x1406,False);ui(uniformloc(P,b'metadata'),1)
 scene=json.loads((ROOT/'data'/'landscape'/'scene.json').read_text())
 for j,name in enumerate(['albedo','normal','geometry','shadow']):
  arr=np.array(Image.open(ROOT/'data'/'landscape'/f'{name}.png').convert('RGB'));texture(j+2,0x0DE1,arr,0x8051,0x1907,0x1401,name!='geometry');ui(uniformloc(P,(name+'Map').encode()),j+2)
 presets={'clear':(0,0,0,1500,0,0),'scattered':(.36,6,0,1500,0,0),'overcast':(1,24,3.912/35000,1500,0,0),'haze':(.06,3,3.912/8000,850,0,0),'fog':(.78,10,3.912/220,65,0,.15),'rain':(1,34,3.912/2400,450,.8,.85)}
 cover,tau,fog,fh,rain,wet=presets[weather];angles=d['suns'];idx=np.clip(np.searchsorted(angles,angle)-1,0,len(angles)-2);blend=np.clip((angle-angles[idx])/(angles[idx+1]-angles[idx]),0,1)
 uniforms=dict(layer0=idx,layer1=idx+1,blend=blend,exposureValue=(1/8500)*2**ev,sunElevation=math.radians(angle),aspect=width/height,fov=math.radians(68),yaw=yaw,pitch=pitch,radius=atm['radius_m'],sunSide=1,weatherTime=84,cloudCover=cover,cloudTau=tau,fogK=fog,fogHeight=fh,rainAmount=rain,wetness=wet,cameraY=scene['camera'][1],quality=1)
 for name,val in uniforms.items():uf(uniformloc(P,name.encode()),float(val))
 gl('glUniform2f',None,[i,f,f])(uniformloc(P,b'texSize'),nw,nh)
 beta=1.24062e-6*(np.array([680.,550.,440.])/1000)**-4*atm['density_scale'];gl('glUniform3f',None,[i,f,f,f])(uniformloc(P,b'rayBeta'),*beta)
 gl('glViewport',None,[i,i,i,i])(0,0,width,height);gl('glDrawArrays',None,[u,i,i])(0x0004,0,6);gl('glFinish',None,[])()
 im=np.empty((height,width,4),np.uint8);gl('glReadPixels',None,[i,i,i,i,u,u,p])(0,0,width,height,0x1908,0x1401,im.ctypes.data);err=gl('glGetError',u,[])();assert err==0,hex(err)
 dest=ROOT/'validation'/out;Image.fromarray(im[::-1]).save(dest);print(dest,flush=True)
 return dest
if __name__=='__main__':
 a=argparse.ArgumentParser();a.add_argument('--world',default='moon');a.add_argument('--angle',type=float,default=90);a.add_argument('--weather',default='scattered');a.add_argument('--yaw',type=float,default=0);a.add_argument('--pitch',type=float,default=.04);a.add_argument('--ev',type=float,default=0);a.add_argument('--out',default='preview_noon.png');x=a.parse_args();run(**vars(x))
