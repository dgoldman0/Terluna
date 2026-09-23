"""Compile and link the exact viewer shaders in Mesa OpenGL ES 3 via EGL.
This is a shader check, not a claim of browser WebGL execution.
"""
import ctypes as C,json,re
from pathlib import Path
E=C.CDLL('libEGL.so.1')
def egl(name,restype,args):
 f=getattr(E,name);f.restype=restype;f.argtypes=args;return f
ptr=C.c_void_p;integer=C.c_int;uint=C.c_uint
getproc=egl('eglGetProcAddress',ptr,[C.c_char_p])
getplatform=C.CFUNCTYPE(ptr,uint,ptr,ptr)(getproc(b'eglGetPlatformDisplayEXT'))
display=getplatform(0x31DD,None,None)
major=integer();minor=integer();assert egl('eglInitialize',uint,[ptr,C.POINTER(integer),C.POINTER(integer)])(display,C.byref(major),C.byref(minor))
assert egl('eglBindAPI',uint,[uint])(0x30A0)
config=ptr();number=integer();attrs=(integer*11)(0x3033,1,0x3040,0x40,0x3024,8,0x3023,8,0x3022,8,0x3038)
assert egl('eglChooseConfig',uint,[ptr,C.POINTER(integer),C.POINTER(ptr),integer,C.POINTER(integer)])(display,attrs,C.byref(config),1,C.byref(number)) and number.value
context=egl('eglCreateContext',ptr,[ptr,ptr,ptr,C.POINTER(integer)])(display,config,None,(integer*3)(0x3098,3,0x3038))
surface=egl('eglCreatePbufferSurface',ptr,[ptr,ptr,C.POINTER(integer)])(display,config,(integer*5)(0x3057,1200,0x3056,750,0x3038))
assert context and surface
assert egl('eglMakeCurrent',uint,[ptr,ptr,ptr,ptr])(display,surface,surface,context)
def gl(name,restype,args):
 p=getproc(name.encode());assert p,name;return C.CFUNCTYPE(restype,*args)(p)
getstring=gl('glGetString',C.c_char_p,[uint]);record={'renderer':getstring(0x1F01).decode(),'version':getstring(0x1F02).decode()}
createshader=gl('glCreateShader',uint,[uint]);source=gl('glShaderSource',None,[uint,integer,C.POINTER(C.c_char_p),C.POINTER(integer)])
compile_=gl('glCompileShader',None,[uint]);getsh=gl('glGetShaderiv',None,[uint,uint,C.POINTER(integer)]);getlog=gl('glGetShaderInfoLog',None,[uint,integer,C.POINTER(integer),ptr])
shaders=[]
for id_,type_ in [('vert',0x8B31),('frag',0x8B30)]:
 text=(Path(__file__).parent/('scene.'+id_)).read_bytes();ss=C.c_char_p(text)
 s=createshader(type_);source(s,1,C.byref(ss),None);compile_(s);ok=integer();getsh(s,0x8B81,C.byref(ok));buf=C.create_string_buffer(16384);getlog(s,len(buf),None,buf);record[id_]={'compiled':bool(ok.value),'log':buf.value.decode()};shaders.append(s)
 if not ok.value:raise RuntimeError(record)
program=gl('glCreateProgram',uint,[])()
for s in shaders:gl('glAttachShader',None,[uint,uint])(program,s)
gl('glLinkProgram',None,[uint])(program);ok=integer();gl('glGetProgramiv',None,[uint,uint,C.POINTER(integer)])(program,0x8B82,C.byref(ok));buf=C.create_string_buffer(16384);gl('glGetProgramInfoLog',None,[uint,integer,C.POINTER(integer),ptr])(program,len(buf),None,buf);record['link']={'success':bool(ok.value),'log':buf.value.decode()}
print(json.dumps(record,indent=2));(Path(__file__).parent/'validation'/'shader_checks.json').write_text(json.dumps(record,indent=2));assert ok.value
