"""Compile our GLSL with GLES 3 through Mesa EGL.

Three.js engine includes are replaced by a minimal declared syntax harness.
This checks our shader bodies and varying interfaces, not full Three.js linking,
material integration, WebGL portability, rendering correctness or image quality.
"""
import ctypes as C
import hashlib
import json
import re
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
E=C.CDLL('libEGL.so.1'); ptr=C.c_void_p; integer=C.c_int; uint=C.c_uint

def egl(name,restype,args):
    f=getattr(E,name);f.restype=restype;f.argtypes=args;return f
getproc=egl('eglGetProcAddress',ptr,[C.c_char_p])
getplatform=C.CFUNCTYPE(ptr,uint,ptr,ptr)(getproc(b'eglGetPlatformDisplayEXT'))
display=getplatform(0x31DD,None,None);major=integer();minor=integer()
assert egl('eglInitialize',uint,[ptr,C.POINTER(integer),C.POINTER(integer)])(display,C.byref(major),C.byref(minor))
assert egl('eglBindAPI',uint,[uint])(0x30A0)
config=ptr();number=integer();attrs=(integer*11)(0x3033,1,0x3040,0x40,0x3024,8,0x3023,8,0x3022,8,0x3038)
assert egl('eglChooseConfig',uint,[ptr,C.POINTER(integer),C.POINTER(ptr),integer,C.POINTER(integer)])(display,attrs,C.byref(config),1,C.byref(number)) and number.value
context=egl('eglCreateContext',ptr,[ptr,ptr,ptr,C.POINTER(integer)])(display,config,None,(integer*3)(0x3098,3,0x3038))
surface=egl('eglCreatePbufferSurface',ptr,[ptr,ptr,C.POINTER(integer)])(display,config,(integer*5)(0x3057,32,0x3056,32,0x3038))
assert context and surface
assert egl('eglMakeCurrent',uint,[ptr,ptr,ptr,ptr])(display,surface,surface,context)
def gl(name,restype,args):
    p=getproc(name.encode());assert p,name;return C.CFUNCTYPE(restype,*args)(p)
getstring=gl('glGetString',C.c_char_p,[uint]);record={'scope':__doc__,'renderer':getstring(0x1F01).decode(),'version':getstring(0x1F02).decode(),'programs':[]}
createshader=gl('glCreateShader',uint,[uint]);source=gl('glShaderSource',None,[uint,integer,C.POINTER(C.c_char_p),C.POINTER(integer)])
compile_=gl('glCompileShader',None,[uint]);getsh=gl('glGetShaderiv',None,[uint,uint,C.POINTER(integer)]);getlog=gl('glGetShaderInfoLog',None,[uint,integer,C.POINTER(integer),ptr])
common='#version 300 es\nprecision highp float;precision highp int;\n#define texture2D texture\nuniform vec3 cameraPosition;\n'
vertex='''#define attribute in
#define varying out
#define USE_INSTANCING
in mat4 instanceMatrix;
in vec3 position,normal;in vec2 uv;
uniform mat4 projectionMatrix,viewMatrix,modelMatrix,modelViewMatrix;uniform mat3 normalMatrix;
'''
fragment='''#define varying in
out vec4 outColor;
#define gl_FragColor outColor
'''
material='''struct DirectLight {vec3 color;};struct DirectionalLight {vec3 color;};
struct ReflectedLight {vec3 indirectDiffuse;vec3 indirectSpecular;};
uniform DirectionalLight directionalLight;
DirectLight directLight;ReflectedLight reflectedLight;
void getDirectionalLightInfo(DirectionalLight a,out DirectLight b){b.color=a.color;}
'''
includes={'tonemapping_fragment':'/* identity tone map for syntax check */',
          'colorspace_fragment':'/* identity colour space for syntax check */',
          'begin_vertex':'vec3 transformed=position;',
          'project_vertex':'gl_Position=projectionMatrix*modelViewMatrix*instanceMatrix*vec4(transformed,1.);',
          'roughnessmap_fragment':'float roughnessFactor=.8;',
          'lights_fragment_end':'reflectedLight.indirectDiffuse=vec3(1.);reflectedLight.indirectSpecular=vec3(1.);'}
def expand(text):
    def repl(m):
        if m[1] not in includes:raise ValueError('Unhandled engine include '+m[1])
        return includes[m[1]]
    return re.sub(r'#include\s+<([^>]+)>',repl,text)
for p in json.loads((ROOT/'validation/custom_shaders.json').read_text())['programs']:
    rec={'name':p['name']};sh=[]
    for stage,kind in [('vertex',0x8B31),('fragment',0x8B30)]:
        text=common+(vertex if stage=='vertex' else fragment)
        if stage=='fragment' and 'material-additions' in p['name']:text+=material
        text+=expand(p[stage]);ss=C.c_char_p(text.encode());s=createshader(kind)
        source(s,1,C.byref(ss),None);compile_(s);ok=integer();getsh(s,0x8B81,C.byref(ok));buf=C.create_string_buffer(16384);getlog(s,len(buf),None,buf)
        rec[stage]={'compiled':bool(ok.value),'log':buf.value.decode(),'checked_glsl_sha256':hashlib.sha256(text.encode()).hexdigest()};sh.append(s)
    program=gl('glCreateProgram',uint,[])()
    for s in sh:gl('glAttachShader',None,[uint,uint])(program,s)
    gl('glLinkProgram',None,[uint])(program);ok=integer();gl('glGetProgramiv',None,[uint,uint,C.POINTER(integer)])(program,0x8B82,C.byref(ok));buf=C.create_string_buffer(16384);gl('glGetProgramInfoLog',None,[uint,integer,C.POINTER(integer),ptr])(program,len(buf),None,buf)
    rec['link']={'success':bool(ok.value),'log':buf.value.decode()};record['programs'].append(rec)
    gl('glDeleteProgram',None,[uint])(program)
    for s in sh:gl('glDeleteShader',None,[uint])(s)
record['success']=all(p['link']['success'] and p['vertex']['compiled'] and p['fragment']['compiled'] for p in record['programs'])
(ROOT/'validation/shader-checks.json').write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps(record,indent=2))
assert record['success']
