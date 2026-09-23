#version 300 es
in vec2 a;out vec2 uv;void main(){uv=.5*a+.5;gl_Position=vec4(a,0.,1.);}
