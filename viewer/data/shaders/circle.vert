~include version;
~include camera;
uniform mat4 model=mat4(1.0f);
layout (location=0) in vec4 pos;
layout (location=1) in vec4 normal;
out vec3 vtx_normal;
out vec3 vtx_frg_pos;
void main()												
{
	gl_Position=pvm*model*vec4(pos.xyz,1.f);
	vtx_normal=vec3(pvm*model*vec4(normal.xyz,1.f));
	vtx_frg_pos=vec3(model*vec4(pos.xyz,1.f));
}	