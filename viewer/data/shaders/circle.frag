~include version;
~include material;
~include lights;
~include phong_dl_fast_func;
in vec3 vtx_normal;
in vec3 vtx_frg_pos;
out vec4 frag_color;
void main()
{
    vec3 norm=normalize(vtx_normal);

	float dif_coef=abs(dot(norm,vec3(1.f,1.f,1.f)));
    vec4 dif=dif_coef*vec4(.5f)*mat_dif+vec4(.3f);
	vec3 color=dif.rgb;

	if(gl_FrontFacing) frag_color=vec4(color,1.f);
	else frag_color=vec4(1.f-color.r,1.f-color.g,1.f-color.b,1.f);
}