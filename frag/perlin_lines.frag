#ifdef GL_ES
precision mediump float;
#endif

uniform vec2 u_resolution;
uniform vec2 u_mouse;
uniform float u_time;
uniform sampler2D u_tex0;
uniform sampler2D u_tex1;

#define iResolution u_resolution
#define iTime u_time
#define iMouse u_mouse
#define iChannel0 u_tex0
#define iChannel1 u_tex1

#define texture texture2D

/////=====================================================================================
/*

*/
#define PI  3.14159265359
#define TAU 6.28318530718

float hash(vec2 p) {
	return fract(1e4 * sin(117.0 * p.x + p.y * 0.1) * (0.1 + abs(sin(p.y * 133.0 + p.x))));
} 

float noise(vec2 x) {
	vec2 i = floor(x);
	vec2 f = fract(x);
	// Four corners in 2D of a tile
	float a = hash(i);
	float b = hash(i + vec2(1.0, 0.0));
	float c = hash(i + vec2(0.0, 1.0));
	float d = hash(i + vec2(1.0, 1.0));

	vec2 u = f * f * (3.0 - 2.0 * f);
	return mix(a, b, u.x) + (c - a) * u.y * (1.0 - u.x) + (d - b) * u.x * u.y;
}


vec3 lines(vec2 p)
{
	float n2 = 5., n3 = 2., n = 4., t = 0.; 

    vec3 col0 = vec3(1., 0.5, 0.5), col1 = vec3(0.5, 1., 0.5), col2 = vec3(0.5, 0.5, 1.),
	col3 = vec3(1., 1., 0.5), col4 = vec3(1., 0.5, 1.), col5 = vec3(0.5, 1., 1.),
	colline = vec3(0.);
	mat2 rot = mat2(cos(0.5), sin(0.5), -sin(0.5), cos(0.50));
    p.x += iTime*0.2;
    
	vec2 x = p, shift = vec2(0.1, 0.2);
	
    x = rot*x + shift;
	t += 1.8*noise(x*n3);
	
    x = rot*x;
	t += 0.2*noise(x*n2);
    
	p.y += t;
	p.y *= n;
	
	vec3 col = col0;
	float ncol = mod(floor(p.y), 6.);
	if (ncol == 1.)
		col = col1;
	if (ncol == 2.)		
		col = col2;
	if (ncol == 3.)		
		col = col3;
	if (ncol == 4.)		
		col = col4;
	if (ncol == 5.)		
		col = col5;			
	
	float y = fract(p.y), h = 0.06, eps = 0.02, s1 = smoothstep(1. - h - eps, 1.-h, y),	
	s2 = smoothstep(h, h-eps, y);
	col = mix(col, colline, s1);
	col = mix(col, colline, s2);
	return col;	
}


void mainImage(out vec4 fragColor, in vec2 fragCoord) {
    //vec2 p = vec2(fragCoord.x/iResolution.x, fragCoord.y/iResolution.y); //(-iResolution.xy + 2.0 * fragCoord) / iResolution.y;
	vec2 p = (-iResolution.xy + 2.0 * fragCoord) / iResolution.y;
	vec3 col =  lines(p);
	fragColor = vec4(col, 1.0);
}
/////=====================================================================================
void main() {
	vec4 fragColor = vec4(0);
	mainImage(fragColor, gl_FragCoord.xy);
	gl_FragColor = fragColor;
}