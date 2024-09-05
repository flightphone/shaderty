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

#define PI  3.14159265359
#define TAU 6.28318530718
#define rot(f) mat2(cos(f), -sin(f), sin(f), cos(f))

const float dist_infin = 7.0;
#define nn 128
float npp =15.;
float level = 0.995;
vec3 hsky = vec3(0., 2., 0);



vec3 point(vec3 p) {
    return floor(p*npp)/npp;
}


#define NOISE fbm
#define NUM_NOISE_OCTAVES 5

// Precision-adjusted variations of https://www.shadertoy.com/view/4djSRW
float hash(float p) { p = fract(p * 0.011); p *= p + 7.5; p *= p + p; return fract(p); }
//float hash(vec2 p) {vec3 p3 = fract(vec3(p.xyx) * 0.13); p3 += dot(p3, p3.yzx + 3.333); return fract((p3.x + p3.y) * p3.z); }
float hash(vec2 p) { return fract(1e4 * sin(17.0 * p.x + p.y * 0.1) * (0.1 + abs(sin(p.y * 13.0 + p.x)))); }  //!!!Best
//	<https://www.shadertoy.com/view/4dS3Wd>
//	By Morgan McGuire @morgan3d, http://graphicscodex.com


float noise(float x) {
	float i = floor(x);
	float f = fract(x);
	float u = f * f * (3.0 - 2.0 * f);
	return mix(hash(i), hash(i + 1.0), u);
}

float noise(vec2 x) {
	vec2 i = floor(x);
	vec2 f = fract(x);
	// Four corners in 2D of a tile
	float a = hash(i);
	float b = hash(i + vec2(1.0, 0.0));
	float c = hash(i + vec2(0.0, 1.0));
	float d = hash(i + vec2(1.0, 1.0));

	// Simple 2D lerp using smoothstep envelope between the values.
	// return mix(mix(a, b, smoothstep(0.0, 1.0, u.x)),
	//			mix(c, d, smoothstep(0.0, 1.0, u.x)),
	//			smoothstep(0.0, 1.0, u.y));

	// Same code, with the clamps in smoothstep and common subexpressions
	// optimized away.
	vec2 u = f * f * (3.0 - 2.0 * f);
	return mix(a, b, u.x) + (c - a) * u.y * (1.0 - u.x) + (d - b) * u.x * u.y;
}

// This one has non-ideal tiling properties that I'm still tuning
float noise(vec3 x) {
	const vec3 step = vec3(110, 241, 171);

	vec3 i = floor(x);
	vec3 f = fract(x);
 
	// For performance, compute the base input to a 1D hash from the integer part of the argument and the 
	// incremental change to the 1D based on the 3D -> 1D wrapping
    float n = dot(i, step);

	vec3 u = f * f * (3.0 - 2.0 * f);
	return mix(mix(mix( hash(n + dot(step, vec3(0, 0, 0))), hash(n + dot(step, vec3(1, 0, 0))), u.x),
                   mix( hash(n + dot(step, vec3(0, 1, 0))), hash(n + dot(step, vec3(1, 1, 0))), u.x), u.y),
               mix(mix( hash(n + dot(step, vec3(0, 0, 1))), hash(n + dot(step, vec3(1, 0, 1))), u.x),
                   mix( hash(n + dot(step, vec3(0, 1, 1))), hash(n + dot(step, vec3(1, 1, 1))), u.x), u.y), u.z);
}


float fbm(float x) {
	float v = 0.0;
	float a = 0.5;
	float shift = float(100);
	for (int i = 0; i < NUM_NOISE_OCTAVES; ++i) {
		v += a * noise(x);
		x = x * 2.0 + shift;
		a *= 0.5;
	}
	return v;
}


float fbm(vec2 x) {
	float v = 0.0;
	float a = 0.5;
	vec2 shift = vec2(100);
	// Rotate to reduce axial bias
    mat2 rot = mat2(cos(0.5), sin(0.5), -sin(0.5), cos(0.50));
	for (int i = 0; i < NUM_NOISE_OCTAVES; ++i) {
		v += a * noise(x);
		x = rot * x * 2. + shift;
		a *= 0.5;
	}
	return v;
}


float fbmA(vec2 x) {
	float v = 0.0;
	float a = 0.5;
	vec2 shift = vec2(100);
	// Rotate to reduce axial bias
    mat2 rot = mat2(cos(0.5), sin(0.5), -sin(0.5), cos(0.50));
	for (int i = 0; i < NUM_NOISE_OCTAVES; ++i) {
        float n = abs((0.5 - noise(x))*2.);
        n = 0.5 - n;
        n = n*n;
		v += a * n;
		x = rot * x * 2. + shift;
		a *= 0.5;
	}
	return v;
}

float fbm(vec3 x) {
	float v = 0.0;
	float a = 0.5;
	vec3 shift = vec3(100);
	for (int i = 0; i < NUM_NOISE_OCTAVES; ++i) {
		v += a * noise(x);
		x = x * 2.0 + shift;
		a *= 0.5;
	}
	return v;
}
vec3 sky = vec3(0.08, 0.42, 0.87);
vec3 line = vec3(1.0);

vec3 wood(vec2 p)
{
	float t = noise(vec2(p.x*2., p.y*6.));
    t = fract(t*10.);
    return vec3(t);
}
vec3 zebra(vec2 p)
{
	//float t = fbm(vec2(p.x*4., p.y*8.));
	//t = smoothstep(0., 0.01, t-0.5);

	//float t = fbmA(vec2(p.x*4., p.y*8.));
	//t = smoothstep(0., 0.01, t-0.1);
	
	float t = noise(vec2(p.x*0.6, p.y*3.));
	t = fract(t*20.);
	t = smoothstep(0., 0.01, t-0.5);
    //t = smoothstep(0., 0.01, t - 0.3)*smoothstep(0., 0.01, 0.5 - t);
	
	
	return vec3(1.-t);
   
}




vec3 bark(vec2 p)
{
	float t = fbm(vec2(p.x*2., p.y*2.));
	//float t = fbmA(vec2(p.x*2., p.y*2.)); 
	t = fract(t*10.);
    return vec3(t);
}

void mainImage(out vec4 fragColor, in vec2 fragCoord) {
    //vec2 p = vec2(fragCoord.x/iResolution.x, fragCoord.y/iResolution.y); //(-iResolution.xy + 2.0 * fragCoord) / iResolution.y;
    vec2 p = (-iResolution.xy + 2.0 * fragCoord) / iResolution.y;
    
    /*
    float t = noise(p.x * 4.);
    float l = smoothstep(0.01, 0.0, abs(t-p.y));
    vec3 col = mix(sky, line, l);
    */
    vec3 col = wood(p);
	//vec3 col = bark(p);
	//vec3 col = zebra(p);
	


    
    fragColor = vec4(col, 1.0);
}
/////=====================================================================================
void main() {
    vec4 fragColor = vec4(0);
    mainImage(fragColor, gl_FragCoord.xy);
    gl_FragColor = fragColor;
}