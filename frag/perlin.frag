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
float npp = 15.;
float level = 0.995;
vec3 hsky = vec3(0., 2., 0);

vec3 point(vec3 p) {
	return floor(p * npp) / npp;
}

#define NOISE fbm
#define NUM_NOISE_OCTAVES 5

// Precision-adjusted variations of https://www.shadertoy.com/view/4djSRW
/*
float hash(float p) {
	p = fract(p * 0.011);
	p *= p + 7.5;
	p *= p + p;
	return fract(p);
}
*/
float hash(float n) { return fract(sin(n) * 437558.5453123); }
//float hash(vec2 p) {vec3 p3 = fract(vec3(p.xyx) * 0.13); p3 += dot(p3, p3.yzx + 3.333); return fract((p3.x + p3.y) * p3.z); }
float hash(vec2 p) {
	return fract(1e4 * sin(17.0 * p.x + p.y * 0.1) * (0.1 + abs(sin(p.y * 13.0 + p.x))));
}  //!!!Best
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
	return mix(mix(mix(hash(n + dot(step, vec3(0, 0, 0))), hash(n + dot(step, vec3(1, 0, 0))), u.x), mix(hash(n + dot(step, vec3(0, 1, 0))), hash(n + dot(step, vec3(1, 1, 0))), u.x), u.y), mix(mix(hash(n + dot(step, vec3(0, 0, 1))), hash(n + dot(step, vec3(1, 0, 1))), u.x), mix(hash(n + dot(step, vec3(0, 1, 1))), hash(n + dot(step, vec3(1, 1, 1))), u.x), u.y), u.z);
}

float fbm(float x) {
	float v = 0.0;
	float a = 0.5;
	float shift = float(100);
	for(int i = 0; i < NUM_NOISE_OCTAVES; ++i) {
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
	for(int i = 0; i < NUM_NOISE_OCTAVES; ++i) {
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
	for(int i = 0; i < NUM_NOISE_OCTAVES; ++i) {
		float n = abs((0.5 - noise(x)) * 2.);
		n = 0.5 - n;
		n = n * n;
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
	for(int i = 0; i < NUM_NOISE_OCTAVES; ++i) {
		v += a * noise(x);
		x = x * 2.0 + shift;
		a *= 0.5;
	}
	return v;
}

float CloudNoise(in vec2 xy) {
	float w = .7;
	float f = 0.0;

	for(int i = 0; i < 4; i++) {
		f += noise(xy) * w;
		w *= 0.45;
		xy *= 2.7;
	}
	return f;
}

vec3 sky = vec3(0.08, 0.42, 0.87);
vec3 pink = vec3(0.9, 0.6, 0.6);
vec3 line = vec3(1.0);
vec3 clouds(vec2 p) {
	p *= 4.;
	p.x += iTime;
	float f = (CloudNoise(p) - 0.55) * 5.;
	// Uses the ray's y component for horizon fade of fixed colour clouds...
	vec3 col = mix(sky, vec3(.75, .75, .72), clamp(f, 0.0, 1.0));
	return col;

}


vec3 ccol(vec3 nor)
{
	vec3 col = sky;
	vec3 light = normalize(vec3(0.0, 1.0, 2.5)); //light
	vec3 rd = normalize(vec3(1., 1., -1.));
	vec3 R = reflect(light, nor);
	float specular = pow(max(abs(dot(R, rd)), 0.), 16.);
	float difu = abs(dot(nor, light));
	col = col * (clamp(difu, 0., 1.0) + 0.5) + vec3(1., .7, .4) * specular;
	float fre = pow(clamp(dot(nor, rd) + 1., .0, 1.), 3.); // Fresnel, for some mild glow.
	col += vec3(.1, .1, 0.1) * fre; //?
	col = sqrt(col);
	return col;
}

vec3 wood(vec2 p) {
	float t = noise(vec2(p.x * 2., p.y * 6.));
	t = fract(t * 10.);
	return vec3(t);
}

vec3 zebra(vec2 p) {
	//float t = fbm(vec2(p.x*4., p.y*8.));
	//t = smoothstep(0., 0.01, t-0.5);

	//float t = fbmA(vec2(p.x*4., p.y*8.));
	//t = smoothstep(0., 0.01, t-0.1);

	float t = noise(vec2(p.x * 0.6, p.y * 3.));
	t = fract(t * 20.);
	t = smoothstep(0., 0.01, t - 0.5);
    //t = smoothstep(0., 0.01, t - 0.3)*smoothstep(0., 0.01, 0.5 - t);

	return vec3(1. - t);

}

vec3 bull(vec2 p) {
	
	float t = noise(vec2((p+3.)*4. ));
	//t =fract(t*5.);
	t = smoothstep(0., 0.01, t - 0.4);
    //t = smoothstep(0., 0.01, t - 0.3)*smoothstep(0., 0.01, 0.5 - t);

	return vec3(t);

}


vec3 tunnel(vec2 p) {
	p.x -= iTime*0.3;
	float t = noise(vec2(p*4. ));
	//t = smoothstep(0., 0.01, t - 0.4);
    t =fract(t*2.);
	t = smoothstep(0., 0.01, t - 0.3)*smoothstep(0., 0.01, 0.5 - t);
	return vec3(t);

}

vec3 bark(vec2 p) {
	float t = fbm(vec2(p.x * 2., p.y * 2.));
	//float t = fbmA(vec2(p.x*2., p.y*2.)); 
	t = fract(t * 10.);
	return vec3(t);
}

float map2(vec3 p) {
	
	float r = length(5.*p) - iTime*0.1;
	r = noise(r) +  .5 * noise(2.*r) +  .1* noise(r * 4.)+0.3*noise(r*8.);
	r = p.x*p.x*p.x + p.y*p.y*p.y - r;
	return r;
		
}

float map3(vec3 p) {
	float lon = atan(p.y, p.x);
	float lat = atan(p.z, length(p.xy));
	float r = sin(4. * lon) + cos(4. * lat) + iTime * 0.3;
	r =  noise(r) + 0.5 * noise(4.*r) +  .3* noise(r * 10.) + 0.1* noise(r*20.);
	r = length(p) - r;
	return r;
		
}

float map1(vec3 p) {
	//wave
	p += 45.;
	//p.x-= iTime*0.5;
	float r = noise(1.5*p.x) + noise(1.5*p.y);// + iTime * 0.1;
	r =  noise(r) + 0.5 * noise(4.*r) +  .3* noise(r * 10.) + 0.1* noise(r*20.);
	r = length(p) - r;
	r = r + 0.5*noise(r) + .2 * noise(2.*r);// +  0.1 * noise(r * 20.) + 0.1*noise(r*40.);
	return r;
	
}

float map(vec3 p) {
	//wave
	p += 15.;
	//p.x-= iTime*0.5;
	float r = noise(1.5*p.xy); //+ iTime * 0.2;
	r =  noise(r) + 0.5 * noise(4.*r) +  .3* noise(r * 10.) + 0.1* noise(r*20.);
	r = length(p) - r;
	r = r + 0.5*noise(r) + .2 * noise(2.*r);// +  0.1 * noise(r * 20.) + 0.1*noise(r*40.);
	return r;
	
}

float fbmW(vec2 x) {
	//x += iTime*0.2;
	float v = 0.0;
	float a = 0.5;
	vec2 shift = vec2(0.);//vec2(iTime*0.1);
	// Rotate to reduce axial bias
	mat2 rot = mat2(cos(0.5), sin(0.5), -sin(0.5), cos(0.50));
	x = rot*x + shift ;
	v += noise(x); 
	x = rot*x + shift;
	v += 1.5 * noise(2.*x);// + noise(iTime*0.3);
	return v;
	
}

float fbmW2(vec2 x) {
	//x += iTime*0.2;
	float v = 0.0;
	float a = 0.5;
	vec2 shift = vec2(0.);//vec2(iTime*0.1);
	// Rotate to reduce axial bias
	mat2 rot = mat2(cos(0.5), sin(0.5), -sin(0.5), cos(0.50));
	x = rot*x + shift;
	v += noise(x); 
	x = rot*x + shift;
	v += .5 * noise(4.*x);
	return v;
	
}


float patern(vec2 p, out vec2 q, out vec2 r ) {
	q = vec2( fbmW( p + vec2(0.0,0.0) ),
                   fbmW( p + vec2(5.2,1.3) ) );

    r = vec2( fbmW( p + 2.0*q + vec2(1.7,9.2) ),
                   fbmW( p + 1.0*q + vec2(8.3,2.8) ) );

    return /*p.x*p.x*p.x*p.x + p.y*p.y*p.y*p.y -*/ fbmW( p + 3.0*r );
	
}

vec3 warp(vec2 pos) 
{
	vec2 q = vec2(0), r = vec2(0);
	vec3 col = vec3(0.64, 0.76, 0.73);
	float res = patern(pos, q, r);
	//q = normalize(q);
	//r = normalize(r);
	
	col = mix(vec3(0.74, 0.65, 0.74), col, smoothstep(fract(length(q)), 0.3, 0.6));
	col = mix(vec3(0.86, .85, 0.82), col,  smoothstep(fract(length(r)), 0.3, 0.7));
	col = mix(vec3(0.5, 0.5, 0.5), col, smoothstep(fract(res), 0.3, 0.7));
	
	return col;
}


vec3 pearl(vec2 pos) 
{
	vec2 q = vec2(0), r = vec2(0);
	vec3 col = vec3(0.64, 0.76, 0.73);
	float res = patern(pos, q, r);
	
	
	col = mix(vec3(0.74, 0.65, 0.74), col, fract(noise(q*10.)));
	col = mix(vec3(0.86, .85, 0.82), col,  fract(noise(r*10.)));
	col = mix(vec3(0.1, 0.1, 0.1), col, fract(noise(res*1.)));
	
	return col;
}

vec3 radial(vec2 p) 
{
	//p.xy *= rot(iTime/3.);
	float q = noise(4.*p.xy);
	//q = length(p) - q;
	float r =  noise(q) + 0.5 * noise(4.*q) +  .3* noise(q * 10.) + 0.1* noise(q*20.);
	r = length(p) - r;//  - iTime * 0.2;
	//r = p.x*p.x /2. + p.y*p.y/4. - r;
	//r = sdBox(p, vec2(r));
	//r = sdEquilateralTriangle(p, r);
	float res = r + 0.5*noise(r) + .2 * noise(2.*r)- iTime * 0.2;;
	vec3 col = vec3(0.64, 0.76, 0.73);
	
	/*
	col = mix(vec3(0.74, 0.65, 0.74), col, q);
	col = mix(vec3(0.86, .85, 0.82), col,  r);
	col = mix(vec3(0.5, 0.5, 0.5), col, res);
	*/
	
	col = mix(vec3(0.74, 0.65, 0.74), col, fract(noise(q*5.)));
	col = mix(vec3(0.86, .85, 0.82), col,  fract(noise(r*5.)));
	col = mix(vec3(0.3, 0.3, 0.3), col, fract(noise(res*5.)));
	
	
	
	return col;
}
float fun(vec2 p)
{
	//float t =  fbmW(p);
	//return t*t;
	return length(wood(p)) - 0.3;
}

float map8(vec3 p)
{
	//return fbm(p);
	//return fun(p.xy) - p.z;
	//return length(pearl(p.xy));
	
	
	mat2 rot = mat2(cos(0.5), sin(0.5), -sin(0.5), cos(0.50));
	float t = noise(p); 
	
	p.xy *= rot;
	p*=2.2;
	t += .5*noise(p) ;
	return t;
	
	//float t = noise(p + 0.2*p*2. + 0.1*p*3.);
	//return t;
	
}



vec3 calcNormal8(in vec3 p) {
	const float eps = 0.0001;
	vec2 q = vec2(0.0, eps);
	vec3 res = vec3(map8(p + q.yxx) - map8(p - q.yxx), map8(p + q.xyx) - map8(p - q.xyx), map8(p + q.xxy) - map8(p - q.xxy));
	return normalize(res);
}
vec3 drops(vec2 p) {
	vec3 pos = vec3(p.xy*5., 3.5);
	vec3 nor = calcNormal8(pos);
	return ccol(nor);
	//float nor = map8(pos);
	//return vec3(fract(nor));
	
}

vec3 calcNormal(in vec3 p) {
	const float eps = 0.0001;
	vec2 q = vec2(0.0, eps);
	vec3 res = vec3(map(p + q.yxx) - map(p - q.yxx), map(p + q.xyx) - map(p - q.xyx), map(p + q.xxy) - map(p - q.xxy));
	return normalize(res);
}





vec3 wave(vec2 p) {
	vec3 pos = vec3(p, 3.);
	vec3 nor = calcNormal(pos);
	return ccol(nor);
}


void mainImage(out vec4 fragColor, in vec2 fragCoord) {
    //vec2 p = vec2(fragCoord.x/iResolution.x, fragCoord.y/iResolution.y); //(-iResolution.xy + 2.0 * fragCoord) / iResolution.y;
	vec2 p = (-iResolution.xy + 2.0 * fragCoord) / iResolution.y;

    /*
    float t = noise(p.x * 4.);
    float l = smoothstep(0.01, 0.0, abs(t-p.y));
    vec3 col = mix(sky, line, l);
    */
	
    //vec3 col = wood(p);
	//vec3 col = bark(p);
	//vec3 col = zebra(p);
	//vec3 col = clouds(p);
	//vec3 col = drops(p);
	vec3 col = pearl(p*0.5 +15.);
	//vec3 col = wave(p);
	//vec3 col = warp(p*0.5);
	//vec3 col = radial(p*4.);
	//vec3 col = bull(p);
	//vec3 col = tunnel(p);
	

	fragColor = vec4(col, 1.0);
}
/////=====================================================================================
void main() {
	vec4 fragColor = vec4(0);
	mainImage(fragColor, gl_FragCoord.xy);
	gl_FragColor = fragColor;
}