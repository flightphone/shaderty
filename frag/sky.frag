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


// Precision-adjusted variations of https://www.shadertoy.com/view/4djSRW
//float hash(float p) { p = fract(p * 0.011); p *= p + 7.5; p *= p + p; return fract(p); }
//float hash(vec2 p) {vec3 p3 = fract(vec3(p.xyx) * 0.13); p3 += dot(p3, p3.yzx + 3.333); return fract((p3.x + p3.y) * p3.z); }
float hash (vec3 p) {
    return fract(sin(dot(p, vec3(127.1,311.7, 74.7))) * 43758.5453123);
}

//	<https://www.shadertoy.com/view/4dS3Wd>
//	By Morgan McGuire @morgan3d, http://graphicscodex.com
//
float hash(float n) { return fract(sin(n) * 43758.5453123); }
//float hash(vec2 p) { return fract(1e4 * sin(17.0 * p.x + p.y * 0.1) * (0.1 + abs(sin(p.y * 13.0 + p.x)))); }
float hash (vec2 p) { return fract(sin(dot(p, vec2(1227.1,3311.7))) * 34558.5453123);}
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
	 return mix(mix(a, b, smoothstep(0.0, 1.0, f.x)),
				mix(c, d, smoothstep(0.0, 1.0, f.x)),
				smoothstep(0.0, 1.0, f.y));

	// Same code, with the clamps in smoothstep and common subexpressions
	// optimized away.
	//vec2 u = f * f * (3.0 - 2.0 * f);
	//return mix(a, b, u.x) + (c - a) * u.y * (1.0 - u.x) + (d - b) * u.x * u.y;
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

float FractalNoise(in vec3 xy)
{
	float w = .7;
	float f = 0.0;

	for (int i = 0; i < 4; i++)
	{
		f += noise(xy) * w;
		w *= 0.45;
		xy *= 2.7;
	}
	return f;
}

float fbm(vec2 n) {
    mat2 m = mat2( 1.6,  1.2, -1.2,  1.6 );
	float total = 0.0, amplitude = 0.1;
	for (int i = 0; i < 7; i++) {
		total += noise(n) * amplitude;
		n = m * n;
		amplitude *= 0.4;
	}
	return total;
}

vec2 lonlat (vec3 p)
{
    
    float lon = atan(p.y, p.x)/TAU;
    float lat = atan(length(p.xy), p.z)/PI;
    return vec2(1.0-lon, lat);
}

vec3 Clouds(vec3 sky, vec3 pos)
{
    float f = (FractalNoise(pos) - 0.55)*5.;
	// Uses the ray's y component for horizon fade of fixed colour clouds...
	sky = mix(sky, vec3(.75, .75, .72), clamp(f, 0.0, 1.0));
    return sky;
}
vec3 GetClouds(vec3 sky, vec3 ro, vec3 rd)
{
    
    float t = (dot(hsky,hsky) - dot(ro,hsky))/dot(rd, hsky);
    if (t > 0. && t < 100. && rd.y > 0.05)
    {
        vec3 pos = ro + t*rd;
        float f = (FractalNoise(pos) - 0.55)*5.;
        //float f = fbm(pos.xz)*5.;
        
	    // Uses the ray's y component for horizon fade of fixed colour clouds...
	    sky = mix(sky, vec3(.75, .75, .72), clamp(f, 0.0, 1.0));
    }
    return sky;
}
/*
vec3 GetClouds(in vec3 sky, in vec3 rd, vec3 cameraPos)
{
	
    if (rd.y < 0.01) return sky;
	float v = (200.0-cameraPos.y)/rd.y;
	rd.xz *= v;
	rd.xz += cameraPos.xz;
	rd.xz *= .010;
	float f = (FractalNoise(rd.xz) -.55) * 5.0;
	// Uses the ray's y component for horizon fade of fixed colour clouds...
	//sky = mix(sky, vec3(.55, .55, .52), clamp(f*rd.y-.1, 0.0, 1.0));
    sky = mix(sky, vec3(.55, .55, .52), clamp(f-0.1, 0.0, 1.0));
    return sky;
}
*/

vec3 GetRayDir(vec2 uv, vec3 p, vec3 l, float z) {
    vec3 f = normalize(l - p), r = normalize(vec3(f.z, 0, -f.x)), u = cross(f, r), c = f * z, i = c + uv.x * r + uv.y * u;
    return normalize(i);
}
/*
#if HW_PERFORMANCE==0
#define AA 1
#else
#define AA 2
#endif
*/
#define AA 1

void mainImage(out vec4 fragColor, in vec2 fragCoord) {
    vec3 light = normalize(vec3(0.0, .0, 1.)); //light
    vec3 light2 = normalize(vec3(0.0, 0.0, -1.)); //light
    vec2 mo = 1.5 * cos(0.5 * iTime + vec2(0, 11));
    //if  (iMouse.z > 0.0)
    {
        mo = (-iResolution.xy + 2.0 * (iMouse.xy)) / iResolution.y;
    }
    vec3 ro = vec3(0.0, 0.0, 2.5); // camera
    //camera rotation
    ro.yz *= rot(mo.y);
    ro.xz *= rot(-mo.x - 1.57);

    const float fl = 1.5; // focal length
    

    vec3 sky = vec3(0.08, 0.42, 0.87);
    //vec3 b1 = vec3(0.0509, 0.2980, 0.4705), b2 = vec3(0.3764, 0.7529, 0.8784), bg = mix(b1, b2, vec3((1.0 - abs(fragCoord.x - iResolution.x / 2.) / iResolution.y * 2.) * fragCoord.y / iResolution.x));   
    vec3 bg = sky;
    //antialiasing
    vec3 tot = vec3(0.0);
    for(int m = 0; m < AA; m++) for(int n = 0; n < AA; n++) {
            vec2 o = vec2(float(m), float(n)) / float(AA) - 0.5;
            vec2 p = (-iResolution.xy + 2.0 * (fragCoord + o)) / iResolution.y;
            vec3 rd = GetRayDir(p, ro, vec3(0, 0., 0), fl); //ray direction
            vec3 col = GetClouds(sky, ro, rd);
            //vec3 col = bg; // background  
            //sky
            float d = length(cross(ro, rd));
            float td = abs(dot(ro, rd));
            d = sqrt(dist_infin * dist_infin - d * d);
            vec3 pos = ro + rd * (td + d);
            //col = Clouds(sky, pos);
            
            
            //vec3 pos = vec3(lonlat(rd), 0.);
            vec3 pp = point(pos);
            float fil = hash(pp);
            if(fil > level) {
                //cyrcle
                if ((length(pos - (pp + vec3(0.5/npp, 0.5/npp, 0.5/npp)))) < 0.5/npp)
                    col = vec3(1.);
            }
            
            tot += col;
        }
    tot = tot / float(AA) / float(AA);
    //tot = tot / float(AA);
    //antialiasing
    fragColor = vec4(tot, 1.0);
}
/////=====================================================================================
void main() {
    vec4 fragColor = vec4(0);
    mainImage(fragColor, gl_FragCoord.xy);
    gl_FragColor = fragColor;
}