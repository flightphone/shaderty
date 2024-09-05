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

//noise

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




/////=====================================================================================
//Collection of implicit surfaces. implicit surfaces, raytracing, binary search
/*
Rendering implicit surfaces. Using raytracing and binary searchy. 
Here, these same surfaces are obtained by creating grids using an algorithm 
3D Marching Cubes: https://flightphone.github.io/paramgeometry.html
*/
#define PI  3.14159265359
#define TAU 6.28318530718
#define rot(f) mat2(cos(f), -sin(f), sin(f), cos(f))
#define nn 128.
#define newton 5
float csurf = 0.;
float scale = 10.;



float glz() {

    float t = iTime / 5.;
    float st = mod(floor(t), 4.);
    float res;
    if(st == 0.)
        res = 1.;
    if(st == 1.)
        res = cos(fract(t) * PI / 2.);//(1.- fract(t))*(1.- fract(t));
    if(st == 2.)
        res = 0.;
    if(st == 3.)
        res = sin(fract(t) * PI / 2.); //fract(t)*fract(t);   
    return res;
}

float rottime()
{
    float t = iTime / 5.;
    
    float st = mod(floor(t), 4.);
    float ct = floor(t/4.)*4.;
    float res = ct/2.0;
    if(st == 0.)
        res += t-ct;
    if(st == 1.)
        res += 1.;    
    if(st == 2.)
        res += t-ct -1.;
    if(st == 3.)
        res += 2.;
    return res*5.;
    
}

float isf(vec3 p) {
    float x = p.x, y = -p.z, z = p.y;
    return (2. * y * (y * y - 3. * x * x) * (1. - z * z) + (x * x + y * y) * (x * x + y * y) - (9. * z * z - 1.) * (1. - z * z));// IMPLICIT SURFACE Function
}



float capsule(vec3 p)
{
    float h = 1.5;
    float x = p.x, y = p.y, z = 0.;
    if (p.z > h)
        z = p.z - h;
    if (p.z < -h)
        z = p.z + h; 

    return x*x + y*y + z*z - 0.4*0.4;       

}

float capsule2(vec3 p)
{
    float h = 1.5;
    float x = p.x, y = p.y, z = 0.;
    if (p.z > h)
        z = p.z - h;
    if (p.z < -h)
        z = p.z + h; 

    float a = length(p - vec3(0., 1., 0.));
    return 1./(x*x + y*y + z*z) - 1./0.2/0.2 + 20./a;       

}
float hyper(vec3 p)
{
    float x = p.x, y = p.y, z = p.z;
    return pow(x, 2./3.) + pow(y, 2./3.) + pow(z, 2./3.) - 1.;
}
float sine(vec3 p)
{
    float x = p.x, y = p.y, z = p.z, a = 1.3;
    return 4.*x*x*y*y*z*z + a*a*(x-y-z)*(x+y-z)*(x-y+z)*(x+y+z)-0.005;
}

float tooth(vec3 p)
{
     float x = p.x*p.x, y = p.y*p.y, z = p.z*p.z; 
     return x*x + y*y + z*z - (x+y+z);
}
float hant(vec3 p)
{
     float x = p.x, y = p.y, z = p.z; 
     float a = x*x+y*y+z*z-13.;
     float b = 3.*x*x + y*y - 4.*z*z - 12.;
     return 4.*a*a*a + 27.*b*b;
}
vec2 lonlat (vec3 p)
{
    
    float lon = atan(p.y, p.x)/2.0;
    float lat = atan(length(p.xy), p.z);
    return vec2(1.0-lon, lat);
}
float sphere(vec3 p) {
    float x = p.x, y = p.y, z = p.z;
    float r = 1.5;
    vec2 ll = lonlat(p);
    return (x * x + y * y + z * z - r*r + 0.4*(sin(20.*ll.x)+sin(20.*ll.y)));
}

float plane(vec3 p)
{
    return dot(p, vec3(0., 0., 1.));
}


float gyroide(vec3 p) {
    float x = p.x*scale, y = p.y*scale, z = p.z*scale; 
    return (cos(x) * sin(y) + cos(y) * sin(z) + cos(z) * sin(x));
}

float eggbox(vec3 p)
{
    return -0.1*(sin(p.x*scale) + sin(p.y*scale)) + p.z;
}

float eggbox2(vec3 p)
{
    return eggbox(p)*(eggbox(p)) - 0.01;
}

float gyroide2(vec3 p)
{
    scale = 6.;
    return gyroide(p)*(gyroide(p)) - 0.1;
}

float sphere3(vec3 p)
{
     return min(sphere(p),capsule(p));
    //return sphere(p);
     //return sphere(p) * capsule(p) - 0.01;
}

float sphere2(vec3 p)
{
    return length(p-vec3(0., 0., .5)) - 0.6;
}



vec3 knot(float t)
{
    return vec3(sin(t) + 2.*sin(2.*t)-0.1, 
    cos(t) - 2.*cos(2.*t), sin(3.*t)-0.1);
}

vec3 curve(float z, float s)
{
    float h = 1.;
    z = clamp(z, -h, h);
    float t = asin(z)/3. + s;
    return knot(t);
    
    //
}

vec3 curve2(float z, float s)
{
    float h = 1.;
    z = clamp(z, -h, h);
    float t = (PI - asin(z))/3. + s;
    return knot(t);
    //
}

vec3 curve0(float z)
{
    float h = 2.1;
    z = clamp(z, -h, h);
    return vec3(0.5*cos(z*PI/2.), 0.5*sin(z*PI/2.), z);
    //
}
float sline0(vec3 p)
{
    float res = 1.;
    vec3 a = curve0(p.z);
    return  (length(p-a) - 0.3);
}
float sline(vec3 p)
{
    
    float res = 1., d = 0.5;
    vec3 a = curve(p.z, 0.);
    res *=  (d*length(p-a) - 0.2);
    a = curve(p.z, TAU/3.);
    res *= (d*length(p-a) - 0.2);
    a = curve(p.z, 2.*TAU/3.);
    res *= (d*length(p-a) - 0.2);
    a = curve2(p.z, 0.);
    res *= (d*length(p-a) - 0.2);
    a = curve2(p.z, TAU/3.);
    res *= (d*length(p-a) - 0.2);
    a = curve2(p.z, 2.0*TAU/3.);
    res *= (d*length(p-a) - 0.2);
    
    
    return res;
}

float floo(vec3 p)
{
    return sphere(p) * capsule(p) - 0.01;
}

float zin(vec3 p)
{
    float x = p.x, y = p.y, z = p.z, a = 1.;
    return z*(x*x - y*y) - 2.*a*x*y;
}

float digdong(vec3 p)
{
    float x = p.x, y = p.y, z = p.z, a = 2.;
    //z = clamp(z, 0., 10.);
    return (x*x + y*y)*2.5 - (a-z)*z*z;
}

float umb(vec3 p)
{
     float x = p.x, y = p.y, z = p.z, a = 2.;
    return x*x - z*y*y;
}

float umb2(vec3 p)
{
     float x = p.x, y = p.y, z = p.z, a = 2.;
     float r =  x*x - z*y*y;
    
    if (z < -0.01)
        r += 1.;
    return r - 0.01;    
}

float eight (vec3 p)
{
    float k = 0.6, x = p.x*k, y = p.y*k, z = p.z*k, a = 1.;
    //z = clamp(z, 0., 10.);
    return 4.*z*z*z*z + a*a*(x*x + y*y - 4.*z*z) ;
}

float eight2 (vec3 p)
{
    float x = p.x, y = p.y, z = p.z, a = 1.1;
    return 4.*z*z*z*z + a*a*(x*x + y*y - 4.*z*z) + sin(4.*x) + sin(4.*y) + sin(4.*z);
}

float cayley(vec3 p)
{
    
    float k = 0.7, x = p.x*k, y = p.y*k, z = p.z*k, a = 0.9;
    if (length(p) > 2.)
        return 1.;
    
    return x*x + y*y - x*x*z + y*y*z + z*z - a;
}

float chair(vec3 p)
{
    
    float s = 0.6, x = p.x*s, y = p.y*s, z = p.z*s, a = 0.8, k = 1.2, b = 0.4, 
    d = x*x + y*y + z*z - a*k*k;
       
    return  d*d - b*((z - k)*(z-k) - 2.*x*x)*((z+k)*(z+k) - 2.*y*y);
}

float miter(vec3 p)
{
    float x = p.x, y = p.y, z = p.z;
    return  4.*x*x*(x*x + y*y + z*z) - y*y*(2.5 - y*y - z*z);
}

float piri(vec3 p)
{
    float x = p.x, y = p.y, z = p.z, a = 1.6;
    return  (x*x*x*x - a*x*x*x) + a*a*(y*y+z*z);
}

float roman(vec3 p)
{
    float x = p.x, y = p.y, z = p.z, a = 1.6;
    return  x*x*y*y + y*y*z*z + z*z*x*x - a*a*x*y*z;
}

float piri2(vec3 p)
{
    float x = p.x, y = p.y, z = p.z, a = 1.;
    return  (x*x*x*x - a*x*x*x) + a*a*(y*y+z*z)+ sin(4.*x) + sin(4.*y) + sin(4.*z);;
}

float barth(vec3 p)
{
    //https://mathworld.wolfram.com/BarthSextic.html
    float x = p.x, y = p.y, z = p.z, f = 1., w = 1.;
    return -4.*(f*f*x*x - y*y)*(f*f*y*y - z*z)*(f*f*z*z - x*x) + (1. + 2.*f)*(x*x + y*y + z*z - w*w)*(x*x + y*y + z*z - w*w)*w*w;
}

float barth2(vec3 p)
{

    float k = 0.8, x = p.x*k, y = p.y*k, z = p.z*k, f = 1.1, w = 1.;
    return -4.*(f*f*x*x - y*y)*(f*f*y*y - z*z)*(f*f*z*z - x*x) + (1. + 2.*f)*(x*x + y*y + z*z - w*w)*(x*x + y*y + z*z - w*w)*w*w - 0.4;
}

float  algebraic(vec3 p)
    {
        float x=p.x*p.x, y=p.y*p.y, z=p.z*p.z,a = 1.*1., r = 0.005;
        return ((x + y - a)*(x + y - a) + (z - 1.)*(z - 1.)) *
                ((z + y - a)*(z + y - a) + (x - 1.)*(x - 1.)) *
                ((x + z - a)*(x + z - a) + (y - 1.)*(y - 1.)) - r;

    }


float combo(vec3 p)
{
    //scale = 10.;
    if (csurf == 0.0)
        return eight(p);
    else
    if (csurf == 1.0)
        return barth2(p);
    else       
        return mix(eight(p), barth2(p), csurf);
}


float map(vec3 p) {
    //return sine(p);
    //return tooth(p);
    return barth2(p);
    //return algebraic(p);
    //return barth(p) - 0.1;
    //return barth(p)*barth(p) - 0.1;
    //return eight(p);
    //return zin(p)*zin(p) - 0.02;
    //return zin(p) - 0.005;
    //return digdong(p);
    //return cayley(p);
    //return umb2(p);
    //return chair( p);
    //return miter(p)*miter(p) - 0.05;
    //return miter(p) - 0.01;
    //return piri(p);
    //return roman(p);
    //return eight2(p);
    //return piri(p);
    //return isf(p);
    //return capsule2(p);
    //return min(sphere(p),sline0(p));
    //return floo(p);
    //return gyroide2(p);
    //return eggbox2(p);
    //return combo(p);
    //return gyroide(p);
    //return sline0(p);
    
    
    
    
}

vec3 calcNormal(in vec3 p) {
    const float eps = 0.0001;
    vec2 q = vec2(0.0, eps);
    vec3 res = vec3(map(p + q.yxx) - map(p - q.yxx), map(p + q.xyx) - map(p - q.xyx), map(p + q.xxy) - map(p - q.xxy));
    return normalize(res);
}

vec3 getPoint(vec3 a, vec3 b, float v0, float v1) {
    vec3 m;
    //binary search with  n iterations, n = newton
    for(int i = 0; i < newton; i++) {
        m = (a + b) * 0.5;
        float v = map(m);
        if(v == 0.)
            break;

        if(sign(v) * sign(v0) <= 0.) {
            v1 = v;
            b = m;
        } else {
            v0 = v;
            a = m;
        }
    }
    return m;
}

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
    csurf = glz();
    float dist_infin = 2.2;
    float hh = 4.2;

    vec3 light = normalize(vec3(0.0, 1.0, -2.5)); //light
    vec2 mo = 1.5*cos(0.5*rottime() + vec2(0,11));
    //vec2 mo = 1.5*cos(0.5*iTime + vec2(0,11));
    //if  (iMouse.z > 0.0)
    {
        mo = (-iResolution.xy + 2.0 * (iMouse.xy)) / iResolution.y;
    }
    vec3 ro = vec3(0.0, 0.0, hh); // camera
    //camera rotation
    ro.yz *= rot(mo.y);
    ro.xz *= rot(-mo.x - 1.57);

    const float fl = 1.5; // focal length
    float dist = dist_infin;

    vec3 bg = vec3(0.7, 0.7, 0.9)*0.6; //vec3(0.); //
    vec3 col1 = vec3(0.73, 0.59, 0.3);
    vec3 col2 = vec3(0.72, 0.01, 0.01);

    //antialiasing
    vec3 tot = vec3(0.0);
    for(int m = 0; m < AA; m++) for(int n = 0; n < AA; n++) {
            vec2 o = vec2(float(m), float(n)) / float(AA) - 0.5;
            vec2 p = (-iResolution.xy + 2.0 * (fragCoord + o)) / iResolution.y;
            vec3 rd = GetRayDir(p, ro, vec3(0, 0., 0), fl); //ray direction
            vec3 col = bg; // background  

            //STEP 1. Calculating bounding sphere
            float d = length(cross(ro, rd));
            if(d >= dist) {
                tot += col;
                continue;
            }
            /*
            STEP 2.
            ray tracing inside the bounding sphere, 
            searching for a segment with different signs of the function value 
            at the ends of the segment
            */
            float td = abs(dot(ro, rd));
            d = sqrt(dist * dist - d * d);
            vec3 pos0 = ro + rd * (td - d);
            vec3 pos1 = ro + rd * (td + d);
            vec3 rd0 = pos1 - pos0;
            vec3 pos = pos0;
            float val0 = map(pos0);
            for(float i = 1.; i < nn; i++) {
                pos1 = pos0 + rd0 * i / (nn - 1.);
                float val1 = map(pos1);
                if(sign(val0) * sign(val1) <= 0.) {
                    //different signs of the function value  at the ends of the segment
                    //STEP 3. binary search to clarify the intersection of a ray with a surface.
                    pos = getPoint(pos, pos1, val0, val1);
                    vec3 nor = calcNormal(pos);
                    col = col2;
                    
                    if(dot(rd, nor) < 0.0)
                        col = col1;
                    
                    //texture
                    
                    float tx = noise(pos*2.);
                    tx = fract(tx*5.);
                    tx = smoothstep(0., 0.01, tx-0.5);
                    col*=tx;
                    

                    //else break;    
                    vec3 R = reflect(light, nor);
                    float specular = pow(max(abs(dot(R, rd)), 0.), 25.);
                    float difu = abs(dot(nor, light));
                    col = col * (col * clamp(difu, 0., 1.0) + 0.5) + vec3(.5) * specular * specular;
                    col = sqrt(col);
                    break;
                }
                //if (sign(val1) < 0.) col = col2;
                val0 = val1;
                pos = pos1;
            }
            tot += col;
        }
    tot = tot / float(AA) / float(AA);
    //antialiasing
    fragColor = vec4(tot, 1.0);
}
/////=====================================================================================
void main() {
    vec4 fragColor = vec4(0);
    mainImage(fragColor, gl_FragCoord.xy);
    gl_FragColor = fragColor;
}