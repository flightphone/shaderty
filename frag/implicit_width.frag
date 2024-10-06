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
float hash(float p) {
    p = fract(p * 0.011);
    p *= p + 7.5;
    p *= p + p;
    return fract(p);
}
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

/////=====================================================================================

#define PI  3.14159265359
#define TAU 6.28318530718
#define rot(f) mat2(cos(f), -sin(f), sin(f), cos(f))
#define nn 128.
#define newton 5

float sdSphere2(vec3 p, float R) {
    float res = length(p) - R;
    if(p.z < 0.) {
        float t = length(p.xy) - R;
        if(t <= 0.) {
            res = -p.z;
        } else {
            res = length(vec2(p.z, t));
        }
    }
    return res;
}

float acr(vec3 p, float R, float h, float l) {
    float res = 0.;
    if(p.z >= 0.)
        res = length(vec2(p.x, max(p.z - h, 0.))) - R;
    else
        res = length(vec2(p.z, abs(p.x) - R));
    res = length(vec2(max(abs(p.y) - l, 0.), res));
    return res;
}

float sdCutHollowSphere(vec3 p, float r, float h, float t) {
  // sampling independent computations (only depend on shape)
    float w = sqrt(r * r - h * h);

  // sampling dependant computations
    vec2 q = vec2(length(p.xz), p.y);
    return ((h * q.x < w * q.y) ? length(q - vec2(w, h)) : abs(length(q) - r)) - t;
}

float pown(vec3 p) {
    float k = 2.5, x = p.x * k, y = p.y * k, z = p.z * k, a = 2.5;
    float t = (x * x + y * y) * 4.5 - (a - z) * z * z;
    float f = .5;
    t += pow(min(f + z, 0.), 4.);
    return t;
}

float ring(vec3 p, float h, float R, float r) {
    float d = abs(length(p.xy) - R) - r;
    vec2 w = vec2(d, abs(p.z) - h);
    return min(max(w.x, w.y), 0.0) + length(max(w, 0.0));
}

float quin(vec3 p) {
    float k = 3., x = p.x * k, y = p.y * k, z = p.z * k;
    float res = (x * x + y * y) - z * z * .1 - 0.5 + pow(min(z + 2., 0.), 2.) + pow(max(z - 4., 0.), 6.) - 2. * max(z - 4., 0.) + 1. * min(z + 2., 0.);
    z += 3.;
    float res2 = ring(vec3(x, y, z), 0.5, 1.2, 0.1) - 0.1;
    return min(res2, res);
}

float pluspole(vec3 p, float r) {
    return pow(5.0*max(length(p) - r, 0.0), 2.);
}

float pluspole(vec2 p, float r) {
    return pow(2. * max(length(p) - r, 0.0), 10.);
}

float pluspole3(vec3 p, float r, float scale) {
    p *= scale;
    return pow(max(abs(p.x) - r, 0.0), 3.) +
        pow(max(abs(p.y) - r, 0.0), 3.) +
        pow(max(abs(p.z) - r, 0.0), 3.);
}

float gyroide(vec3 p, float scale) {
    float x = p.x * scale, y = p.y * scale, z = p.z * scale;
    return (cos(x) * sin(y) + cos(y) * sin(z) + cos(z) * sin(x));
}

float gyroide2(vec3 p) {
    float scale = 4.;
    float x = p.x * scale, y = p.y * scale, z = p.z * scale;
    return gyroide(p, scale) * gyroide(p, scale) - 0.1 + pluspole3(p, PI * 2. - 0.3, scale);
}

float sch(vec3 p, float scale) {
    p *= scale;
    return cos(p.x) + cos(p.y) + cos(p.z);
}
float sch2(vec3 p) {
    float scale = 5.;
    return sch(p, scale) * sch(p, scale) - 0.05 + pluspole3(p, 3. * PI, scale);

}

float goursat(vec3 p) {
    //https://mathcurve.com/surfaces.gb/goursat/goursat.shtml
    float x = p.x, y = p.y, z = p.z;
    float a = 0.5;
    return (2. * (x * x * y * y + x * x * z * z + z * z * y * y) - a * a * (x * x + y * y + z * z) - a * a * a * a);
}

float goursat2(vec3 p) {
    //https://mathcurve.com/surfaces.gb/goursat/goursat.shtml
    return goursat(p) * goursat(p) - 0.5 + +pluspole3(p, 1.8, 1.);
    ;
}

float cylinder(vec3 p, float R, float h) {

    float d = length(p.xy) - R;
    vec2 w = vec2(d, abs(p.z) - h);
    return min(max(w.x, w.y), 0.0) + length(max(w, 0.0));
}

float smin(float d1, float d2, float k) {
    float h = clamp(0.5 + 0.5 * (d2 - d1) / k, 0.0, 1.0);
    return mix(d2, d1, h) - k * h * (1.0 - h);
}

float smin_out(float d1, float d2, float k) {
    float h = clamp(0.5 - 0.5 * (d2 + d1) / k, 0.0, 1.0);
    return mix(d2, -d1, h) + k * h * (1.0 - h);
}

float cylinder2(vec3 p) {
    float t = cylinder(p, 0.7, 2.);
    t = smin(t, cylinder(p.yzx, 0.7, 2.), 0.3);
    t = smin(t, cylinder(p.xzy, 0.7, 2.), 0.3);   
    /*
    t = max(t, -cylinder(p, 0.3, 3.) );
    t = max(t, -cylinder(p.yzx, 0.3, 3.) );
    t = max(t, -cylinder(p.xzy, 0.3, 3.) );
    */
    t = smin_out(cylinder(p, 0.3, 3.), t, 0.05);
    t = smin_out(cylinder(p.yzx, 0.3, 3.), t, 0.05);
    t = smin_out(cylinder(p.xzy, 0.3, 3.), t, 0.05);

    return t - 0.05;
}

vec3 getSg(vec3 p, float nseg) {
    float fi = mod(atan(p.y, p.x), TAU);
    fi = mod(fi + PI / nseg, TAU);
    float n = floor(fi / TAU * nseg);
    p.xy *= rot(-n * TAU / nseg);
    return p;
}
float sdSegment(in vec2 p, in vec2 a, in vec2 b) {
    vec2 pa = p - a, ba = b - a;
    float h = clamp(dot(pa, ba) / dot(ba, ba), 0.0, 1.0);
    return length(pa - ba * h);
}

float gou(vec3 p) {
    float k = 0.7, x = p.x * k, y = p.y * k, z = p.z * k;
    return x * x * x * x + y * y * y * y + z * z * z * z - x * x - y * y - z * z + 0.4;
}

float dish(vec3 p) {
    float k = 0.7, x = p.x * k, y = p.y * k, z = p.z * k;
    return x * x * x * x + 0.25 * y * y * y * y + z * z * z * z + 2. * x * x * z * z - 3. * y * (x * x + z * z) + 2. * y * y - 0.01;
}

float klein(vec3 p) {
    float k = 2., x = p.x * k, y = p.y * k, z = p.z * k, a = (x * x + y * y + z * z + 2. * y - 1.);
    return a * ((a - 4. * y) * (a - 4. * y) - 8. * z * z) + 16. * x * z * (a - 4. * y) - 0.1;
}

float disk(vec3 p) {

    float k = 1., x = p.x * k, y = p.y * k, z = p.z * k;
    float r = 1., a = max(length(p.xy) - r, 0.0);//a = max(x*x + y*y - r*r, 0.); 
    float t = (x + 1.) / 2.0;
    a *= a;
    float b = ((p.z - 2. * (-2. * t * t * t + 3. * t * t)));
    b *= b;

    return a + b - 0.01;
}
float archytas_curve(vec3 p) {
    float k = 0.5, x = p.x * k, y = p.y * k, z = p.z * k;
    float a = 1., 
    //u1 = length(p) - a* length(p.xy), 
    u1 = (x * x + y * y + z * z) * (x * x + y * y + z * z) - a * a * (x * x + y * y), u2 = x * x + y * y - a * x;
    u1 *= u1;
    u2 *= u2;
    //return u1 + 0.1*u2 - 0.01;
    return u1 + u2 - 0.01;
}

float plucker(vec3 p) {
    float k = .5, x = p.x * k, y = p.y * k, z = p.z * k, a = 1.;
    return z * (x * x + y * y) - a * (x * x - y * y);
}
float sphplu(vec3 p) {
    float a = 1., k = 0.5;
    return length(p) * k - a * a;
}

float disk2(vec3 p) {

    float k = 1., x = p.x * k, y = p.y * k, z = p.z * k;
    float r = 2.5; //a = max(length(p.xy) - r, 0.), t = (x+1.)/2.0;
    float a = (length(p.xy) - r);
    a *= a;

    return a + z * z - 0.01;
}

float heelb(vec3 p) {
    float k = 1., x = p.x * k, y = p.y * k, z = p.z * k, d = 0., r = 0.5, w = 1.;
    if(x < 0.) {
        d = length(vec2(x, max(abs(y) - r, 0.)));
    } else {
        d = length(vec2(max(x - w, 0.), abs(y) - r)) - r;
    }

    return d * d + (z - 1.) * (z - 1.) - 0.01;
}

float insole(vec3 p) {
    float k = 1., x = p.x * k, y = p.y * k, z = p.z * k;
    float w = 1.5, r = 0.2;
    x += w;
    float z0 = 1.0 - smoothstep(0.1, w * 0.85, x);
    r = r + clamp(x, 0., w) * 0.15;
    
    float d = max(length(vec2(x - clamp(x, 0., w), y)) - r, 0.);// + pow(1.0 * min(x, 0.), 2.);
    return d * d + (z - z0) * (z - z0) - 0.002;
}

float streep(vec3 p, float h, float r) {
    float d = abs(length(vec2(p.y, p.z)) - r) + pow(10.*max(abs(p.x) - h, 0.), 6.) + 
    pow(10.*min(p.z, 0.), 6.);
    //vec2 w = vec2(d, abs(p.z) - h);
    //float res =  min(max(w.x, w.y), 0.0) + length(max(w, 0.0));
    return d;
}

float galoshes(vec3 p)
{
    p.yz *= rot(PI/2.);
    float k = 0.6, x = p.x * k, y = p.y * k, z = p.z * k;
    float r = 0.8, w = 2., d = 0., h = 0.5;
    z += h*.5;
    x += w*(1. + h)/2.0;
    //r = r * (1. - z*z*0.5);
    w = w * (1. + z*0.3);
    //d = length(vec2(max(abs(x) - w, 0.), y)) - r 
    d = abs(y) - r + pow(5.0 * min(x, 0.), 2.) + pow(max(x-w, 0.), 2.);
    d = d + pow(3. * min(z, 0.), 2.);
    float eps = 0.08;
    float d2 = abs(y) - (r-eps) + pow(5.0 * min(x - eps, 0.), 2.) + pow(max(x-(w-eps), 0.), 2.);
    d2 = d2 + pow(3. * min(z-eps, 0.), 2.);
    //float d2 = length(vec2(max(abs(x) - (w-eps), 0.), y)) - (r-eps) + pow(10.0 * min(x-eps, 0.), 6.);
    //d2 = d2 + pow(10. * min(z-eps, 0.), 2.);
    return d*d2 + pow(1.*max(z-h, 0.), 2.) + pow(10.*min(x-4.*eps, 0.)*max(z-h, 0.), 6.) ;
}

float heel(vec3 p) {

    float k = 0.45, x = p.x * k, y = p.y * k, z = p.z * k;
    float r = 0.05, w = 0.05, d = 0.;
    z += .5;
    x += .64;
    r = r * (2. * z * z + 1.);
    w = w * (2. * z * z + 1.);
    d = length(vec2(max(abs(x) - w, 0.), y)) - r + pow(10.0 * max(x, 0.), 6.);
    d = d + pow(10. * min(z, 0.), 4.);
    d = d + pow(10. * max(z - 1., 0.), 4.);

    float d2 = insole(vec3(x - 1.35, y, z - 0.1));
    float d3 = streep(vec3(x-1.35, y, z-0.15), 0.1, 0.4)-0.02;
    /*
    float d3 = abs(length(vec2(p.y, p.z)) - r) + pow(5.*max(abs(p.x) - h, 0.), 4.) + 
    pow(10.*min(p.z, 0.), 6.);
    */
    //float d = abs(length(vec2(p.y, p.z)) - r) + pow(max(abs(p.x) - h, 0.), 2.);
    d =  d * d2 - 0.0002;
    return d*d3;

}

float coti(vec3 p)
{
    float k = 1., x = p.x * k, y = p.y * k, z = p.z * k, a = 0.6/TAU,
    dz = TAU*a;
    float f = mod(atan(y, x), TAU);
    z -= floor(z/dz)*dz; 
    //return (f*a - z);
    return (f*a - z);
    //return (f/a - z)*(f/a - z) - 0.1;
}

float coti2(vec3 p)
{
    //https://mathcurve.com/surfaces.gb/helicoiddroit/helicoiddroit.shtml
    float a = 0.5*TAU, f = p.z*a;
    return dot(p.xy, vec2(cos(f), sin(f)));
}

float cyl(vec3 p)
{
    float k = 1., x = p.x * k, y = p.y * k, z = p.z;
    return length(vec2(x, y)) - .3;
}

float coti3(vec3 p)
{
    p*=1.;
    return coti2(p)*coti2(p) + pluspole(p, 2.6) - 0.01+ cyl(p)*cyl(p);
}


float map(vec3 p) {
    //https://mathcurve.com/surfaces.gb/orthobicycle/orthobicycle.shtml
    
    
    //return coti2(p);
    //return coti2(p) + cyl(p)*cyl(p)+ pluspole(p, 2.6);
    //return coti2(p)*coti2(p) + pluspole(p, 2.6) - 0.1;
    //return coti2(p)*coti2(p) + pluspole(p.xy, 1.6) + pluspole(p, 2.6) - 0.1;
    //return coti3(p);
    //return galoshes(p);
    //return insole(p);
    //return heel(p);
    //return streep(p, 1., 0.5)-0.1;

    return sch2(p);
    //return dish(p);
    //return gou(p);
    //return disk(p);

    //return archytas_curve(p);
    //return plucker(p)*plucker(p) + sphplu(p)*sphplu(p) - 0.05;
    //return abs(sphplu(p)) - 0.01;
    //return abs(plucker(p)) - 0.01;
    //return max(plucker(p)*plucker(p),sphplu(p)*sphplu(p)) - 0.2;
    //return disk2(p);
    //return gyroide2(p);
    //return goursat2(p);
    //return cylinder2( p);
    //return umb(p);
    //return quin(p);
    //return ring(p)-0.1;
    //return pown(p);
    //return klein(p);
    //return sdSphere2(p, 2.) - 0.1;
    //return (sdSphere2(p, 2.) - 0.01) * (sdSphere2(p, 1.7) - 0.01);
    //return sdCutHollowSphere(p, 2., 0.4, 0.05);

    //return acr(p, 0.5, 1., 1.) - 0.05;
    //return (acr(p, 0.5, 1., 1.)-0.1)*(acr(p.yxz, 0.5, 1., 1.) - 0.1)-0.001;
    //return max((acr(p, 0.5, 1., 1.)-0.01),(acr(p.yxz, 0.5, 1., 1.)-0.01))-0.2;

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
    //csurf = glz();
    float dist_infin = 3.;
    float hh = 6.;

    vec3 light = normalize(vec3(0.0, 1.0, -2.5)); //light
    //vec2 mo = 1.5*cos(0.5*rottime() + vec2(0,11));
    vec2 mo = 1.5 * cos(0.5 * iTime + vec2(0, 11));
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

    vec3 bg = vec3(0.7, 0.7, 0.9) * 0.6; //vec3(0.); //
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
            float sect = 0.;
            for(float i = 1.; i < nn; i++) {
                pos1 = pos0 + rd0 * i / (nn - 1.);
                float val1 = map(pos1);
                if(sign(val0) * sign(val1) <= 0.) {
                    //different signs of the function value  at the ends of the segment
                    //STEP 3. binary search to clarify the intersection of a ray with a surface.
                    pos = getPoint(pos, pos1, val0, val1);
                    sect = 1.0;
                    break;
                }

                val0 = val1;
                pos = pos1;
            }
            if(sect == 1.0) {
                vec3 nor = calcNormal(pos);
                col = col2*col2;

                if(dot(rd, nor) < 0.0)
                    col = col1;

                /*
                if (
                    length(pos.xy - vec2(0.8, 0.3)) < 0.2 ||
                    length(pos.xy - vec2(-0.8, -0.3)) < 0.2
                )
                    col = vec3(0.,0., 1.);
                */    
                    //texture
                    /*
                    float tx = noise(pos*2.);
                    tx = fract(tx*5.);
                    tx = smoothstep(0., 0.01, tx-0.5);
                    col*=tx;
                    */
                    /*
                    float tx = noise(pos*2.);
                    tx = fract(tx*20.);
                    col*=tx;
                    */

                vec3 R = reflect(light, nor);
                float specular = pow(max(abs(dot(R, rd)), 0.), 16.);
                float difu = abs(dot(nor, light));
                col = col * (col * clamp(difu, 0., 1.0) + 0.5) + vec3(.5) * specular;
                float fre = pow(clamp(dot(nor, rd) + 1., .0, 1.), 3.); // Fresnel, for some mild glow.
                col += vec3(.1, .1, 0.1) * fre; //?
                col = sqrt(col);
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