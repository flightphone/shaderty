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
//Collection of implicit surfaces. implicit surfaces, raytracing, binary search
/*
Rendering implicit surfaces. Using raytracing and binary searchy. 
Here, these same surfaces are obtained by creating grids using an algorithm 
3D Marching Cubes: https://flightphone.github.io/paramgeometry.html
*/
#define PI  3.14159265359
#define TAU 6.28318530718
#define rot(f) mat2(cos(f), -sin(f), sin(f), cos(f))
#define nn 64.
#define newton 6
float csurf = 0.;
float scale = 10.;
float npp = 60.;
float level = 0.9;

vec3 point(vec3 p) {
    return floor(p*npp)/npp;
}


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

float hesh (vec3 p) {
    return fract(sin(dot(p, vec3(127.1,311.7, 74.7))) * 43758.5453123);
}

float isf(vec3 p) {
    float x = p.x, y = -p.z, z = p.y;
    return (2. * y * (y * y - 3. * x * x) * (1. - z * z) + (x * x + y * y) * (x * x + y * y) - (9. * z * z - 1.) * (1. - z * z));// IMPLICIT SURFACE Function
}



float capsule(vec3 p)
{
    float h = 1.8;
    float x = p.x, y = p.y, z = 0.;
    if (p.z > h)
        z = p.z - h;
    if (p.z < -h)
        z = p.z + h; 

    return x*x + y*y + z*z - 0.4*0.4;       

}

float sphere(vec3 p) {
    float x = p.x, y = p.y, z = p.z;
    float r = 1.5;
    return (x * x + y * y + z * z - r*r + 0.06*(sin(35.*x)+sin(20.*y)+sin(25.*z)));
}

float plane(vec3 p)
{
    return dot(p, vec3(0., 0., 1.));
}

float sphere2(vec3 p)
{
    return length(p-vec3(0., 0., .5)) - 0.6;
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
    scale = 5.;
    return gyroide(p)*(gyroide(p)) - 0.3;
}

float combo(vec3 p)
{
    //scale = 10.;
    if (csurf == 0.0)
        return eggbox(p);
    else
    if (csurf == 1.0)
        return sphere2(p);
    else       
        return mix( eggbox(p), sphere2(p), csurf);
}

float map(vec3 p) {
    return isf(p);
    //return capsule(p);
    //return sphere(p);
    //return gyroide2(p);
    //return eggbox(p);
    //return combo(p);
    //return gyroide(p);
    //return sphere2(p)*eggbox(p) - 0.01;
    
    
    
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

    vec3 bg = vec3(0.);
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
                    //points
                    vec3 pp = point(pos);
                    float fil = hesh(pp);
                    if (fil > level)
                    {
                        col = vec3(1.);
                        break;
                    }
                    
                    
                    /*
                    vec3 nor = calcNormal(pos);
                    col = col2;
                    if(dot(rd, nor) < 0.0)
                        col = col1;
                    vec3 R = reflect(light, nor);
                    float specular = pow(max(abs(dot(R, rd)), 0.), 25.);
                    float difu = abs(dot(nor, light));
                    col = col * (col * clamp(difu, 0., 1.0) + 0.5) + vec3(.5) * specular * specular;
                    col = sqrt(col);
                    break;
                    */
                }
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