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
#define newton 10

float dist_infin =2.2;
float csurf = 0.;
mat3 rotateX(float f) {
    return mat3(vec3(1.0, 0.0, 0.0), vec3(0.0, cos(f), -sin(f)), vec3(.0, sin(f), cos(f)));
}
mat3 rotateY(float f) {
    return mat3(vec3(cos(f), 0.0, sin(f)), vec3(0.0, 1.0, 0.0), vec3(-sin(f), 0.0, cos(f)));
}

float sdRound(vec3 p, float r, float f)
{
    float d = abs(length(p.xy) - r*cos(f));
    d = length(vec2(p.z-r*sin(f), d));
    return d;
}
//https://www.shadertoy.com/view/DtdBR4
float sinewave(vec3 p, float a, float b, float m, float n)
{
    float fi = atan(p.y, p.x); //aafi(p.xy)
    float w = dist_infin;
    for (float i = 0.; i < 10.; i++ )
    {
        if (mod(i,m) == 0.0 && i > 0.)
            break;
        float t = (fi + TAU*i)/m;    
        float wt = abs(p.z - b*sin(n*t));
        w = min(w, wt);
    }
    float r = length(vec2(length(p.xy) - a, w))/2.0;
    return r - 0.03;
}

//https://www.shadertoy.com/view/mttfzl
float sdLonLat(vec3 p, float r)
{
        #define ll 20.
        float fi = atan(p.x, p.y);
        fi += step(p.y, 0.0)*TAU;
        float ln = floor(fi/TAU*ll);
        float l1 = ln * TAU/ll;
        float l2 = l1 + TAU/ll;
        float d = min(
            sdRound(rotateX(l1)*rotateY(PI/2.)*p, r, 0.), 
            sdRound(rotateX(l2)*rotateY(PI/2.)*p, r, 0.));
        
        fi = atan(p.z, length(p.xy));
        float mm = ll/4.0;
        ln = floor(abs(fi)/PI*2.0*mm);
        l1 = ln*PI/2.0/mm;
        l2 = l1 + PI/2.0/mm;
        float d2 = min(sdRound(p, r, l1*sign(p.z)), sdRound(p, r, l2*sign(p.z)));
        d = min(d2, d);
        return d - 0.05;
}


float glz() {
    float t = iTime / 4.;
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
vec3 S = vec3(3,2,0); 
float knot(vec3 q) {
        q = vec3( length(q.xz)-.5, q.y, atan(q.z,q.x) ), // cylindrical coordinates
        q = vec3( length(q.xy), atan(q.y,q.x), q.z );             // torus coordinates
        return length( vec2(q.x-.2, ( mod(S.x*q.z-S.y*q.y,6.28)-3.14 ) /16. )) - 0.1; // solenoïd
}
float sdTorus( vec3 p)
{
  vec2 t = vec2(0.6, 0.15);
  float d = (p.x*p.x + p.y*p.y + p.z*p.z + t.x*t.x - t.y*t.y);
  return d*d - 4.*t.x*t.x*(p.x*p.x + p.y*p.y); 
  //vec2 q = vec2(length(p.xz)-t.x,p.y);
  //return length(q)-t.y;
}
float map( in vec3 p )
{
    //return knot(p.xzy);
    //return sdTorus(p);
  if (csurf == 0.0)
        return knot(p.xzy);
    else
    if (csurf == 1.0)
        return sdTorus(p);
    else       
        return sdTorus(p)*csurf + (1.-csurf)*knot(p.xzy);

    /*
    if (csurf == 0.0)
        return sinewave(p, .9, .4, 3.5, 9.);
    else
    if (csurf == 1.0)
        return sdLonLat(p, 1.0);
    else       
        return sdLonLat(p, 1.0)*csurf + (1.-csurf)*sinewave(p, .9, .4, 3.5, 9.);
    */    
    //return pown(pos);
    //return queen(pos);
}


vec3 calcNormal(in vec3 p) {
    const float eps = 0.0001;
    vec2 q = vec2(0.0, eps);
	vec3 res =  vec3(map(p+q.yxx) - map(p-q.yxx), 
			    map(p+q.xyx) - map(p-q.xyx),
			    map(p+q.xxy) - map(p-q.xxy));
    return normalize(res);
}

vec3 getPoint(vec3 a, vec3 b, float v0, float v1) {
            vec3 m;
            //binary search with  n iterations, n = newton
            for (int i = 0; i < newton; i++) {
                m = (a+b)*0.5;
                float v = map(m);
                if (v == 0.)
                    break;

                if (sign(v) * sign(v0) <= 0.) {
                    v1 = v;
                    b = m;
                }
                else {
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
    dist_infin = 1.2;
    float    hh =1.6;
        
    vec3 light = normalize(vec3(0.0, 1.0, 1.0)); //light
    vec3 light2 = normalize(vec3(0.0, -1.0, 2.5)); //light
    vec2 mo = 1.5*cos(0.3*iTime + vec2(0,11));
    //if  (iMouse.z > 0.0)
    {
        mo = (-iResolution.xy + 2.0 * (iMouse.xy)) / iResolution.y;
    }
    vec3 ro = vec3(0.0, 0.0, hh ); // camera
    //camera rotation
    ro.yz *= rot(mo.y);
    ro.xz *= rot(-mo.x-1.57);

    const float fl = 1.5; // focal length
    float dist = dist_infin;

    //antialiasing
    vec3 tot = vec3(0.0);
    for(int m = 0; m < AA; m++) 
    for(int n = 0; n < AA; n++) {
            vec2 o = vec2(float(m), float(n)) / float(AA) - 0.5;
            vec2 p = (-iResolution.xy + 2.0 * (fragCoord + o)) / iResolution.y;
            vec3 rd = GetRayDir(p, ro, vec3(0, 0., 0), fl); //ray direction
            vec3 col = vec3(0.0);
            
            
            //STEP 1. Calculating bounding sphere
            float d = length(cross(ro, rd));
            if (d >= dist)
            {
                 continue;
            }
            /*
            STEP 2.
            ray tracing inside the bounding sphere, 
            searching for a segment with different signs of the function value 
            at the ends of the segment
            */
            float td = abs(dot(ro, rd));
            d = sqrt(dist*dist - d*d);
            vec3 pos0 = ro + rd * (td - d);
            vec3 pos1 = ro + rd * (td + d);
            vec3 rd0 = pos1 - pos0;
            vec3 pos = pos0;
            float val0 = map(pos0);
            for(float i = 1.; i < nn; i++) {
                pos1 = pos0 + rd0 * i / (nn - 1.);
                float val1 = map(pos1);
                if (sign(val0)*sign(val1) <= 0.)
                {
                    //different signs of the function value  at the ends of the segment
                    //STEP 3. binary search to clarify the intersection of a ray with a surface.
                    
                    pos = getPoint(pos, pos1, val0, val1);
                    vec3 nor = calcNormal(pos);
                    float dif = clamp( dot(nor,light), 0.2, 1.0 );
                    float amb = 0.5 + 0.5*dot(nor,light2);
                    col = vec3(0.2,0.3,0.4)*amb + vec3(0.85,0.75,0.65)*dif;
                    break;
                }
                val0 = val1;
                pos = pos1;
            }
            // gamma        
            col = sqrt( col );
            tot += col;
        }
    tot = tot / float(AA)/ float(AA);
    //antialiasing
    fragColor = vec4(tot, 1.0);
}
/////=====================================================================================
void main() {
    vec4 fragColor = vec4(0);
    mainImage(fragColor, gl_FragCoord.xy);
    gl_FragColor = fragColor;
}