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
#define nn 50.
#define newton 5

float dist_infin =100.2;
float csurf = 0.;
float tmax = TAU;

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


vec3 trefoil(float t) {
    float a = 2.;
    return a * vec3(sin(t) + 2. * sin(2. * t), cos(t) - 2. * cos(2. * t), -0.5 * sin(3. * t));
}

vec3 curve(float t)
{
    return trefoil(t);
}
struct CU_LINE {
    vec3 ro;
    vec3 rd;
};

CU_LINE curve_line(float t)
{
    float h = 0.01;
    vec3 vt = curve(t);
    vec3 vth = curve(t + h);
    vec3 rd = normalize(vth - vt);
    return CU_LINE(vt, rd);
}

float cdist(float t, vec3 ro, vec3 rd, float r)
{
    CU_LINE val = curve_line(t);
    float x = (dot(val.ro, val.rd) - dot(ro, val.rd))/dot(rd, val.rd);
    if (x < 0.)
        return dist_infin;
    return length(ro + x*rd - val.ro) - r;
}

float rodist(float t, vec3 ro, vec3 rd)
{
    CU_LINE val = curve_line(t);
    return (dot(val.ro, val.rd) - dot(ro, val.rd))/dot(rd, val.rd);
}

vec3 calcNormal(float t, vec3 ro, vec3 rd)
{
    CU_LINE val = curve_line(t);
    float x = (dot(val.ro, val.rd) - dot(ro, val.rd))/dot(rd, val.rd);
    return normalize(ro + x*rd - val.ro);
}




float getPoint(float a, float b, float v0, float v1, vec3 ro, vec3 rd, float r) {
            float m;
            //binary search with  n iterations, n = newton
            
            for (int i = 0; i < newton; i++) {
                m = (a+b)*0.5;
                float v = cdist(m, ro, rd, r);
                //if (v == 0.)
                //    break;

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
    //csurf = glz();
    float    hh = 10.;
    float r = 0.6;
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
 

    //antialiasing
    vec3 tot = vec3(0.0);
    for(int m = 0; m < AA; m++) 
    for(int n = 0; n < AA; n++) {
            vec2 o = vec2(float(m), float(n)) / float(AA) - 0.5;
            vec2 p = (-iResolution.xy + 2.0 * (fragCoord + o)) / iResolution.y;
            vec3 rd = GetRayDir(p, ro, vec3(0, 0., 0), fl); //ray direction
            vec3 col = vec3(0.0);
            
            
            float d = dist_infin + 1., t = 0.;
            float pos0 = 0.;
            float val0 = cdist(pos0, ro, rd, r);
            for(float i = 1.; i < nn; i++) {
                float pos1 = tmax* i / (nn - 1.);
                float val1 = cdist(pos1, ro, rd, r);
                if (sign(val0)*sign(val1) <= 0.)
                {
                    //different signs of the function value  at the ends of the segment
                    //binary search to clarify the intersection of a ray with a surface.
                    float tt = (pos1 + pos0)/2.;
                    //float tt = getPoint(pos0, pos1, val0, val1, ro, rd, r);
                    float dd = rodist(tt, ro, rd);
                    if (dd < d)
                    {
                        d = dd;
                        t = tt;
                    }
                    
                }
                val0 = val1;
                pos0 = pos1;
            }
            if (d < dist_infin)
            {
                    vec3 nor = calcNormal(t, ro, rd);
                    float dif = clamp( dot(nor,light), 0.2, 1.0 );
                    float amb = 0.5 + 0.5*dot(nor,light2);
                    col = vec3(0.2,0.3,0.4)*amb + vec3(0.85,0.75,0.65)*dif;
                    // gamma        
                    col = sqrt( col );
            }
            
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