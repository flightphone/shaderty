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


#define nn 128.
#define newton 8
#define AA 1
float csurf = 0.;
const float eps = 0.0001;
float dist_infin =2.2;

vec3 sdfColor;
vec3 resColor;
vec3 col1 = vec3(0.3764, 0.8196, 0.3725);
vec3 col2 = vec3(0.8117, 0.1764, 0.8078);


vec3 q1, q2, q3, q4, q5;
float potential(vec3 p) {
        float res = 0.;
        res += 0.55/length(p - q1);
        res += 0.55/length(p - q2);
        res += 0.55/length(p - q3);
        res += 0.55/length(p - q4);
        res += 0.55/length(p - q5);
        return 4.5 - res;
    }
float heart(vec3 p)
    //https://www.shadertoy.com/view/XtXGR8
    {
        float x = p.x, y = p.z, z = p.y;
        float s = sin(iTime * 5.0);
        s *= s;
        s *= s;
        s *= 0.2;
        float a = (x*x + 2.25*y*y + z*z - (1.+s));
        return a*a*a - (x*x + 0.1125*y*y)*z*z*z;
    }    

float riemann(vec3 p) {
    float x = p.x, y = p.y, z = p.z;
    float a = 1.;
    float d = (x * x + y * y - a * a) * z - a * x;
    return d;
}


float holed2(vec3 p) {
    float x = p.x, y = p.y, z = p.z;
    float x2 = x * x, y2 = y * y, z2 = z * z;
    float u = (x2) * (1. - x2) - y2;
    return u * u + z2 - 0.01;
}

float kummerj(vec3 p) {
        float x = p.x, y = p.y, z = p.z;
        return x * x * x * x + y * y * y * y + z * z * z * z - 5. * (x * x * y * y + y * y * z * z + z * z * x * x) + 56. * x * y * z -
            20. * (x * x + y * y + z * z) + 16.;
        
    }

float gayley(vec3 p) {
    //https://mathcurve.com/surfaces.gb/cayley/cayley.shtml
    float x = p.x, y = p.y, z = p.z,  k = 2., a = 1.;
    return ((x + y + z - a) * (x * y + y * z + z * x) - k * x * y * z);
}


float goursat(vec3 p) {
    //https://mathcurve.com/surfaces.gb/goursat/goursat.shtml
    float x = p.x, y = p.y, z = p.z; 
    float a = 1.;
    return (2. * (x * x * y * y + x * x * z * z + z * z * y * y) - a * a * (x * x + y * y + z * z) - a * a * a * a);
}

float gyroide(vec3 p) {
    float x = p.x, y = p.y, z = p.z; 
    return cos(x) * sin(y) + cos(y) * sin(z) + cos(z) * sin(x);
}

float isf(vec3 p) {
   float x = p.x, y = p.y, z = p.z; 
   return (2. * y * (y * y - 3. * x * x) * (1. - z * z) + (x * x + y * y) * (x * x + y * y) - (9. * z * z - 1.) * (1. - z * z));// IMPLICIT SURFACE Function
}
float desimp(vec3 p)
{
    float res = p.x*p.x + p.y*p.y + p.z*p.z + sin(4.*p.x) + sin(4.*p.y) + sin(4.*p.z) - 1.11;
    return (res);
}
float  algebraic(vec3 p)
    {
        float x=p.x*p.x, y=p.y*p.y, z=p.z*p.z,a = 1.*1., r = 0.005;
        return ((x + y - a)*(x + y - a) + (z - 1.)*(z - 1.)) *
                ((z + y - a)*(z + y - a) + (x - 1.)*(x - 1.)) *
                ((x + z - a)*(x + z - a) + (y - 1.)*(y - 1.)) - r;

    }


float cassinian(vec3 p) {
    
        float res = 1., b = 4.15, h0 = 3., h1 = 5.0, r = (h0 + h1)/2. - (h1 - h0)/2.*cos(2.*iTime);
        
        for (float i = 0.; i < 2.; i++)
            for (float j = 0.; j < 2.; j++)
                for (float k = 0.; k < 2.; k++)
                    res *= length(p - vec3(r * (i - 0.5), r * (j - 0.5), r * (k - 0.5)));
        return res - b * b * b * b * r * r * r * r;
    }    


float map(vec3 p) {
    if (csurf == 0.)
        return desimp(p);
    if (csurf == 1.)
        return algebraic(p);
    if (csurf == 2.)
        return isf(p);
    if (csurf == 3.)
        return cassinian(p);
    if (csurf == 4.)
        return gyroide(p);
    if (csurf == 5.)
        return goursat(p);
    if (csurf == 6.)
        return gayley(p);
    if (csurf == 7.)
        return potential(p);
    if (csurf == 8.)
        return kummerj(p);
    if (csurf == 9.)
        return holed2(p);
    if (csurf == 10.)
        return riemann(p);  
    if (csurf == 11.)
        return heart(p);  
}


vec3 calcNormal(in vec3 p) {
    vec2 q = vec2(0.0, eps);
	vec3 res =  vec3(map(p+q.yxx) - map(p-q.yxx), 
			    map(p+q.xyx) - map(p-q.xyx),
			    map(p+q.xxy) - map(p-q.xxy));
    return normalize(res);
}

vec3 getPoint(vec3 a, vec3 b, float v0, float v1) {
            vec3 m;
            //binary find with  n iterations, n = newton
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


vec3 calccolor(vec3 col_in, vec3 backcol, vec3 rd, vec3 light1, vec3 light2, vec3 nor) {
    vec3 col = col_in;
    float d = dot(rd, nor);
    if(d < 0.0)
        col = backcol;

    nor *= -sign(d);
    float difu1 = dot(nor, light1);
    float difu2 = dot(nor, light2);
    float difu = max(difu1, difu2);
    col *= clamp(difu, 0.3, 1.0);
    return col;
}

void initQ()
{
    float x, y, z, i, time = 2.*iTime;
    i = 1.;
    x = sin(i + 0.62 * time * (1.03 + 0.5 * cos(1.51 * i)))*0.3;
    y = (cos(i + 1.17 * time * cos(1.22 + 1.1424 * i)))*0.8; // dip into the floor //Math.abs
    z = cos(i + 0.51 * time * 0.1 * sin((0.92 + 1.43 * i)))*0.3;
    q1 = vec3(x, y, z);

    i = 2.;
    x = sin(i + 0.62 * time * (1.03 + 0.5 * cos(1.51 * i)))*0.3;
    y = (cos(i + 1.17 * time * cos(1.22 + 1.1424 * i)))*0.8; // dip into the floor //Math.abs
    z = cos(i + 0.51 * time * 0.1 * sin((0.92 + 1.43 * i)))*0.3;
    q2 = vec3(x, y, z);

    i = 3.;
    x = sin(i + 0.62 * time * (1.03 + 0.5 * cos(1.51 * i)))*0.3;
    y = (cos(i + 1.17 * time * cos(1.22 + 1.1424 * i)))*0.8; // dip into the floor //Math.abs
    z = cos(i + 0.51 * time * 0.1 * sin((0.92 + 1.43 * i)))*0.3;
    q3 = vec3(x, y, z);

    i = 4.;
    x = sin(i + 0.62 * time * (1.03 + 0.5 * cos(1.51 * i)))*0.3;
    y = (cos(i + 1.17 * time * cos(1.22 + 1.1424 * i)))*0.8; // dip into the floor //Math.abs
    z = cos(i + 0.51 * time * 0.1 * sin((0.92 + 1.43 * i)))*0.3;
    q4 = vec3(x, y, z);

    i = 5.;
    x = sin(i + 0.62 * time * (1.03 + 0.5 * cos(1.51 * i)))*0.3;
    y = (cos(i + 1.17 * time * cos(1.22 + 1.1424 * i)))*0.8; // dip into the floor //Math.abs
    z = cos(i + 0.51 * time * 0.1 * sin((0.92 + 1.43 * i)))*0.3;
    q5 = vec3(x, y, z);

}

void mainImage(out vec4 fragColor, in vec2 fragCoord) {
    
    csurf = mod(floor(iTime/6.) , 12.);
    //csurf  = 1.;
    float hh = 5.;
    if (csurf == 0.)
    {
        dist_infin = 2.2;
        hh = 4.;
        //desimp(p);
    }    
    if (csurf == 1.)
    {
        dist_infin = 2.;
        hh = 4.5;
        //algebraic(p);
    }    
        
    if (csurf == 2.)
    {
        dist_infin = 2.;
        hh = 4.;
        //isf(p);
    }    
        
    if (csurf == 3.)
    {
        dist_infin = 5.;
        hh = 9.;
        //cassinian(p);
    }    
        
    if (csurf == 4.)
    {
        dist_infin = 10.;
        hh = 26.;
        //gyroide(p);
    }    
        
    if (csurf == 5.)
    {
        dist_infin = 3.2;
        hh = 7.;
        //goursat(p);
    }    
        
    if (csurf == 6.)
    {
        dist_infin = 3.2;
        hh = 7.;
        //gayley(p);
    }    
        
    if (csurf == 7.)
    {
        initQ();
        dist_infin = 1.5;
        hh = 3.;
        //potential(p);
    }   

    if (csurf == 8.)
    {
        dist_infin = 7.5;
        hh = 15.;
        //kummerj(p); 
    }
        
    if (csurf == 9.)
    {
        dist_infin = 1.2;
        hh = 2.4;
        //holed2
    }   

    if (csurf == 10.)
    {
        dist_infin = 4.;
        hh = 8.;
        //holed2
    }   
    if (csurf == 11.)
    {
        dist_infin = 1.8;
        hh = 4.;
        //heart
    }   
        
    vec3 light = normalize(vec3(0.0, 1.0, -2.5)); //light
    vec3 light2 = normalize(vec3(0.0, -1.0, 2.5)); //light
    vec2 mo = vec2(0.0, 0.0);
    //if  (iMouse.z > 0.0)
    {
        mo = (-iResolution.xy + 2.0 * (iMouse.xy)) / iResolution.y;
    }
    vec3 ro = vec3(0.0, 0.0, hh ); // camera
    //camera rotation
    ro.yz *= rot(mo.y * PI);
    ro.xz *= rot(-mo.x * TAU);

    const float fl = 1.5; // focal length
    float dist = dist_infin;

    vec3 
    b1 = vec3(0.0509, 0.2980, 0.4705), 
    b2 = vec3(0.3764, 0.7529, 0.8784), 
    bg = mix(b1, b2, vec3((1.0 - abs(fragCoord.x - iResolution.x / 2.) / iResolution.y * 2.) * fragCoord.y / iResolution.x));   
    //antialiasing
    vec3 tot = vec3(0.0);
    for(int m = 0; m < AA; m++) 
    for(int n = 0; n < AA; n++) {
            vec2 o = vec2(float(m), float(n)) / float(AA) - 0.5;
            vec2 p = (-iResolution.xy + 2.0 * (fragCoord + o)) / iResolution.y;
            vec3 rd = GetRayDir(p, ro, vec3(0, 0., 0), fl); //ray direction
            vec3 col = bg * bg; // background  
            
            //================isSphere==================
            float d = length(cross(ro, rd));
            if (d >= dist)
            {
                 tot += col;
                 continue;
            }
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
                    col = col1;
                    //bisect
                    //pos = (pos1 + pos)*0.5;
                    pos = getPoint(pos, pos1, val0, val1);
                    
                    vec3 nor = calcNormal(pos);
                    col = calccolor(col, col2, -rd, light, light2, nor);
                    break;
                }
                val0 = val1;
                pos = pos1;
            }
            tot += col;
        }
    tot = sqrt(tot) / float(AA);
    //antialiasing
    fragColor = vec4(tot, 1.0);
}
/////=====================================================================================
void main() {
    vec4 fragColor = vec4(0);
    mainImage(fragColor, gl_FragCoord.xy);
    gl_FragColor = fragColor;
}