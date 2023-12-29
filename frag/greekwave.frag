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

vec3 greekwave(vec2 U)
{
    U = fract(U)-.5;                           // local tile coords in [-.5,.5]²
    
    float v, d = length(U),                    // distance to tile center
          t = 2.5*PI,                          // note the time delay with position
          a = t * max(0.,.5-d );               // angle ~time, and decrease with d
    U *= rot(a);   // rotate frame by angle(t,d)
    v = U.y;   
    vec3 O = vec3(smoothstep(-1.,1.,v*500. ));  
    return O;     
}

vec3 men1(vec2 p)
{
    float n = 4.;
    float w = 0.5;
    p = (p - 0.5)*n;
    vec2 a = floor(p);
    vec2 f = fract(p);
    vec3 res = vec3(1.0);
    float mx, mi;
    if (f.x < w)
    {
        if (a.x > 0.)
        {
            mx = a.x;
            mi = -a.x;
            
        }
        else
        {
            mx = -a.x + 1.;
            mi = a.x;
        }
        if ((a.y < mx && a.y >= mi) || (a.y == mx && f.y < w)   )
                res = vec3(0.);
    }

    if (f.y < w)
    {
        if (a.y >= 0.)
        {
            mx = a.y;
            mi = -a.y + 1.;
            
        }
        else
        {
            mx = -a.y;
            mi = a.y;
        }
        if (a.x < mx && a.x >= mi)
                res = vec3(0.);
    }
    
    return res;    
}
vec3 men3(vec2 p)
{
    p*=vec2(6.,7.);
    vec2 a = floor(p);
    vec3 res = vec3(1.0);
    if (a.y == 0. || a.y == 6.)
        res = vec3(0.0);
    if (a.x == 0. && a.y < 6. && a.y > 2.)    
        res = vec3(0.0);
    if (a.y == 2. && (a.x < 2. || a.x == 5.))    
        res = vec3(0.0);
    if (a.y == 4. && a.x >= 2. && a.x <=4.)    
        res = vec3(0.0);
    if (a.x == 3. && a.y < 4.)    
        res = vec3(0.0);    
    return 1. - res;    
}
vec3 men2(vec2 p)
{
    p*=vec2(8., 11.);
    vec2 a = floor(p);
    vec3 res = vec3(1.0);
    if (a.y == 0. || a.y == 10.)
        res = vec3(0.0);
    a.y = a.y - 2.;
    float mx, mi;
    
    if (mod(a.y, 2.) == 0.)
    {
        float h = abs(3.-a.y);
        float d = (a.y <=2.)? 0.:1.;
        mi = 5.-h-d;
        mx = 5.+h-2.*d;
        if (a.x <= mx && a.x >= mi) 
            res = vec3(0.);
    }    
    if (mod(a.x, 2.) == 0.)
    {
        mx = 0.;
        mi = 1.;
        if (a.x <= 4.)
        {
            float h = 3.-a.x;
            mi = clamp(2.-h, 0., 6.);
            mx = 3.+ h;

        }

        if (a.x == 6.0)
        {
            mi = 3.;
            mx = 6.;
        }

        if (a.y <= mx && a.y >= mi) 
            res = vec3(0.);
        
    }
    
    return res;
}

void mainImage(out vec4 fragColor, in vec2 fragCoord) 
{
    vec2 p = fragCoord/iResolution.xy;
    p*= vec2(6., 3.);
    float k = floor(p.y);
    p = fract(p);
    vec3 tot = vec3(1.);
    if (k == 0.)
        tot = greekwave(p);
    if (k == 1.)
        tot = men2(p);   
    if (k == 2.)
        tot = men3(p);           
         
    //men2(p);//
    fragColor = vec4(tot, 1.0);

    //https://indasil.club/1821-grecheskie-ornamenty-i-uzory.html
}
    
/////=====================================================================================
void main() {
    vec4 fragColor = vec4(0);
    mainImage(fragColor, gl_FragCoord.xy);
    gl_FragColor = fragColor;
}