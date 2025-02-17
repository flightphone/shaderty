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


float sdArc(vec2 p, float r, float a0, float a1)
{
    vec2 p0 = vec2(r*cos(a0), r*sin(a0));
    vec2 p1 = vec2(r*cos(a1), r*sin(a1));

    float a = mod(atan(p.y, p.x), TAU);
    if (a > a0 && a < a1)
        return abs(length(p) - r);   
    return min(length(p0-p), length(p1-p)); 
}

float sdArc2(vec2 p, float r, float a0, float a1)
{
    vec2 p0 = vec2(r*cos(a0), r*sin(a0));
    vec2 p1 = vec2(r*cos(a1), r*sin(a1));

    float a = mod(atan(p.y, p.x), TAU);
    if (!(a > a0 && a < a1))
        return abs(length(p) - r);   
    return min(length(p0-p), length(p1-p)); 
}


float sdSegment( in vec2 p, in vec2 a, in vec2 b )
{
    vec2 pa = p-a, ba = b-a;
    float h = clamp( dot(pa,ba)/dot(ba,ba), 0.0, 1.0 );
    return length( pa - ba*h );
}

// signed distance to a 2D triangle
float sdTriangleIsosceles( in vec2 p, in vec2 q )
{
    p.x = abs(p.x);
	vec2 a = p - q*clamp( dot(p,q)/dot(q,q), 0.0, 1.0 );
    vec2 b = p - q*vec2( clamp( p.x/q.x, 0.0, 1.0 ), 1.0 );
    float k = sign( q.y );
    float d = min(dot( a, a ),dot(b, b));
    float s = max( k*(p.x*q.y-p.y*q.x),k*(p.y-q.y)  );
	return sqrt(d)*sign(s);
}

vec3 curColor(float d1, vec3 col, vec3 col1, float dd)
{
    float s1 = smoothstep(0., -0.003, d1);
    float cs = dd*dd - (dd-abs(d1))*(dd-abs(d1));
    if (cs < 0.)
        cs = 1.;
    else
        cs = sqrt(cs)/dd;  
    vec3 colccs = col1*cs;      
    //return  (s1>0.) ? s1*col1*cs: col;
    return mix(col, colccs, vec3(s1));
}

vec3 vi4(vec2 p, float w1, float h1, float k)
{
     
    p.x *=k;
   
    float dd = 4./w1*k;
    vec3 col = vec3(0.78, 0.78, 0.53);
    vec3 col1 = vec3(0.95, 0.87, 0.6);
    vec3 col2 = vec3(0.35, 0.96, 0.97);
    float x1 = 14./w1*k;
    
    //=====================lines==========================
    float d1 = abs(p.x - x1) - dd;
    col = curColor(d1, col, col1, dd);

    float x2 = k - x1;
    d1 = abs(p.x - x2) - dd;
    col = curColor(d1, col, col1, dd);


    float y3 = 258./h1;
    d1 = abs(p.y - y3) - dd;
    col = curColor(d1, col, col1, dd);
    
    float y4 = 39./h1;
    d1 = sdSegment(p, vec2(x1, y4), vec2(x2, y4)) - dd;
    col = curColor(d1, col, col1, dd);

    float x5 = 0.5*k;
    d1 = sdSegment(p, vec2(x5, 0), vec2(x5, y3)) - dd;
    col = curColor(d1, col, col1, dd);

    float y6 = 188./h1;
    d1 = sdSegment(p, vec2(x1, y6), vec2(x2, y6)) - dd;
    col = curColor(d1, col, col1, dd);

    float y7 = 106./h1;
    d1 = sdSegment(p, vec2(x1, y7), vec2(0, y7)) - dd;
    col = curColor(d1, col, col1, dd);

    d1 = sdSegment(p, vec2(x2, y7), vec2(k, y7)) - dd;
    col = curColor(d1, col, col1, dd);
    //=====================lines==========================


    //===================romb=====================
    float r0 = (h1/2.- 88.)/h1*k;
   
    float dd2 = dd*0.5;
    vec2 pp = p - vec2(k/2., y4);
    float n = 4., ta = mod(atan(pp.y, pp.x), TAU), i = floor(ta/TAU*n);
    vec2 cnt = vec2(r0, r0);
    cnt.xy *= rot(i*PI/2.);
    cnt += vec2(k/2., y4);
    float a0 = PI + i*PI/2., a1 = PI + (i+1.)*PI/2.;
    if (a0 >= TAU)
    {
        a0 = mod(a0, TAU);
        a1 = mod(a1, TAU);
    }
    
    
    d1 = length(p - cnt) - r0;
    if (d1 > 0. && abs(pp.x) < r0 && abs(pp.y) < r0)
    {
        float cs = length(pp);
        float r02 = r0;
        cs = sqrt(r02*r02 - cs*cs)/r02;
        col = col2*cs*cs;
    }        
    
    
    d1 = sdArc(p - cnt, r0, a0, a1) - dd2;
    col = curColor(d1, col, col1, dd2);
    
    
    //===================romb=====================

    //====================floor 0==================
    float mir = 0.;
    if (p.x > k/2.)
    {
        mir = 1.;
        p.x = k - p.x;
    }
    
    float xr0 = 15./h1, yr0 = 150./h1, rr0 = 45./h1, xr1 = 75./h1, yr1 = 143./h1,
    rr1 = length(vec2(xr0, yr0) - vec2(xr1, yr1)) - rr0,
    xr2 = 60./h1, yr2 = 176./h1;
    //177 - 146
    a0 = asin((y6-yr0)/rr0);
    cnt = vec2(xr0, yr0) + rr0*vec2(cos(a0), sin(a0));
    float rr2 = length(vec2(xr2, yr2) - cnt);
    if (length(p - vec2(xr1, yr1)) - rr1 < 0.)
    {
        float cs = length(p - vec2(xr1, yr1));
        float rr = rr1*2.;
        cs = clamp(sqrt(rr*rr - cs*cs)/rr, 0.2, 1.);
        col = col2*cs;
        col = col2;
        
    }
    if (length(p - vec2(xr2, yr2)) - rr2 < 0. 
    && length(p - vec2(xr0, yr0)) - rr0 > 0.)
    {
        float cs = length(p - vec2(xr2, yr2));
        float rr = rr2*2.;
        cs = clamp(sqrt(rr*rr - cs*cs)/rr, 0.2, 1.);
        col = col2*cs;
        col = col2;
        
    }

    
    a1 = TAU - atan(yr0-yr1, xr1 - xr0);
    d1 = sdArc2(p - vec2(xr0, yr0), rr0, a0, a1) - dd2;
    
    a0 = PI - atan(yr0-yr1, xr1 - xr0);
    a1 = TAU;
    d1 = min(sdArc(p - vec2(xr1, yr1), rr1, a0, a1) - dd2, d1);
    
    a1 = PI - atan(cnt.y - yr2, xr2 - cnt.x);
    a0 = 0.;
    d1 = min(sdArc(p - vec2(xr2, yr2), rr2, a0, a1) - dd2, d1);
    
    col = curColor(d1, col, col1, dd2);

    if (mir == 1.)
        p.x = k - p.x;


    //=====================floor===================
    float yl = 163./h1, rl = 21./h1, al = PI/11.;
    vec2 ver = vec2(k/2., yl) + vec2(0., rl/cos(PI/2. - al));
    d1 = length(p - vec2(k/2., yl)) - rl;
    if (d1 < 0.)
    {
        float cs = length(p - vec2(k/2., yl));
        cs = sqrt(rl*rl - cs*cs)/rl;
        col = col2*cs;
        
    }
    vec2 tri = vec2(rl*cos(al), - ver.y +  rl*sin(al) + yl); // width, height    
    d1 = sdTriangleIsosceles(p- ver, 
    tri); 
    vec3 nor = vec3(-cos(al), 0., sin(al));
    if (d1 < 0.)
    {
        float rl2 = abs((p.y - ver.y)/tri.y*tri.x);
        float aa = PI/2. + asin((p.x - k/2.)/rl2);
        nor.xy *= rot(aa);
        float cs = dot(nor, vec3(0., -1., 0.));
        col = col2 * cs;
    }
   
    d1 = sdSegment(p, vec2(k/2., yl) + vec2(rl*cos(al), rl*sin(al)), ver) - dd2;
    d1 = min(sdSegment(p, vec2(k/2., yl) + vec2(-rl*cos(al), rl*sin(al)), ver) - dd2, d1);


    pp = p - vec2(k/2., yl);
    pp.xy *= rot(-al);
    d1 = min(sdArc(pp, rl, PI - 2.*al, 2.*PI) - dd2, d1);
    col = curColor(d1, col, col1, dd2);

    return col;
}

//https://www.xposz.shop/?ggcid=336928
void mainImage(out vec4 fragColor, in vec2 fragCoord) {
	//vec2 p = (-iResolution.xy + 2.0 * fragCoord) / iResolution.xy;
    vec2 p = fragCoord / iResolution.xy;
    float w1 = 184., h1 = 274., k = w1/h1;
    float n = iResolution.x/iResolution.y/k;
    p.x = fract(p.x*n);
    vec3 col = vi4(p, w1, h1, k);
   
	fragColor = vec4(col, 1.0);
}
/////=====================================================================================
void main() {
	vec4 fragColor = vec4(0);
	mainImage(fragColor, gl_FragCoord.xy);
	gl_FragColor = fragColor;
}