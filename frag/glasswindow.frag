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
/*
stained glass windows
2D-SDF, tile, noise
stained glass windows
*/
#define PI  3.14159265359
#define TAU 6.28318530718
#define rot(f) mat2(cos(f), -sin(f), sin(f), cos(f))

// Precision-adjusted variations of https://www.shadertoy.com/view/4djSRW
float hash(float p) { p = fract(p * 0.011); p *= p + 7.5; p *= p + p; return fract(p); }
//	<https://www.shadertoy.com/view/4dS3Wd>
//	By Morgan McGuire @morgan3d, http://graphicscodex.com

// This one has non-ideal tiling properties that I'm still tuning
float noise(float x) {
	float i = floor(x);
	float f = fract(x);
	float u = f * f * (3.0 - 2.0 * f);
	return mix(hash(i), hash(i + 1.0), u);
}

float sdArc(vec2 p, float r, float a0, float a1)
{
    vec2 p0 = vec2(r*cos(a0), r*sin(a0));
    vec2 p1 = vec2(r*cos(a1), r*sin(a1));

    float a = mod(atan(p.y, p.x), TAU);
    if (a > a0 && a < a1)
        return abs(length(p) - r);   
    //return min(length(p0-p), length(p1-p)); 
    return 100.;
}

float sdArc3(vec2 p, float r, float a0, float a1)
{
    vec2 p0 = vec2(r*cos(a0), r*sin(a0));
    vec2 p1 = vec2(r*cos(a1), r*sin(a1));

    float a = mod(atan(p.y, p.x), TAU);
    if (a > a0 && a < a1)
        return abs(length(p) - r);   
    return min(length(p0-p), length(p1-p)); 
    //return 100.;
}



float sdArcElips(vec2 p, float a, float b, float t0, float t1)
{
    float t = mod(atan(p.y/b, p.x/a), TAU);
    if (t > t0 && t < t1)
        return abs(p.x*p.x/a/a + p.y*p.y/b/b - 1.);   
    return +100.;
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

vec3 curColor(float d1, vec3 col, vec3 col1)
{
    float s1 = smoothstep(0., -0.003, d1);
    return mix(col, col1, vec3(s1));
}

vec3 curColorN(float d1, vec3 col, vec3 col1, float dd, vec3 rd)
{
    float s1 = smoothstep(0., -0.003, d1);
    float cs = dd*dd - (dd-abs(d1))*(dd-abs(d1));
    if (cs < 0.)
    {
        cs = 1.;
        d1 = 0.;
    }
    else
    {
        cs = sqrt(cs)/dd;  
        d1 = d1/dd;
    }
    vec3 norm = vec3(d1, 0., cs);
    vec3 light = vec3(0., 0., 1.);
    //vec3 rd = vec3(0., 0., -1.);
    float difu = dot(norm , light);   
    vec3 R1 = reflect (light, norm);
    float shininess=10.0;
    float specular    =  pow(max(dot(R1, rd), 0.), shininess);
    vec3 colccs = col1*difu + 0.8*specular;      
    //return  (s1>0.) ? s1*col1*cs: col;
    return mix(col, colccs, vec3(s1));
}

float ticolor()
{
    return (1. + (0.5 - noise(iTime))*1.18);
}

vec3 draw_cyrcle(vec2 p, float r, float dd, vec3 col, vec3 col1)
{
    float d1 = abs(length(p) - r) - dd;
    return curColor(d1, col, col1);   
    
}

float prism(vec2 p)
{
    p-=0.5;
    float t = max(abs(p.x), abs(p.y));
    return 1.0 - t;
}

float glass_h(vec2 p) {
    
    
    float n = 25., k = 1./n*0.1;
	p = fract(p*n);
    float t = prism(p)*k;
	return t;
    
 
}

vec3 glass_norm(vec2 p)
{
	float h = 0.00001, dx = glass_h(p + vec2(h, 0)) - glass_h(p - vec2(h, 0)),
	dy = glass_h(p + vec2(0, h)) - glass_h(p - vec2(0, h)), dz = -2.*h;
	return -normalize(vec3(dx, dy, dz));
}


vec3 draw_disc(vec2 p, float r, vec3 col, vec3 col1)
{
    float d1 = length(p) - r;
    col1  =  curColor(d1, col, col1);
    return col1;
}


vec3 draw_disc_glass(vec2 p, float r, vec3 col, vec3 col1)
{
    float d1 = length(p) - r;
       
    vec3 norm = glass_norm(p);
    vec3 light = normalize(vec3(sin(iTime/3.), cos(iTime/3.), 1.));
	//vec3 light = normalize(vec3(1., 1., 1.));
    vec3 rd = vec3(0., 0., -1.);
    float difu = dot(norm , light);   
    vec3 R1 = reflect (light, norm);
    float shininess=5.0;
    float specular    =  pow(max(dot(R1, rd), 0.), shininess);
    col1 = col1*difu + 0.8*specular;      
    col1 =  pow(col1, vec3(1.));
    col1  =  curColor(d1, col, col1);
    return col1;
}

vec3 vi4(vec2 p, float k)
{
    p.x *=k;
    float h = 737., r1 = 135./h, dd = 35./h/2.;

    //vec3 col = vec3(0.3, 0.3, 0.3);
    vec3 col0 = vec3(0.3);
    vec3 col1 = vec3(0.05, 0.06, 0.09);
    vec3 col2 = vec3(0.5, 0.21, 0.1);
    vec3 col3 = vec3(0.92, 0.81, 0.51);
    vec3 col4 = vec3(0.64, 0.49, 0.37);
    vec3 col5 = vec3(0.91, 0.24, 0.054);
    vec3 col6 = vec3(0.76, 0.9, 0.98);
    vec3 col7 = vec3(0.57, 0.52, 0.48);
    vec3 col8 = vec3(0.5372, 0.4980, 0.2980);
    vec3 col9 = vec3(0.7215, 0.1724, 0.0352);
    vec3 col10 = vec3(0.9921, 0.7098, 0.1333);

    vec3 col =col4*col4;//texture(iChannel0, p).rgb;
    vec2 pp = p - 0.5;
    float rr0 = 0.5, rr1 = rr0-dd, dd2 = dd*0.65, rr2 = rr1 - dd - dd2,
    dd3 = 0.75*dd2, rr3 = rr2 - dd2 - dd3, rr4 = rr3-dd3;
    col = draw_disc(pp, rr0, col, col1);
    
    col = draw_cyrcle(pp, rr2, dd2, col, col3);
    //================segments===================
    float a = mod(atan(pp.y, pp.x), TAU); 
    float ns1 = 12., dt1 = TAU/ns1, 
    as0 = floor(a/dt1)*dt1, as1 = as0 + dt1;
    float d1 = sdSegment(pp, (rr2-dd2)*vec2(cos(as0), sin(as0)), (rr2+dd2)*vec2(cos(as0), sin(as0))) - 0.003;
    d1 = min(sdSegment(pp, (rr2-dd2)*vec2(cos(as1), sin(as1)), (rr2+dd2)*vec2(cos(as1), sin(as1)))- 0.003, d1);
    col = curColor(d1, col, col1);


    //===========================================
    col = draw_cyrcle(pp, rr3, dd3, col, col4);
    col = draw_disc_glass(pp, rr4, col, col5);
    
    
    //===========================lett cyrcle
    d1 = abs(length(pp) - (0.5 - dd)) - dd;
    float n = 4., aa = a + PI/n;
    float nn = floor(aa/TAU*n);
    pp.xy *= rot(-nn*TAU/n);
    pp.xy *= rot(PI);
    pp = pp + vec2(0.5) - vec2(0.5 - r1*sqrt(2.), 0.5);
    
    
    col = draw_disc(pp, r1+dd, col, col2);
    d1 = min(sdArc3(pp, r1, PI/4., TAU - PI/4.) - dd, d1);
    vec3 rd = normalize(vec3(0., 0., -1.));
    col = curColorN(d1, col, col0, dd, rd);
    float r2 = r1 - dd - dd2, r3 = r2 - dd2;
    col = draw_disc(pp, r2+0.018, col, col1);
    col = draw_cyrcle(pp, r2, dd2, col, col6);
    col = draw_disc_glass(pp, r3, col, col7);
    //======================segments======================
    ns1 = 10.; 
    float a10 = mod(atan(pp.y, pp.x), TAU);
    as0 = floor(a10/dt1)*dt1;
    as1 = as0 + dt1;
    
    d1 = sdSegment(pp, (r2-dd2)*vec2(cos(as0), sin(as0)), (r2+dd2)*vec2(cos(as0), sin(as0))) - 0.0025;
    d1 = min(sdSegment(pp, (r2-dd2)*vec2(cos(as1), sin(as1)), (r2+dd2)*vec2(cos(as1), sin(as1)))- 0.0025, d1);
    col = curColor(d1, col, col1);


    /*=======================hrizo====================*/
    float nh = 12., ah = mod(atan(pp.y, pp.x), TAU), rh = r3*0.55, r6h = rh*sin(PI/nh);
    col = draw_disc(pp, rh, col, col3);
    ah += PI/nh;
    float ih = floor(ah/TAU*nh);
    vec2 ppp = pp;
    ppp.xy *= rot(-ih*TAU/nh);
    ppp.xy *= rot(PI);
    vec2 cnl = vec2(-rh*cos(PI/nh), 0.);
    float dh = 0.003;
    col = draw_disc(ppp - cnl, r6h, col, col3);
    d1 = sdArc(ppp - cnl, r6h, PI/2., TAU - PI/2.) - 0.003;
    col = curColor(d1, col, col1);
    d1 = sdSegment(ppp, cnl + vec2(0., r6h), vec2(0., 0.)) - dh;
    d1 = min(sdSegment(ppp, cnl - vec2(0., r6h), vec2(0., 0.)) - dh, d1);
    col = curColor(d1, col, col1);


   /*============center================================*/
    r2 = r2, r3 = r2 - dd2;
    pp = p - 0.5;
    col = draw_disc(pp, r2+0.018, col, col1);
    col = draw_cyrcle(pp, r2, dd2, col, col6);
    col = draw_disc(pp, r3, col, col7);
    float dd4 = 33./h/2., r4 = r3 - dd4, r5 = r4 - dd4;
    col = draw_cyrcle(pp, r4, dd4, col, col8);
    col = draw_disc(pp, r5, col, col1);
    col = draw_disc_glass(pp, r5-0.0018, col, col9);

    float r6 = r5/2.- 0.002;

    //================segments===================
    ns1 = 8.;
    dt1 = TAU/ns1; 
    as0 = floor(a/dt1)*dt1;
    as1 = as0 + dt1;
    d1 = sdSegment(pp, (r4-dd4)*vec2(cos(as0), sin(as0)), (r4+dd4)*vec2(cos(as0), sin(as0))) - 0.0025;
    d1 = min(sdSegment(pp, (r4-dd4)*vec2(cos(as1), sin(as1)), (r4+dd4)*vec2(cos(as1), sin(as1)))- 0.0025, d1);
    col = curColor(d1, col, col1);

    
    as0 = floor((a + PI/ns1)/dt1)*dt1 - PI/ns1;
    as1 = as0 + dt1;
    
    d1 = sdSegment(pp, (r2-dd2)*vec2(cos(as0), sin(as0)), (r2+dd2)*vec2(cos(as0), sin(as0))) - 0.0025;
    d1 = min(sdSegment(pp, (r2-dd2)*vec2(cos(as1), sin(as1)), (r2+dd2)*vec2(cos(as1), sin(as1)))- 0.0025, d1);
    col = curColor(d1, col, col1);


    //===========================================

    //pp = p-0.5;
    pp.xy *= rot(-nn*TAU/n);
    pp.xy *= rot(PI);
    pp = pp + vec2(0.5) - vec2(0.5 - r6, 0.5);
    col = draw_disc(pp, r6, col, col10);
    d1 = sdArc(pp, r6, PI/2., TAU - PI/2.) - 0.003;
    col = curColor(d1, col, col1);
    d1 = sdSegment(pp, vec2(0., r6), vec2(r6, 0.)) - 0.003;
    d1 = min(sdSegment(pp, vec2(0., -r6), vec2(r6, 0.)) - 0.003, d1);
    col = curColor(d1, col, col1);

    
    /*============romb================*/
    pp = p - vec2(0.5);
    float ta = mod(atan(pp.y, pp.x), TAU), i = floor(ta/TAU*n);
    vec2 cnt = vec2(r5, r5);
    cnt.xy *= rot(i*PI/2.);
    float a0 = PI + i*PI/2., a1 = PI + (i+1.)*PI/2.;
    if (a0 >= TAU)
    {
        a0 = mod(a0, TAU);
        a1 = mod(a1, TAU);
    }
    
    
    d1 = length(pp - cnt) - r5;
    if (d1 > 0. && abs(pp.x) < r5 && abs(pp.y) < r5)
    {
        col = vec3(1.);
    }        
    d1 = sdArc(pp - cnt, r5, a0, a1) - 0.003;
    col = curColor(d1, col, col1);
    /*============romb================*/

    float r7 = r5*sqrt(2.) - r5, dd5 = 5./h, r8 = r7-dd5 - 0.002;
    pp = p - 0.5;
    col = draw_disc(pp, r7, col, col1);

    col = draw_cyrcle(pp, r8, dd5, col, col9);
    col = draw_disc(pp, r8-dd5, col, col10);


    return col;
}

//https://www.xposz.shop/?ggcid=336928
//https://ca.pinterest.com/pin/557109416405805140/
void mainImage(out vec4 fragColor, in vec2 fragCoord) {
	//vec2 p = (-iResolution.xy + 2.0 * fragCoord) / iResolution.xy;
    vec2 p = fragCoord / iResolution.xy;
    float k = 1.0;
    float n = iResolution.x/iResolution.y/k;
    p.x = fract(p.x*n);
    vec3 col = vi4(p, k);
   
	fragColor = vec4(col, 1.0);
}
/////=====================================================================================
void main() {
	vec4 fragColor = vec4(0);
	mainImage(fragColor, gl_FragCoord.xy);
	gl_FragColor = fragColor;
}