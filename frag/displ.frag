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



float sdBox( in vec2 p, in vec2 b )
{
    vec2 d = abs(p)-b;
    return length(max(d,0.0)) + min(max(d.x,d.y),0.0);

}
float sdBox2(vec2 p, float x0, float y0, float x1, float y1)
{
    return sdBox(p - vec2((x1+x0)/2., (y1 + y0)/2.), vec2((x1-x0)/2., (y1-y0)/2.));
}

float sdArc(vec2 p, float r, float a0, float a1)
{
    vec2 p0 = vec2(r*cos(a0), r*sin(a0));
    vec2 p1 = vec2(r*cos(a1), r*sin(a1));

    float a = mod(atan(p.y, p.x), TAU);
    /*
    float res = length(p0-p)*step(a0, a) + length(p1-p)*step(a, a1)
    + abs(length(p) - r)*step(a, a0)*step(a1, a);
    */
    if (a > a0 && a < a1)
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

vec3 grid0(vec2 p, float k)
{
    p.x *=k;
    vec2 cn = vec2(k/2., 0.5);
    vec3 col = vec3(1.);
    float d = sdBox(p-cn, vec2(0.45*k, 0.45));
    vec3 col1 = vec3(smoothstep(0., 0.001, d));
    col -= col1;

    float d2 = 0.2 - length(p-cn);
    vec3 col2 = vec3(smoothstep(0., 0.001, d2));
    col -= col2;

    return col;
}

vec3 grid2(vec2 p, float k)
{
    p.x *=k;
    float h1 = 167., w1= 68.; 
    float y0 = 101./h1;
    float x0 = 5./w1*k;
    float x1 = (w1 - 5.)/w1*k;
    float y1 = 130./h1;

    float dx = 5./w1*k;
    float dy = 5./h1;

    float bound = 0.005;
    float bound2 = 0.01;


    vec3 col = vec3(1.);

    float d0 = sdBox2(p, x0, y0, x1, y1);
    vec3 col0 = vec3(smoothstep(-bound, -bound+0.001, d0)*smoothstep(bound, bound-0.001, d0));
    col -= col0;

    x0+=dx;
    y0+=dy;
    x1-=dx;
    y1-=dy;
    float d1 = sdBox2(p, x0, y0, x1, y1);
    vec3 col1 = vec3(smoothstep(bound2, bound2-0.001, d1));
    col -= col1;

    float y2 = 5./h1;
    float x2 = 5./w1*k;
    float x3 = (w1 - 5.)/w1*k;
    float y3 = 91./h1;
    float d2 = sdBox2(p, x2, y2, x3, y3);
    vec3 col2 = vec3(smoothstep(-bound, -bound+0.001, d2)*smoothstep(bound, bound-0.001, d2));
    col -= col2;

    x2+=dx;
    y2+=dy;
    x3-=dx;
    y3-=dy;
    float d3 = sdBox2(p, x2, y2, x3, y3);
    vec3 col3 = vec3(smoothstep(bound2, bound2-0.001, d3));
    col -= col3;

    float r = 5./h1;
    float yr = y3 - 6./h1;
    float xr = x2 + 6./w1*k;
    float dr0 = r - length(p - vec2(xr, yr));
    vec3 colr = vec3(smoothstep(-0.001, 0., dr0));
    col += colr;


    xr = k - xr;
    dr0 = r - length(p - vec2(xr, yr));
    colr = vec3(smoothstep(-0.001, 0., dr0));
    col += colr;

    float h6 = 153./h1;
    float d6 = abs(p.y - h6);
    float bound6 = 5./h1;
    vec3 col6 = vec3(smoothstep(bound6 + 0.001, bound6, d6));
    col -= col6;

    float dx2 = (dx - bound)/4.;
    float y7 = 51./h1;
    float r7 = 14./h1;
    float x7 = x2 + r7 - 2.5*dx2;
    float d7 = sdArc(p - vec2(x7, y7), r7, PI/2., PI);
    vec3 col7 = vec3(smoothstep(dx2, dx2-0.001, d7));
    col += col7;

    float x8 = k - x7;
    float d8 = sdArc(p - vec2(x8, y7), r7, 0., PI/2.);
    vec3 col8 = vec3(smoothstep(dx2, dx2-0.001, d8));
    col += col8;

    float x9 = (x7 + x8)/2., r9 = (x8-x9 - 2.*dx2)/2./k,
    d9 = sdArc(p - vec2(x9, y7 + r7), r9, 0., PI);
    vec3 col9 = vec3(smoothstep(dx2, dx2-0.001, d9));
    col += col9;



    return col;
}

vec3 ruf(vec2 p) {
    p.x *= TAU*25.;
    vec3 col = vec3(clamp(sin(p.x), 0., 1.));
    return col;
}

vec3 testArc(vec2 p)
{
    vec3 col = vec3(0);
    float dx = 0.003;
    float d = sdArc(p - vec2(0.5, 0.5), 0.4, 0., PI/2.) - dx;
    vec3 col8 = vec3(smoothstep(0., -0.001, d));
    col += col8;
    return col;

}
vec3 curColor(float d1, vec3 col, vec3 col1, float dd)
{
    float s1 = smoothstep(0., -0.001, d1);
    float cs = dd*dd - (dd-abs(d1))*(dd-abs(d1));
    if (cs < 0.)
        cs = 1.;
    else
        cs = sqrt(cs)/dd;    
    return  (s1>0.) ? s1*col1*cs: col;
}

vec3 vi4(vec2 p, float k)
{
    
    p.x *=k;
    float h1 = 274., w1= 184.; 
    float dd = 3./w1*k;
    vec3 col = vec3(0.78, 0.78, 0.53);
    vec3 col1 = vec3(0.95, 0.87, 0.6);
    vec3 col2 = vec3(0.35, 0.96, 0.97);
    float x1 = 14./w1*k;
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

    float y6 = 177./h1;
    d1 = sdSegment(p, vec2(x1, y6), vec2(x2, y6)) - dd;
    col = curColor(d1, col, col1, dd);

    float y7 = 106./h1;
    d1 = sdSegment(p, vec2(x1, y7), vec2(0, y7)) - dd;
    col = curColor(d1, col, col1, dd);

    d1 = sdSegment(p, vec2(x2, y7), vec2(k, y7)) - dd;
    col = curColor(d1, col, col1, dd);


    float r0 = (h1/2.- 88.)/h1*k;
    float xr0 = k/2. - r0, yr0 = y4 - r0, xr1 = xr0, yr1 = y4 + r0, xr2 = k/2. + r0, 
    yr2 = yr1, xr3 = xr2, yr3 = y4 - r0;

    float dr0 = length(p-vec2(xr0, yr0)) - r0,
    dr1 = length(p-vec2(xr1, yr1)) - r0,
    dr2 = length(p-vec2(xr2, yr2)) - r0,
    dr3 = length(p-vec2(xr3, yr3)) - r0;
    if (dr0>0. && dr1 > 0. && dr2 > 0. && dr3 > 0. 
    && abs(p.x - k/2.) < r0 && abs(p.y - y4) < r0 )
    {
        float cs = length(p - vec2(k/2., y4));
        float r02 = r0/1.5;
        cs = sqrt(r02*r02 - cs*cs)/r02;
        col = col2*cs;
    }

    

    d1 = sdArc(p - vec2(xr0, yr0), r0, 0., PI/2.) - dd;
    col = curColor(d1, col, col1, dd);

    d1 = sdArc(p - vec2(xr1, yr1), r0, PI*1.5, TAU) - dd;
    col = curColor(d1, col, col1, dd);

    d1 = sdArc(p - vec2(xr2, yr2), r0, PI, PI*1.5) - dd;
    col = curColor(d1, col, col1, dd);

    d1 = sdArc(p - vec2(xr3, yr3), r0, PI/2., PI) - dd;
    col = curColor(d1, col, col1, dd);

    float yl = 163./h1, rl = 20./h1, al = PI/8.;
    vec2 ver = vec2(k/2., yl) + vec2(0., rl/cos(PI/2. - al));
    dr0 = length(p - vec2(k/2., yl)) - rl;
    if (dr0 < 0.)
    {
        float cs = length(p - vec2(k/2., yl));
        cs = sqrt(rl*rl - cs*cs)/rl;
        col = col2*cs;
        
    }
    vec2 tri = vec2(rl*cos(al), - ver.y +  rl*sin(al) + yl); // width, height    
    dr0 = sdTriangleIsosceles(p- ver, 
    tri); 
    vec3 nor = vec3(-cos(al), 0., sin(al));
    if (dr0 < 0.)
    {
        float rl2 = abs((p.y - ver.y)/tri.y*tri.x);
        float aa = PI/2. + asin((p.x - k/2.)/rl2);
        nor.xy *= rot(aa);
        float cs = dot(nor, vec3(0., -1., 0.));
        col = col2 * cs;
    }


    d1 = sdArc(p - vec2(k/2., yl), rl, PI - al, 2.*PI) - dd;
    col = curColor(d1, col, col1, dd);
    d1 = sdArc(p - vec2(k/2., yl), rl, 0., al) - dd;
    col = curColor(d1, col, col1, dd);
    
    d1 = sdSegment(p, vec2(k/2., yl) + vec2(rl*cos(al), rl*sin(al)), ver) - dd;
    col = curColor(d1, col, col1, dd);
    d1 = sdSegment(p, vec2(k/2., yl) + vec2(-rl*cos(al), rl*sin(al)), ver) - dd;
    col = curColor(d1, col, col1, dd);

    return col;
}
//https://www.xposz.shop/?ggcid=336928
void mainImage(out vec4 fragColor, in vec2 fragCoord) {
	//vec2 p = (-iResolution.xy + 2.0 * fragCoord) / iResolution.xy;
    vec2 p = fragCoord / iResolution.xy;
    
    //vec3 col = ruf(p);
        
    //float n = iResolution.x/iResolution.y;
    //p.x = fract(p.x*n);
    //vec3 col = testArc(p);

    float k = 184./274.;
    float n = iResolution.x/iResolution.y/k;
    p.x = fract(p.x*n);
    vec3 col = vi4(p, k);
    
    //float k = 68./167.;
    //float n = iResolution.x/iResolution.y/k;
    //p.x = fract(p.x*n);
    //vec3 col = grid2(p, k);

    //p.x = fract(p.x*iResolution.x/iResolution.y);
    //vec3 col = grid0(p, k);
    
	fragColor = vec4(col, 1.0);
}
/////=====================================================================================
void main() {
	vec4 fragColor = vec4(0);
	mainImage(fragColor, gl_FragCoord.xy);
	gl_FragColor = fragColor;
}