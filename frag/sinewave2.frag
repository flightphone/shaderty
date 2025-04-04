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

float ticolor()
{
    return (1. + (0.5 - noise(iTime))*1.18);
}

vec3 vi4(vec2 p)
{
     
    p.x *=TAU;
    p.y -= 0.5;
   
    
    vec3 col = vec3(0.78, 0.78, 0.53);
    vec3 col1 = vec3(0.95, 0.87, 0.6);
    vec3 col2 = vec3(0.35, 0.96, 0.97);
    float dd = 0.05;
    float sign = (p.x < PI)?1.:-1.;
    float y = 0.3*cos(p.x)*sign, alf = PI/2. - (atan(-sin(p.x)*sign, 1.));
    float d1 = abs(p.y - y) - dd/sin(alf);
    //float d1 = abs(abs(p.x) - abs(x)) - dd;
    col = curColor(d1, col, col1, dd/sin(alf));

    y = -0.3*cos(p.x) * sign , alf = PI/2. - (atan(sin(p.x)*sign, 1.));
    d1 = abs(p.y - y) - dd/sin(alf);
    col = curColor(d1, col, col1, dd/sin(alf));

    

    return col;
}

//https://www.xposz.shop/?ggcid=336928
//https://ca.pinterest.com/pin/557109416405805140/
void mainImage(out vec4 fragColor, in vec2 fragCoord) {
	//vec2 p = (-iResolution.xy + 2.0 * fragCoord) / iResolution.xy;
    vec2 p = fragCoord / iResolution.xy;
    //p = fract(p*vec2(4., 2.));
    vec3 col = vi4(p);
   
	fragColor = vec4(col, 1.0);
}
/////=====================================================================================
void main() {
	vec4 fragColor = vec4(0);
	mainImage(fragColor, gl_FragCoord.xy);
	gl_FragColor = fragColor;
}