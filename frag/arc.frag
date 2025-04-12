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
float hash(float p) {
    p = fract(p * 0.011);
    p *= p + 7.5;
    p *= p + p;
    return fract(p);
}
//	<https://www.shadertoy.com/view/4dS3Wd>
//	By Morgan McGuire @morgan3d, http://graphicscodex.com

// This one has non-ideal tiling properties that I'm still tuning
float noise(float x) {
    float i = floor(x);
    float f = fract(x);
    float u = f * f * (3.0 - 2.0 * f);
    return mix(hash(i), hash(i + 1.0), u);
}


float sdSegment(in vec2 p, in vec2 a, in vec2 b) {
    vec2 pa = p - a, ba = b - a;
    float h = clamp(dot(pa, ba) / dot(ba, ba), 0.0, 1.0);
    return length(pa - ba * h);
}


vec3 curColorN(float d1, vec3 col, vec3 col1, float dd)
{
    float s1 = smoothstep(0., -5./iResolution.y, d1);
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
    vec3 rd = vec3(0., 0., -1.);
    float difu = dot(norm , light);   
    vec3 R1 = reflect (light, norm);
    float shininess=5.0;
    float specular    =  pow(max(dot(R1, rd), 0.), shininess);
    vec3 colccs = col1*difu + 0.8*specular;      
    
    return mix(col, colccs, vec3(s1));
}

vec2 drog(vec2 p) {
    float x0 = 0.04, w0 = 0.3, w1 = 0.31, y1 = 0.14, k = 1.0;
    float d = sdSegment(p, vec2(x0, -x0), vec2(x0 + w0, -x0));

    float r0 = 0.025, alf0 = atan(y1, w0 - w1) + TAU, v = (length(vec2(w0 - w1, y1)) - r0) / alf0, dlt = v * TAU;
    vec2 pp = p - vec2(x0 + w1, -x0 - y1);
    float alf1 = mod(atan(pp.y, pp.x), TAU), alf2 = min(alf0, alf1 + TAU);
    vec2 p0 = vec2(r0, 0.), p1 = (r0 + v * alf1) * vec2(cos(alf1), sin(alf1)), p2 = (r0 + v * alf2) * vec2(cos(alf2), sin(alf2));
    float d2 = length(pp - p0);
    float eps = 0.4;
    k = eps + (1.-eps)*alf1/alf0;
    if (length(pp - p2) < length(pp - p1))
        k = eps + (1.-eps)*alf2/alf0;
    if (d2 <  length(pp - p2) && d2 < length(pp - p1))
        k = eps;

    d2 = min(d2, length(pp - p1));
    d2 = min(d2, length(pp - p2));
    if (d < d2)
        k = 1.;
    d = min(d, d2);
    return vec2(d, k);
}

vec3 ya(vec2 p, vec3 col) {
    //p -= 0.5;
    //vec3 col = vec3(1.);
    vec2 dd = drog(p);
    float d = dd[0], k = dd[1];
    vec2 dd2 = drog(vec2(-p.y, -p.x));
    float d2 = dd2[0], k2 = dd2[1];
    if (d2 < d)
    {
        d = d2;
        k = k2;
    }
    float h = 0.04;
    //float eps = 15./iResolution.y, h = 0.015;
    //float s = smoothstep(eps, 0., d-h*sin(k*PI/2.));
    //col = mix(col, vec3(0.), s);
    
    vec3 col1 = vec3(1., 0.2, 0.2);
    float dh = h*sin(k*PI/2.);
    col = curColorN(d-dh, col, col1, dh); 
    return col;
}

vec3 ya2(vec2 p) {
    p -= 0.5;
    p*=1.2;
    //vec3 col = draw_cyrcle(p, 0.49, 0.01, vec3(1.), vec3(0.));
    vec3 col = vec3(1.);
    float a = mod(atan(p.y, p.x), TAU);
    float n = 4.; 
    float nn = floor(a/TAU*n);
    p.xy *= rot(-nn*TAU/n);
    p.xy *= rot(-PI/2.);

    col = ya(p, col);
    return col;
}


//https://www.xposz.shop/?ggcid=336928
//https://ca.pinterest.com/pin/557109416405805140/
void mainImage(out vec4 fragColor, in vec2 fragCoord) {
	//vec2 p = (-iResolution.xy + fragCoord) / iResolution.y /2.;
    vec2 p = fragCoord / iResolution.y;
    p.x -= (iResolution.x-iResolution.y)/iResolution.y/2.;
    
    //vec2 p = fragCoord / iResolution.xy;
    //float k = 1.0;
    //float n = iResolution.x / iResolution.y / k;
    //p.x = fract(p.x * n);

    vec3 col = ya2(p);

    fragColor = vec4(col, 1.0);
}
/////=====================================================================================
void main() {
    vec4 fragColor = vec4(0);
    mainImage(fragColor, gl_FragCoord.xy);
    gl_FragColor = fragColor;
}