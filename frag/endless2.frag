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
endless knot
2D_SDF,knot,endless, tibet, buddhism 
[url]https://en.wikipedia.org/wiki/Endless_knot[/url]
*/
#define PI  3.14159265359
#define TAU 6.28318530718
#define rot(f) mat2(cos(f), -sin(f), sin(f), cos(f))

float sdSegment(in vec2 p, in vec2 a, in vec2 b) {
    vec2 pa = p - a, ba = b - a;
    float h = clamp(dot(pa, ba) / dot(ba, ba), 0.0, 1.0);
    return length(pa - ba * h);
}


vec3 curColorN(float d1, vec3 col, vec3 col1, float dd) {
    float s1 = smoothstep(0., -15. / iResolution.y, d1);
    float cs = dd * dd - (dd - abs(d1)) * (dd - abs(d1));
    if(cs < 0.) {
        cs = 1.;
        d1 = 0.;
    } else {
        cs = sqrt(cs) / dd;
        d1 = d1 / dd;
    }
    vec3 norm = vec3(d1, 0., cs);
    vec3 light = vec3(0., 0., 1.);
    vec3 rd = vec3(0., 0., -1.);
    float difu = dot(norm, light);
    vec3 R1 = reflect(light, norm);
    float shininess = 15.0;
    float specular = pow(max(dot(R1, rd), 0.), shininess);
    vec3 colccs = col1 * difu + 0.8 * specular;      
    return mix(col, colccs, vec3(s1));
}



vec3 draw_cyrcle(vec2 p, float r, float dd, vec3 col, vec3 col1)
{
    float d1 = abs(length(p) - r) - dd;
    return curColorN(d1, col, col1, dd);   
    
}

vec3 curColor(float d1, vec3 col, vec3 col1)
{
    float s1 = smoothstep(0., -0.003, d1);
    return mix(col, col1, vec3(s1));
}

vec3 draw_disc(vec2 p, float r, vec3 col, vec3 col1)
{
    float d1 = length(p) - r;
    col1  =  curColor(d1, col, col1);
    return col1;
}


float bridge(vec2 p0, float n, float dd)
{
    vec2 c = vec2(0.5/n, (1. + dd) / n);
    float d = abs(length(p0 - c) - 0.5/n)*n;
    return d;
}

float bridge2(vec2 p0, float n, float dd)
{
    float d = sdSegment(p0, vec2(-1./n, 1. / n), vec2(-1. / n, 1.5 / n))*n;
    d = min(sdSegment(p0, vec2(-1./n, 1. / n), vec2(-1.5 / n, 1. / n))*n, d);
    
    vec2 c = p0 - vec2(-1.5/n, 1.5/n);
    float a = mod(atan(c.y, c.x), TAU);
    if (a <= 1.5*PI)
        d = min(d, abs(length(c) - 0.5/n)*n);
    
    return d;

}

vec3 endless(vec2 p0) {
    float n = 3.1, d = 10.;
    vec3 b1 = vec3(0.23529411764705882, 0.4235294117647059, 0.7725490196078432), b2 = vec3(0.3686274509803922, 0.5725490196078431, 0.8941176470588236);
    vec3 bg = mix(b2, b1, p0.y);  
    vec3 col = bg;
    vec3 col0 = vec3(222./255., 208./255., 159./255.);
    vec3 col1 = vec3(0.3, 0.3, 1.);;
    p0.xy *= rot(-PI/4.);
    vec2 p = p0;
    p *= n;
    float numx = floor(p.x), numy = floor(p.y), num = floor(p.x) + floor(p.y);
    p = fract(p);
    float dd = 0.2, d2 = 10., d1 = 10.;
    col = draw_disc(p0, 0.99, col, col0);  
    col = draw_cyrcle(p0, 0.98, 0.02, col, col1);
    if(abs(p0.x) < (1. + dd) / n && abs(p0.y) < (1. + dd) / n) {
        d1 = min(p.x, 1. - p.x);
        d2 = min(p.y, 1. - p.y);
        d = min(d1, d2);
        if(d1 < dd && d2 < dd) {
            if(p.x < dd && p.y < dd) {
                d2 = p.y;
                d1 = p.x;
            }
            if(p.x < dd && 1. - p.y < dd) {
                d2 = 1. - p.y;
                d1 = p.x;
                num += 1.;
            }
            if(1. - p.x < dd && 1. - p.y < dd) {
                d2 = 1. - p.y;
                d1 = 1. - p.x;
                num += 2.;
            }
            if(1. - p.x < dd && p.y < dd) {
                d2 = p.y;
                d1 = 1. - p.x;
                num += 1.;
            }
            if(mod(num, 2.0) == 0.)
            {
                //d = d1;
                col = curColorN(d1 - dd, col, col1, dd);
                col = curColorN(d2 - dd, col, col1, dd);
            }
            else
            {
                //d = d2;
                col = curColorN(d2 - dd, col, col1, dd);
                col = curColorN(d1 - dd, col, col1, dd);
            }
        }
        else
        {
            col = curColorN(d - dd, col, col1, dd);
        }
    } else {
        if (p0.y > 0.)
            d = min(d, bridge(p0, n, dd));
        if (p0.y < 0.)    
            d = min(d, bridge(vec2(-p0.x, -p0.y), n, dd));
        if (p0.x > 0.)
            d = min(d, bridge(vec2(p0.y, p0.x), n, dd));    
        if (p0.x < 0.)
            d = min(d, bridge(vec2(-p0.y, -p0.x), n, dd));  

        
        if (p0.x < 0. && p0.y > 0.)          
        {
            d = min(bridge2(p0, n, dd), d);
        }
        if (p0.x > 0. && p0.y < 0.)          
        {
            d = min(bridge2(-p0, n, dd), d);
        }
        col = curColorN(d - dd, col, col1, dd);
        
    }
    
    return col;
}

void mainImage(out vec4 fragColor, in vec2 fragCoord) {

    vec2 p = (2.0 * fragCoord - iResolution.xy) / iResolution.y;
    vec3 col = endless(p);
    fragColor = vec4(col, 1.0);
}
/////=====================================================================================
void main() {
    vec4 fragColor = vec4(0);
    mainImage(fragColor, gl_FragCoord.xy);
    gl_FragColor = fragColor;
}