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
braid rings
2D-SDF, tile, braid
simple braid pattern
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
    float shininess = 5.0;
    float specular = pow(max(dot(R1, rd), 0.), shininess);
    vec3 colccs = col1 * difu + 0.8 * specular;      
    //return  (s1>0.) ? s1*col1*cs: col;
    return mix(col, colccs, vec3(s1));
}
float bridge(vec2 p0, float n, float dd)
{
    /*
    float d = sdSegment(p0, vec2(0., 2. / n), vec2(1. / n, 2. / n))*n;
    d = min(sdSegment(p0, vec2(0., 2. / n), vec2(0. / n, 1.5 / n))*n, d);
    d = min(sdSegment(p0, vec2(1./n, 2. / n), vec2(1. / n, 1.5 / n))*n, d);
    */
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
    float n = 3., d = 10.;
    vec3 col = vec3(1., 1., 1.);
    vec3 col0 = vec3(0.2, 0.2, 1.);;
    p0.xy *= rot(-PI/4.);
    vec2 p = p0;
    p *= n;
    float numx = floor(p.x), numy = floor(p.y), num = floor(p.x) + floor(p.y);
    p = fract(p);

    //vec3 col1 = vec3(1., 0.5, 0.5);
    //vec3 col2 = vec3(0.5, 1., 0.5);
    float dd = 0.2, d2 = 10., d1 = 10.;
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
                d = d1;
            else
                d = d2;
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
        
    }

    col = curColorN(d - dd, col, col0, dd);
    //col = curColorN(d2, col, col0, dd );
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