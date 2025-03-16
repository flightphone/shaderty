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

float sdArc(vec2 p, float r, float a0, float a1) {
    vec2 p0 = vec2(r * cos(a0), r * sin(a0));
    vec2 p1 = vec2(r * cos(a1), r * sin(a1));

    float a = mod(atan(p.y, p.x), TAU);
    if(a > a0 && a < a1)
        return abs(length(p) - r);   
    //return min(length(p0-p), length(p1-p)); 
    return 100.;
}

float sdArc3(vec2 p, float r, float a0, float a1) {
    vec2 p0 = vec2(r * cos(a0), r * sin(a0));
    vec2 p1 = vec2(r * cos(a1), r * sin(a1));

    float a = mod(atan(p.y, p.x), TAU);
    if(a > a0 && a < a1)
        return abs(length(p) - r);
    return min(length(p0 - p), length(p1 - p)); 
    //return 100.;
}

float sdArcElips(vec2 p, float a, float b, float t0, float t1) {
    float t = mod(atan(p.y / b, p.x / a), TAU);
    if(t > t0 && t < t1)
        return abs(p.x * p.x / a / a + p.y * p.y / b / b - 1.);
    return +100.;
}

float sdArc2(vec2 p, float r, float a0, float a1) {
    vec2 p0 = vec2(r * cos(a0), r * sin(a0));
    vec2 p1 = vec2(r * cos(a1), r * sin(a1));

    float a = mod(atan(p.y, p.x), TAU);
    if(!(a > a0 && a < a1))
        return abs(length(p) - r);
    return min(length(p0 - p), length(p1 - p));
}

float sdSegment(in vec2 p, in vec2 a, in vec2 b) {
    vec2 pa = p - a, ba = b - a;
    float h = clamp(dot(pa, ba) / dot(ba, ba), 0.0, 1.0);
    return length(pa - ba * h);
}

// signed distance to a 2D triangle
float sdTriangleIsosceles(in vec2 p, in vec2 q) {
    p.x = abs(p.x);
    vec2 a = p - q * clamp(dot(p, q) / dot(q, q), 0.0, 1.0);
    vec2 b = p - q * vec2(clamp(p.x / q.x, 0.0, 1.0), 1.0);
    float k = sign(q.y);
    float d = min(dot(a, a), dot(b, b));
    float s = max(k * (p.x * q.y - p.y * q.x), k * (p.y - q.y));
    return sqrt(d) * sign(s);
}

vec3 curColor(float d1, vec3 col, vec3 col1) {
    float s1 = smoothstep(0., -0.003, d1);
    return mix(col, col1, vec3(s1));
}

vec3 curColorN(float d1, vec3 col, vec3 col1, float dd, vec3 rd) {
    float s1 = smoothstep(0., -0.003, d1);
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
    //vec3 rd = vec3(0., 0., -1.);
    float difu = dot(norm, light);
    vec3 R1 = reflect(light, norm);
    float shininess = 10.0;
    float specular = pow(max(dot(R1, rd), 0.), shininess);
    vec3 colccs = col1 * difu + 0.8 * specular;      
    //return  (s1>0.) ? s1*col1*cs: col;
    return mix(col, colccs, vec3(s1));
}

float ticolor() {
    return (1. + (0.5 - noise(iTime)) * 1.18);
}

vec3 draw_cyrcle(vec2 p, float r, float dd, vec3 col, vec3 col1) {
    float d1 = abs(length(p) - r) - dd;
    return curColor(d1, col, col1);

}

float prism(vec2 p) {
    p -= 0.5;
    float t = max(abs(p.x), abs(p.y));
    return 1.0 - t;
}

float glass_h(vec2 p) {

    float n = 25., k = 1. / n * 0.1;
    p = fract(p * n);
    float t = prism(p) * k;
    return t;

}

vec3 glass_norm(vec2 p) {
    float h = 0.00001, dx = glass_h(p + vec2(h, 0)) - glass_h(p - vec2(h, 0)), dy = glass_h(p + vec2(0, h)) - glass_h(p - vec2(0, h)), dz = -2. * h;
    return -normalize(vec3(dx, dy, dz));
}

vec3 draw_disc(vec2 p, float r, vec3 col, vec3 col1) {
    float d1 = length(p) - r;
    col1 = curColor(d1, col, col1);
    return col1;
}

vec3 draw_disc_glass(vec2 p, float r, vec3 col, vec3 col1) {
    float d1 = length(p) - r;

    vec3 norm = glass_norm(p);
    vec3 light = normalize(vec3(sin(iTime / 3.), cos(iTime / 3.), 1.));
	//vec3 light = normalize(vec3(1., 1., 1.));
    vec3 rd = vec3(0., 0., -1.);
    float difu = dot(norm, light);
    vec3 R1 = reflect(light, norm);
    float shininess = 5.0;
    float specular = pow(max(dot(R1, rd), 0.), shininess);
    col1 = col1 * difu + 0.8 * specular;
    col1 = pow(col1, vec3(1.));
    col1 = curColor(d1, col, col1);
    return col1;
}

float drog(vec2 p) {
    float x0 = 0.04, w0 = 0.3, w1 = 0.27, y1 = 0.125;
    float d = sdSegment(p, vec2(x0, -x0), vec2(x0 + w0, -x0));

    float r0 = 0.025, alf0 = atan(y1, w0 - w1) + TAU, v = (length(vec2(w0 - w1, y1)) - r0) / alf0, dlt = v * TAU;
    vec2 pp = p - vec2(x0 + w1, -x0 - y1);
    float alf1 = mod(atan(pp.y, pp.x), TAU), alf2 = min(alf0, alf1 + TAU);
    vec2 p0 = vec2(r0, 0.), p1 = (r0 + v * alf1) * vec2(cos(alf1), sin(alf1)), p2 = (r0 + v * alf2) * vec2(cos(alf2), sin(alf2));
    float d2 = length(pp - p0);
    d2 = min(d2, length(pp - p1));
    d2 = min(d2, length(pp - p2));
    d = min(d, d2);
    return d;
}

vec3 ya(vec2 p, vec3 col) {
    //p -= 0.5;
    //vec3 col = vec3(1.);
    float d = drog(p);
    d = min(d, drog(vec2(-p.y, -p.x)));
    float eps = 0.001, h = 0.02;
    float s = smoothstep(eps, 0., d-h);
    col = mix(col, vec3(0.), s);
    return col;
}

vec3 ya2(vec2 p) {
    p -= 0.5;
    vec3 col = draw_cyrcle(p, 0.49, 0.01, vec3(1.), vec3(0.));
    float a = mod(atan(p.y, p.x), TAU);
    float n = 4.; 
    float nn = floor(a/TAU*n);
    p.xy *= rot(-nn*TAU/n);
    p.xy *= rot(-PI/2.);

    col = ya(p, col);
    return col;
}

vec3 vi4(vec2 p, float k) {
    p.x *= k;

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

    vec3 col = col4 * col4;//texture(iChannel0, p).rgb;
    vec2 pp = p - 0.5;
    float r0 = 0.1, dlt = 0.05, v = dlt / TAU, h = 0.015, eps = 0.001, alf = mod(atan(pp.y, pp.x), TAU), r = r0 + v * alf, d1 = floor(abs((length(pp) - r)) / dlt) * dlt + r, d2 = d1 + dlt, d = min(abs(length(pp) - d1), abs(length(pp) - d2));
    d = min(d, length(pp - vec2(r0, 0.)));
    float s = smoothstep(0., -eps, d - h);
    col = mix(col, col6, s);

    s = smoothstep(3. * eps, 0., abs(d - h));
    col = mix(col, vec3(0.), s);

    return col;
}
vec3 draw_pixelLine(vec2 p, vec2 a, vec2 b, vec3 col, vec3 col1)
{
    float s = 0.;
    

    if ((
        abs((p.x - a.x)/(p.y - a.y) - (b.x - a.x)/(b.y - a.y)) < 0.01 && 
    length(p-a) <= length(b-a) && dot(p-a, b-a) > 0.
    )
    || p == a || p == b
    ) 
        s = 1.;
    col = mix(col, col1, s);
    return col;    
}

vec3 moz(vec2 p)
{
    vec3 col0 = vec3(1.), col1 = vec3(1., 0., 0.), col = col0;
    float dlt = 1./51., x = floor(p.x/dlt), y = floor(p.y/dlt);
    
    col = draw_pixelLine(vec2(x, y), vec2(1., 25.), vec2(25., 49.), col, col1);
    col = draw_pixelLine(vec2(x, y), vec2(1., 25.), vec2(25., 1.), col, col1);
    col = draw_pixelLine(vec2(x, y), vec2(25., 1.), vec2(49., 25.), col, col1);
    col = draw_pixelLine(vec2(x, y), vec2(25., 49.), vec2(49., 25.), col, col1);
    //col = draw_pixelLine(vec2(x, y), vec2(0., 9.), vec2(1., 8.), col, col1);
    return col;

}
//https://www.xposz.shop/?ggcid=336928
//https://ca.pinterest.com/pin/557109416405805140/
void mainImage(out vec4 fragColor, in vec2 fragCoord) {
	//vec2 p = (-iResolution.xy + 2.0 * fragCoord) / iResolution.xy;
    vec2 p = fragCoord / iResolution.xy;
    float k = 1.0;
    float n = iResolution.x / iResolution.y / k;
    p.x = fract(p.x * n);

    vec3 col = ya2(p);

    fragColor = vec4(col, 1.0);
}
/////=====================================================================================
void main() {
    vec4 fragColor = vec4(0);
    mainImage(fragColor, gl_FragCoord.xy);
    gl_FragColor = fragColor;
}