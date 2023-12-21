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

const float dist_infin = 10.0;
#define nn 128
const float eps = 0.001;

vec3 sdfColor;
vec3 resColor;
//vec3 col1 = vec3(0.5019607843137255, 0.6705882352941176, 0.34509803921568627);
vec3 col1 = vec3(0.3137254901960784, 0.7843137254901961, 0.47058823529411764);
vec3 col2 = vec3(0.7686274509803922, 0.8235294117647058, 0.8745098039215686);
float resReflect = 0.5;

float smin(float a, float b, float k) {
    float h = clamp(0.5 + 0.5 * (b - a) / k, 0.0, 1.0);
    return mix(b, a, h) - k * h * (1.0 - h);
}
float sdCosNp(vec2 p, float a, float n) {
    float df = PI / n / 2.;
    float fi = atan(p.y, p.x);
    float L = length(p.xy);

    float d = dist_infin;
    
    

    float r = a * cos(n * fi);
    if(abs(fi) > df)
        r = -0.01;
    d = min(abs(L - r), d);

    if(L <= a) {
        //  float f = acos(L / a) / n;
        //  d = min(2.0 * abs(sin((fi - f) / 2.0)) * L, d);
        //  d = min(2.0 * abs(sin((fi + f) / 2.0)) * L, d);

        float f = acos(p.x)/n, y = p.x * tan(f);
        if (p.y > 0.)
            d = min(abs(p.y - y), d);
        else    
            d = min(abs(-p.y - y), d);

    }
    return d;

}

float sdCosNp2(vec2 p, float a, float n) {
    float df = PI / n / 2.;
    float fi = atan(p.y, p.x);
    float L = length(p.xy);

    float d = dist_infin;
    
    float ff = asin(p.y/a);
    d = min(abs(a*cos(ff) - p.x), d);
/*
    float r = a * cos(n * fi);
    if(abs(fi) > df)
        r = -0.01;
    d = min(abs(L - r), d);
*/
    //if(L <= a) {
        //  float f = acos(L / a) / n;
        //  d = min(2.0 * abs(sin((fi - f) / 2.0)) * L, d);
        //  d = min(2.0 * abs(sin((fi + f) / 2.0)) * L, d);

        float f = acos(p.x)/n, y = p.x * tan(f);
        if (p.y > 0.)
            d = min(abs(p.y - y), d);
        else    
            d = min(abs(-p.y - y), d);

    //}
    return d;

}

float exmp(vec3 p) {
    /*
    1.:0.6
    1.5:0.5
    2:0.45
    2.5:0.4
    3:0.3
    3.5:0.3
    4.:0.25
    4.5:0.25
    5.:
    */

    float d = sdCosNp2(p.xy, 1., 4.);
    d = length(vec2(d, p.z))*0.4  - 0.01;

    /*
    float d2 = sdCosNp(p.xy, 1.2, 1.5);
    d2 = length(vec2(d2, p.z)) * 0.5 - 0.01;
    d = min(d, d2);

    d2 = sdCosNp(p.xy, 1.2, 1.2);
    d2 = length(vec2(d2, p.z)) * 0.5 - 0.01;
    d = min(d, d2);

    d2 = sdCosNp(p.xy, 1.2, 2.5);
    d2 = length(vec2(d2, p.z)) * 0.4 - 0.01;
    d = min(d, d2);

    d2 = sdCosNp(p.xy, 1.2, 3.5);
    d2 = length(vec2(d2, p.z)) * 0.3 - 0.01;
    d = min(d, d2);
    */

    sdfColor = col1;
    return d;

}

vec3 csky(vec3 p) {
    float n = 5., m = 5., dlat = PI / n, dlon = TAU / m;
    float lon = mod(atan(p.y, p.x), TAU), lat = atan(length(p.xy), p.z);
    float fo = fract(lon / dlon), fa = fract(lat / dlat);

    float pst = fo * fa * (1. - fo) * (1. - fa);
    pst = smoothstep(0.0, 0.0625, pst);
    pst = clamp(pst, 0.1, 1.0);
    return vec3(pst);
}

float map(vec3 p) {

    float d = exmp(p);
    resColor = sdfColor;
    return d;
}

// https://iquilezles.org/articles/normalsSDF
vec3 calcNormal(in vec3 pos) {
    const float h = 0.0001; // replace by an appropriate value
    const vec2 k = vec2(1, -1);
    return normalize(k.xyy * map(pos + k.xyy * h) +
        k.yyx * map(pos + k.yyx * h) +
        k.yxy * map(pos + k.yxy * h) +
        k.xxx * map(pos + k.xxx * h));
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
#define AA 1

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

void mainImage(out vec4 fragColor, in vec2 fragCoord) {
    vec3 light = normalize(vec3(0.0, 1.0, -2.5)); //light
    vec3 light2 = normalize(vec3(0.0, -1.0, 2.5)); //light
    vec2 mo = vec2(0.0, 0.0);
    //if  (iMouse.z > 0.0)
    {
        mo = (-iResolution.xy + 2.0 * (iMouse.xy)) / iResolution.y;
    }
    vec3 ro = vec3(0.0, 0.0, 2.5); // camera
    //camera rotation
    ro.yz *= rot(mo.y * PI);
    ro.xz *= rot(-mo.x * TAU);

    const float fl = 1.5; // focal length
    float dist = dist_infin;

    vec3 b1 = vec3(0.23529411764705882, 0.4235294117647059, 0.7725490196078432), b2 = vec3(0.3686274509803922, 0.5725490196078431, 0.8941176470588236);
    //vec3 bg = mix(b1, b2, vec3((1.0 - abs(fragCoord.x - iResolution.x / 2.) / iResolution.y * 2.) * fragCoord.y / iResolution.x));   
    vec3 bg = mix(b2, b1, fragCoord.y / iResolution.y);   
    //antialiasing
    vec3 tot = vec3(0.0);
    for(int m = 0; m < AA; m++) for(int n = 0; n < AA; n++) {
            vec2 o = vec2(float(m), float(n)) / float(AA) - 0.5;
            vec2 p = (-iResolution.xy + 2.0 * (fragCoord + o)) / iResolution.y;
            vec3 rd = GetRayDir(p, ro, vec3(0, 0., 0), fl); //ray direction
            vec3 col = bg * bg; // background  
            //==========================raymatch=============================
            float td = 0.;
            vec3 pos = vec3(0.);
            for(int i = 0; i < nn; i++) {
                pos = ro + rd * td;
                float h = map(pos);
                if(h < eps || td >= dist_infin)
                    break;
                td += h;
            }
            if(td < dist_infin) {
                col = resColor;
                vec3 nor = calcNormal(pos);

                //reflection

                //vec3 psk = reflect(rd, nor);
                //vec3 c2 = csky(psk);

                col = calccolor(col, col, -rd, light, light2, nor);
                //col = mix(col, c2, resReflect);

                //col += c2*0.1;

            }
            //==========================raymatch=============================
            tot += col;
        }
    tot = sqrt(tot) / float(AA);
    //tot = pow(tot, vec3(0.7)) / float(AA);
    //antialiasing
    fragColor = vec4(tot, 1.0);
}
/////=====================================================================================
void main() {
    vec4 fragColor = vec4(0);
    mainImage(fragColor, gl_FragCoord.xy);
    gl_FragColor = fragColor;
}