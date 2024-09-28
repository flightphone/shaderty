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
umbrella 2
sdf,raymarching,pattern,umbrella,noise
simple umbrella
*/

float hash (vec3 p) {
    return fract(sin(dot(p, vec3(127.1,311.7, 74.7))) * 43758.5453123);
}

#define PI  3.14159265359
#define TAU 6.28318530718
#define rot(f) mat2(cos(f), -sin(f), sin(f), cos(f))

const float dist_infin = 10.0;
#define nn 128
const float eps = 0.001;



vec3 getSg(vec3 p, float nseg) {
    float fi = mod(atan(p.y, p.x), TAU);
    fi = mod(fi + PI / nseg, TAU);
    float n = floor(fi / TAU * nseg);
    p.xy *= rot(-n * TAU / nseg);
    return p;
}
float sdSegment(in vec2 p, in vec2 a, in vec2 b) {
    vec2 pa = p - a, ba = b - a;
    float h = clamp(dot(pa, ba) / dot(ba, ba), 0.0, 1.0);
    return length(pa - ba * h);
}

float detail = 0.;

float umb(vec3 p) {
    p.yz*=rot(PI/2.);
    float n1 = 14., R = 2., h = 2.1;
    //============================handle=======================
    float d = length(vec2(length(vec2(p.xy)), max(abs(p.z) - h, 0.))) - 0.07;
    detail = 0.;
    float d1 = length(vec2(length(vec2(p.xy)), 
    p.z - clamp(p.z, -h, -0.8*h))) - 0.12;
    if (d1 < d)
    {
        d = d1;
        detail = 1.;
    }
    
    p = getSg(p, n1);
    //============================spokes=====================
    float dz = 0.5, dr = sqrt(R*R - dz*dz);
    vec2 a0 = vec2(dr * cos(PI / n1), dr * sin(PI / n1)*sign(p.y));
    float d2 = sdSegment(p.xy, a0, vec2(0.));
    d2 = length(vec2(d2, p.z-dz)) - 0.02;
    if (d2 < d)
    {
        d = d2;
        detail = 1.;
    }
    //=======================================================

    //===================tent=========================
    vec3 p0 = normalize(p);
    vec2 a = vec2(R * cos(PI / n1), R * sin(PI / n1)), b = vec2(R * cos(PI / n1), -R * sin(PI / n1));
    float t = 0.;
    if(p.z > 0.) {
        t = sqrt(R * R / (p0.z * p0.z + p0.x * p0.x / cos(PI / n1) / cos(PI / n1)));
        float l = p0.x * t / cos(PI / n1);
        a = vec2(l * cos(PI / n1), l * sin(PI / n1));
        b = vec2(l * cos(PI / n1), -l * sin(PI / n1));
    }
    float d3 = sdSegment(p.xy, a, b);
    d3 = 0.95*length(vec2(d3, p.z - p0.z * t))-0.02;
    if (d3 < d)
    {
        d = d3;
        detail = 2.;
    }
    //=================================================

    //=====================edges========================
    p.xy *= rot(-PI / n1 * sign(p.y));
    float d4 = abs(length(p.xz) - R);
    if(p.z < 0.)
        d4 = length(p.xz - vec2(R, 0));
    d4 = length(vec2(d4, p.y)) - 0.03;
    if (d4 < d)
    {
        d = d4;
        detail = 2.;
    }
    return d;

}


float map(vec3 p) {
    return umb(p);
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
vec3 col0 = vec3(0.73, 0.59, 0.3);
vec3 col1 = vec3(0.75, 0.75, 0.75);
vec3 col2 = vec3(0.7, 0.7, 1.);
vec3 col4 = vec3(0.9, 0.1, 0.1);
//vec3 resColor = vec3(0.73, 0.59, 0.3);

void mainImage(out vec4 fragColor, in vec2 fragCoord) {
    vec3 light = normalize(vec3(0.0, .0, 1.)); //light
    vec3 light2 = normalize(vec3(0.0, 0.0, -1.)); //light
    vec2 mo = 1.5 * cos(0.5 * iTime + vec2(0, 11));
    //if  (iMouse.z > 0.0)
    {
        mo = (-iResolution.xy + 2.0 * (iMouse.xy)) / iResolution.y;
        mo*=1.7;
    }
    vec3 ro = vec3(0.0, 0.0, 5.); // camera
    //camera rotation
    ro.yz *= rot(mo.y);
    ro.xz *= rot(-mo.x - 1.57);

    const float fl = 1.5; // focal length
    float dist = dist_infin;

    vec3 b1 = vec3(0.23529411764705882, 0.4235294117647059, 0.7725490196078432), b2 = vec3(0.3686274509803922, 0.5725490196078431, 0.8941176470588236);
    vec3 bg = mix(b2, b1, fragCoord.y / iResolution.y);  
    //antialiasing
    vec3 tot = vec3(0.0);
    for(int m = 0; m < AA; m++) for(int n = 0; n < AA; n++) {
            vec2 o = vec2(float(m), float(n)) / float(AA) - 0.5;
            vec2 p = (-iResolution.xy + 2.0 * (fragCoord + o)) / iResolution.y;
            vec3 rd = GetRayDir(p, ro, vec3(0, 0., 0), fl); //ray direction
            vec3 col = bg; // background  
            
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
                if (detail == 0.) col = col0;
                if (detail == 1.) col = col1;
                if (detail == 2.) 
                {
                    col = col2*col2;
                    vec3 bp = floor(pos*3./3.);
                    float h = hash(bp);
                    if (h < 0.1)
                        col = col4;
                }
                vec3 nor = calcNormal(pos);
                vec3 R = reflect(light, nor);
                float specular = pow(max(abs(dot(R, rd)), 0.), 16.);
                float difu = abs(dot(nor, light));
                col = col * (clamp(difu, 0., 1.0) + 0.5) + vec3(1., .7, .4) * specular;
                float fre = pow(clamp(dot(nor, rd) + 1., .0, 1.), 3.); // Fresnel, for some mild glow.
                col += vec3(.1, .1, 0.1) * fre; //?
                col = sqrt(col);
            }
            //==========================raymatch=============================
            tot += col;
        }
    tot = tot / float(AA) / float(AA);
    //tot = tot / float(AA);
    //antialiasing
    fragColor = vec4(tot, 1.0);
}
/////=====================================================================================
void main() {
    vec4 fragColor = vec4(0);
    mainImage(fragColor, gl_FragCoord.xy);
    gl_FragColor = fragColor;
}