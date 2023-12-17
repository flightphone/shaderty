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

const float dist_infin = 20.0;
#define nn 128
const float eps = 0.001;

vec3 sdfColor;
vec3 resColor;
vec3 col1 = vec3(0.13725490196078433, 0.4823529411764706, 0.28627450980392155);
vec3 col2 = vec3(0.3607843137254902, 0.16470588235294117, 0.027450980392156862);

float branch2(vec2 p, float l) {
    p *= rot(PI / 3.);
    float x = clamp(p.x, 0., l);
    float d = length(p - vec2(x, 0.));
    return d;

}

float branch1(vec2 p, float l) {
    p *= rot(PI / 3.);
    float n = 2.;
    float l2 = 0.3;
    float x = clamp(p.x, 0., l);
    float d = length(p - vec2(x, 0.));
    
    float dlat = l / n, k = x / dlat, lat = clamp(floor(k) * dlat + dlat / 2.0, dlat / 2., l - dlat / 2.), dz = p.y;
    float lt = clamp(lat + sign(x - lat) * dlat, dlat / 2., l - dlat / 2.), lt1 = min(lat, lt), lt2 = max(lat, lt);
    float d2 = dist_infin;
    
    for(float i = 0.; i < 1.; i++) {
        float llt = clamp(lt1 - dlat * i, dlat / 2., l - dlat / 2.);
        d2 = branch2(vec2(p.x - llt, dz), l2);
        d = min(d, d2);
        d2 = branch2(vec2(p.x - llt, -dz), l2);
        d = min(d, d2);

        llt = clamp(lt2 + dlat * i, dlat / 2., l - dlat / 2.);
        d2 = branch2(vec2(p.x - llt, dz), l2);
        d = min(d, d2);
        d2 = branch2(vec2(p.x - llt, -dz), l2);
        d = min(d, d2);
    }
    return d;
}

float branch0(vec2 p, float l) {

    float n = 3.;
    float l2 = 1.5;
    float x = clamp(p.x, 0., l);
    float d = length(p - vec2(x, 0.));
   

    float dlat = l / n, k = x / dlat, lat = clamp(floor(k) * dlat + dlat / 2.0, dlat / 2., l - dlat / 2.), dz = p.y;
    float lt = clamp(lat + sign(x - lat) * dlat, dlat / 2., l - dlat / 2.), lt1 = min(lat, lt), lt2 = max(lat, lt);
    float d2 = dist_infin;
    
    for(float i = 0.; i < 2.; i++) {
        float llt = clamp(lt1 - dlat * i, dlat / 2., l - dlat / 2.);
        d2 = branch1(vec2(p.x - llt, dz), l2);
        d = min(d, d2);
        d2 = branch1(vec2(p.x - llt, -dz), l2);
        d = min(d, d2);

        llt = clamp(lt2 + dlat * i, dlat / 2., l - dlat / 2.);
        d2 = branch1(vec2(p.x - llt, dz), l2);
        d = min(d, d2);
        d2 = branch1(vec2(p.x - llt, -dz), l2);
        d = min(d, d2);
    }
    return d;
}

float map(vec3 p) {
    float d = branch0(p.xy, 3.);
    resColor = col2;
    d = length(vec2(p.z, d));
    return d - 0.03;
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
    vec3 ro = vec3(0.0, 0.0, 3.5); // camera
    //camera rotation
    ro.yz *= rot(mo.y * PI);
    ro.xz *= rot(-mo.x * TAU);

    const float fl = 1.5; // focal length
    float dist = dist_infin;

    vec3 b1 = vec3(0.6235294117647059, 0.8, 0.9803921568627451), b2 = vec3(0.49019607843137253, 0.6980392156862745, 0.9568627450980393);
    vec3 bg = mix(b1, b2, vec3(fragCoord.y / iResolution.y));   
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
                col = calccolor(col, col, -rd, light, light2, nor);

            }
            //==========================raymatch=============================
            tot += col;
        }
    tot = sqrt(tot) / float(AA);
    //antialiasing
    fragColor = vec4(tot, 1.0);
}
/////=====================================================================================
void main() {
    vec4 fragColor = vec4(0);
    mainImage(fragColor, gl_FragCoord.xy);
    gl_FragColor = fragColor;
}