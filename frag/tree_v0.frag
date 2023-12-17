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

float spine(vec3 p, float l, float x, float z, float th) {
    vec3 s = vec3(sin(z) * cos(x), sin(z) * sin(x), cos(z));
    float d = clamp(dot(p, s), 0., l);
    float r1 = 0.03, r2 = 0.005, r = r1 + (r2 - r1) * d / l;
    r *= th;
    return length(p - d * s) - r;

}

float getlon(float lon, float n, float shift) {
    lon = mod(lon - shift, TAU);
    float dlon = TAU / n, lon1 = floor(lon / dlon) * dlon;
    if((lon - lon1) > dlon / 2.)
        lon1 += dlon;
    return lon1 + shift;
}

float branch(vec3 p, float l, float r, float ls, float lss, float fi, float n, float th, int prec) {
    float z = clamp(p.z, 0., l);
    r -= r * z / l*0.6;
    float d = length(p - vec3(0, 0., z)) - r;
    sdfColor = col2;
    //float ls = 1.7, lss = 1., n = 30., fi = PI / 6.; 
    float lon = mod(atan(p.y, p.x), TAU), dlon = TAU / n;
    float j = floor(z / lss);
    float h1 = j * lss, h2 = h1 + lss, shift1 = mod(j, 2.) * dlon / 2., shift2 = mod((j + 1.), 2.) * dlon / 2.;
    float h3 = h1 - lss, h4 = h2 + lss, shift3 = mod(j - 1., 2.) * dlon / 2., shift4 = mod((j + 2.), 2.) * dlon / 2.;

    float lon1 = getlon(lon, n, shift1), lon2 = getlon(lon, n, shift2);
    float lon3 = getlon(lon, n, shift3), lon4 = getlon(lon, n, shift4);

    float d2 = dist_infin;

    d2 = spine(p - vec3(r * cos(lon1), r * sin(lon1), h1), ls, lon1, fi, th);
    if(d2 < d) {
        d = d2;
        sdfColor = col1;
    }
    if(h2 < l) {
        d2 = spine(p - vec3(r * cos(lon2), r * sin(lon2), h2), ls, lon2, fi, th);
        if(d2 < d) {
            d = d2;
            sdfColor = col1;
        }
    }
    
    if(prec == 1) {
        if(h3 >= 0.) {
            d2 = spine(p - vec3(r * cos(lon3), r * sin(lon3), h3), ls, lon3, fi, th);
            if(d2 < d) {
                d = d2;
                sdfColor = col1;
            }
        }
        if(h4 < l) {
            d2 = spine(p - vec3(r * cos(lon4), r * sin(lon4), h4), ls, lon4, fi, th);
            if(d2 < d) {
                d = d2;
                sdfColor = col1;
            }
        }
    }
    

    return d;

}

float branch2(vec3 p, float l) {
    //float l = 1.2;
    //float f = PI / 8. * p.z / l;
    //p.yz *= rot(f);
    float d = branch(p, l, 0.025, 0.25, 0.15, PI / 6., 20., 0.3, 1);
    vec3 col = sdfColor;
    p -= vec3(0., 0., l / 2.);
    vec3 p2 = p;
    p2.xz *= rot(PI / 4.);
    float d2 = branch(p2, l / 2., 0.025, 0.25, 0.15, PI / 6., 20., 0.3, 1);
    if(d2 < d) {
        col = sdfColor;
        d = d2;
    }
    p2 = p;
    p2.xz *= rot(-PI / 4.);
    d2 = branch(p2, l / 2., 0.025, 0.25, 0.15, PI / 6., 20., 0.3, 1);
    if(d2 < d) {
        col = sdfColor;
        d = d2;
    }
    sdfColor = col;

    return d;
}

float labranch(vec3 p, float l, float h, float lon)
{
    p -= vec3(0., 0., h);
    p.xy *= rot(lon);
    p.yz *= rot(-PI / 2.);
    
    return branch2(p, l);
}

float tree(vec3 p)
{
    float l = 5.;
    float d = branch(p, l, 0.05, 0.25, 0.3, PI/20., 6., .5, 0);
    vec3 col = sdfColor;
    
    float lss = 0.4;

    float n = 8.0, lon = mod(atan(p.y, p.x), TAU), dlon = TAU / n;
    float z = clamp(p.z, 0., l), j = floor(z / lss);
    float h1 = j * lss, h2 = h1 + lss, shift1 = mod(j, 2.) * dlon / 2., shift2 = mod((j + 1.), 2.) * dlon / 2.;
    //float h3 = h1 - lss, h4 = h2 + lss, shift3 = mod(j - 1., 2.) * dlon / 2., shift4 = mod((j + 2.), 2.) * dlon / 2.;

    float lon1 = getlon(lon, n, shift1), lon2 = getlon(lon, n, shift2);
    //float lon3 = getlon(lon, n, shift3), lon4 = getlon(lon, n, shift4);

    float d2 = dist_infin;
    d2 = labranch(p, 1.5, h1, lon1); 
    if(d2 < d) {
        d = d2;
        sdfColor = col1;
    }
    if(h2 < l) {
        d2 = labranch(p, 1.5, h2, lon2); 
        if(d2 < d) {
            d = d2;
            sdfColor = col1;
        }
    }
    
    sdfColor = col;
    return d;
}

float map(vec3 p) {
    //float d = branch(p + vec3(0.0, 0.0, .5), 1.2, 0.025, 0.25, 0.15, PI/6., 20., 0.3);

    float d = branch2(p + vec3(0.0, 0.0, 2.), 4.);
    //p +=vec3(0.0, 2.5, 0.);
    //float d = tree(p.xzy);
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
    vec3 ro = vec3(0.0, 0.0, 3.5); // camera
    //camera rotation
    ro.yz *= rot(mo.y * PI);
    //ro.z += 4.*mo.y;
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