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

vec3 col1 = vec3(0., 0.878, 0.568);
vec3 col2 = vec3(0.7686274509803922, 0.8235294117647058, 0.8745098039215686);
vec3 col3 = vec3(1., 0.8431, 0.);
float sdfReflect = 0.5;
float resReflect = 0.5;

vec3 csky(vec3 p) {
    float n = 5., m = 5., dlat = PI / n, dlon = TAU / m;
    float lon = mod(atan(p.y, p.x), TAU), lat = atan(length(p.xy), p.z);
    float fo = fract(lon / dlon), fa = fract(lat / dlat);

    float pst = fo * fa * (1. - fo) * (1. - fa);
    pst = smoothstep(0.0, 0.0625, pst);
    pst = clamp(pst, 0.1, 1.0);
    return vec3(pst);
}

float sdSolidAngle(vec3 p, vec2 c, float ra) {
  // c is the sin/cos of the angle
    vec2 q = vec2(length(p.xz), p.y);
    float l = length(q) - ra;
    float m = length(q - c * clamp(dot(q, c), 0.0, ra));
    return max(l, m * sign(c.y * q.x - c.x * q.y));
}

float heigthBranch(vec2 p) {
    float n = 2.5;
    float df = PI / n / 2.5;
    float fi = atan(p.y, p.x);
    float L = length(p.xy);
    float r = cos(n * fi);
    if(abs(fi) > df)
        r = 0.;
    float d = r - L;
    float h = smoothstep(0., 0.3, d * L * L);
    if(h > 0.) {

        sdfColor = col1;
        float pst = smoothstep(0.2, 0., abs(L - 0.6));
        sdfColor = mix(col1, col2, pst);
        sdfReflect = mix(0.2, 0., pst);
    }
    return h;
}

float getlon(float lon, float n, float shift) {
    lon = lon - shift;
    float dlon = TAU / n, lon1 = floor(lon / dlon) * dlon;
    if((lon - lon1) >= dlon / 2.)
        lon1 += dlon;
    return lon1 + shift; ////mod(lon1 + shift, TAU);
}

float sdTree(vec3 p, float l, float r) {
    float mfi = PI / 8.;
    float d = sdSolidAngle(p, vec2(sin(mfi), cos(mfi)), l) - r;
    if(p.y < 0. || p.y > l * cos(mfi)) {
        sdfColor = col2;
        sdfReflect = 0.;
        return d;
    }
    sdfColor = col1;
    sdfReflect = 0.1;

    float n = 8., m = 5., nc = 6., dnc = l / nc;
    float lss = l / 2. / m, ls = 2. * lss;
    float z = clamp(p.y, 0., l);
    float lon = mod(atan(p.z, p.x), TAU), dlon = TAU / n;

    float j = floor(z / lss);
    float h1 = j * lss, shift1 = mod(j, 2.) * dlon / 2.;//,h2 = h1 + lss, shift2 = mod((j + 1.), 2.) * dlon / 2.;
    float h3 = h1 - lss, shift3 = mod(j - 1., 2.) * dlon / 2.;//h4 = h1 - 2.*lss, shift4 = mod((j - 2.), 2.) * dlon / 2.;

    float lon1 = getlon(lon, n, shift1);//, lon2 = getlon(lon, n, shift2);
    float lon3 = getlon(lon, n, shift3);//, lon4 = getlon(lon, n, shift4);

    float h = 0.;
    if(j < n && h1 > 0.)
        h = max(heigthBranch(vec2((p.y - h1) / ls, (lon - lon1) / dlon * 0.5)) * l, h);
    if(h3 > 0. && h3 + ls < l)
        h = max(heigthBranch(vec2((p.y - h3) / ls, (lon - lon3) / dlon * 0.5)) * l, h);

    //shpere

    float tj = j;
    float h2 = h1 + lss, shift2 = mod((j + 1.), 2.) * dlon / 2.;
    float lon2 = getlon(lon, n, shift2);
    float hp = h1, lonsp = lon1;
    if(h2 - z < z - h1) {
        hp = h2;
        lonsp = lon2;
        tj = j + 1.;
    }
    if(tj > 2.) {
        float dx = hp * tan(mfi) * (lon - lonsp), dy = (z - hp) / cos(mfi), dr = l / 15.;
        float ra = length(vec2(dx, dy) / dr);
        if(ra < 0.4) {
            h = max(sqrt(0.16 - ra * ra), h);
            sdfColor = vec3(0.698, 0.098, 0.176);
            sdfReflect = 0.4;

        }
    }

    float pst = smoothstep(0.1, 0., fract(z / dnc));
    sdfColor = mix(sdfColor, col3, pst);
    sdfReflect = mix(sdfReflect, 0.8, pst);

    return d * 0.3 - h * 0.06 * sqrt(z / l);

}

float map(vec3 p) {
    float l = 2.3;
    p.xy *= rot(PI);
    p += vec3(0., l / 2., 0.);
    p.xz *= rot(iTime / 2.);
    float d = sdTree(p, l, 0.05);
    resColor = sdfColor;
    resReflect = sdfReflect;
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
    vec3 light = normalize(vec3(1.0, .0, -2.5)); //light
    vec3 light2 = normalize(vec3(-1.0, -.0, 2.5)); //light
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
                col = resColor;
                vec3 nor = calcNormal(pos);

                //reflection

                vec3 psk = reflect(rd, nor);
                vec3 c2 = csky(psk);

                col = calccolor(col, col, -rd, light, light2, nor);
                col = mix(col, c2, resReflect);

                //col += c2*0.1;

            }
            //==========================raymatch=============================
            tot += col;
        }
    //tot = sqrt(tot) / float(AA);
    tot = tot / float(AA);
    //antialiasing
    fragColor = vec4(tot, 1.0);
}
/////=====================================================================================
void main() {
    vec4 fragColor = vec4(0);
    mainImage(fragColor, gl_FragCoord.xy);
    gl_FragColor = fragColor;
}