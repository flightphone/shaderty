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
vec3 col1 = vec3(0.3764705882352941, 0.8196078431372549, 0.37254901960784315);
vec3 col2 = vec3(0.8117647058823529, 0.17647058823529413, 0.807843137254902);

//https://iquilezles.org/articles/smin/
// polynomial smooth min 1 (k=0.1)
float smin(float a, float b, float k) {
    float h = clamp(0.5 + 0.5 * (b - a) / k, 0.0, 1.0);
    return mix(b, a, h) - k * h * (1.0 - h);
}

float rand(float t) {
    return fract(sin(t * 213456.12234));
}

float heirw(vec3 p, float h, float r, float r0, float fi) {
    h = h * (1.0 + 0.05 * cos(iTime * 4. + fi));
    float z = clamp(p.z, 0., h), // radius pimple 
    x = sin(z * PI * 2.0 + iTime * 4. + fi) * h * 0.3 * z * (h - z) / h / h, y = sin(z * PI * 2.0 - iTime * 4. + fi) * h * 0.3 * z * (h - z) / h / h;

    vec3 p2 = vec3(x, y, z);
    //Color
    sdfColor = mix(col1, col2, pow(vec3(p.z / h), vec3(3.)));
    return length(p - p2) * 0.5 - r * (h - z) / h - r0;
}

float bbody(vec3 p) {

    float h = 1.1, //height pimple
    r = 0.3, n = 15., m = 6., z = clamp(p.z, r, h + r);
    vec3 p2 = vec3(0., 0., z);
    float dz = length(p - p2) - r;

    float dlon = TAU / n, dlat = h / m, l = length(p.xy), lon = mod(atan(p.y, p.x), TAU), lat = p.z - r, //longitude and latitude
    i = floor(lon / dlon), j = clamp(floor(lat / dlat), 0., m), dp = dz, x1 = lon / TAU, y1 = clamp(p.z, 0., h + 2. * r) / (h + 2. * r);

    //calc three row

    for(float kj = 1.; kj < 2.; kj++) for(float ki = 1.; ki < 2.; ki++) {
            float j1 = clamp(j + kj - 1., 0., m), i1 = mod(i + ki - 1., n), lon1 = i1 * dlon + 0.5 * dlon, dx = (lon - lon1) * r, lat1 = j1 * dlat + 0.5 * dlat, //longitude and latitude nearest pimple
            dy = lat - lat1, num = (i1 + 1.) * m + (j1 + 1.), fi = rand(num) * PI, d = heirw(vec3(dx, dy, dz), 0.5, 0.01, 0.001, fi);
            if(d < dp) {
                resColor = sdfColor;
                dp = d;
            }
        }

    //texture
    vec3 cl = texture(iChannel1, vec2(x1, y1)).rgb;
    float disp = dot(cl, vec3(0.3, 0.59, 0.11));
    disp *= r * 0.1;
    dz -= disp;
    if(dz < dp)
        resColor = cl;

    return smin(dp, dz, 0.01);

}

float sdBact(vec3 p) {
    float d = dist_infin;
    if(p.z < 0.) {
        p.z *= -1.;
        d = heirw(p, 1.8, 0.03, 0.003, 0.);
        resColor = sdfColor;
    } else {
        d = bbody(p);
    }
    return d;
}

float map(in vec3 pos) {
    pos.xz *= rot(PI / 2.);
    return sdBact(pos);
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

struct HIT {
    float dist;
    vec3 nor;
    vec3 pos;
};
const HIT hit_inf = HIT(dist_infin, vec3(0.0), vec3(0.0));

HIT giper3D(vec3 ro, vec3 rd) {
    float t = 0.;
    for(int i = 0; i < nn; i++) {
        vec3 pos = ro + rd * t;
        float h = map(pos);
        if(h < eps || t >= dist_infin)
            break;
        t += h;
    }

    if(t >= dist_infin)
        return hit_inf;

    vec3 pos = ro + t * rd;
    vec3 nor = calcNormal(pos);
    return HIT(t, nor, pos);
}

vec3 GetRayDir(vec2 uv, vec3 p, vec3 l, float z) {
    vec3 f = normalize(l - p), r = normalize(vec3(f.z,0,-f.x)), u = cross(f, r), c = f * z, i = c + uv.x * r + uv.y * u;
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
   
    vec3 b1 = vec3(0.050980392156862744, 0.2980392156862745, 0.47058823529411764), b2 = vec3(0.3764705882352941, 0.7529411764705882, 0.8784313725490196), bg = mix(b1, b2, vec3((1.0 - abs(fragCoord.x - iResolution.x / 2.) / iResolution.y * 2.) * fragCoord.y / iResolution.x));   

    //antialiasing
    vec3 tot = vec3(0.0);
    for(int m = 0; m < AA; m++) for(int n = 0; n < AA; n++) {
            vec2 o = vec2(float(m), float(n)) / float(AA) - 0.5;
            vec2 p = (-iResolution.xy + 2.0 * (fragCoord + o)) / iResolution.y;
            vec3 rd = GetRayDir(p, ro, vec3(0, 0., 0), fl); //ray direction
            
            
            vec3 col = bg * bg; // background  
            HIT giper = giper3D(ro, rd);
            if(giper.dist < dist) {
                col = resColor;
                col = calccolor(col, col, -rd, light, light2, giper.nor);
            }
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