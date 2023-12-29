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
vec3 col1 = vec3(0.3764, 0.8196, 0.3725);
vec3 col2 = vec3(0.8117, 0.1764, 0.8078);

//https://iquilezles.org/articles/smin/
// polynomial smooth min 1 (k=0.1)
float smin(float a, float b, float k) {
    float h = clamp(0.5 + 0.5 * (b - a) / k, 0.0, 1.0);
    return mix(b, a, h) - k * h * (1.0 - h);
}

float rand(float t) {
    return fract(sin(t * 213.12234));
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

    float lon1 = i * dlon + 0.5 * dlon, dx = (lon - lon1) * r, lat1 = j * dlat + 0.5 * dlat, //longitude and latitude nearest pimple
    dy = lat - lat1, num = (i + 1.) * m + (j + 1.), fi = rand(num) * PI, d = heirw(vec3(dx, dy, dz), 0.5, 0.01, 0.001, fi);
    if(d < dp) {
        resColor = sdfColor;
        dp = d;
    }

    //texture
    vec3 cl = texture(iChannel0, vec2(x1, y1)).rgb;
    float disp = dot(cl, vec3(0.3, 0.59, 0.11));
    disp *= r * 0.1;
    dz -= disp;
    if(dz < dp)
        resColor = cl;

    return smin(dp, dz, 0.01);

}

float map(vec3 p) {
    p.xz *=rot(PI/2.); //rotate object
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
#define AA 2

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

//https://www.shadertoy.com/view/dlKBRc
vec3 lightingv3(vec3 lightColor, vec3 rd, vec3 L, vec3 normal) 
{   
    vec3 V = rd;
    vec3 N = normal;
    vec3 R = reflect (-L, N);
    float shadow = 1.;
    float occ = 0.7;
    float Ka = 0.5;
    vec3 ambient = Ka + Ka * dot(normal, vec3(0., 1., 0.))*lightColor;
    ambient*=0.5;
    vec3 fresnel =  lightColor *  pow(clamp(1.0 + dot(rd, N), 0.0, 1.0), 2.0);;
    float diff= clamp(dot(N, L), 0., 1.0);
    vec3 diffuse =  lightColor * diff;
    float shininess=10.0;
    float specular    = pow(max(dot(R, V), 0.0), shininess);
    vec3 back = 0.5 * lightColor * clamp(dot(N, -L), 0.0, 1.0); // back
    vec3 colOut = occ*lightColor*(ambient+diffuse*shadow+.25 +back) + vec3(.5)*specular*specular;
    return colOut;
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

    vec3 b1 = vec3(0.0509, 0.2980, 0.4705), b2 = vec3(0.3764, 0.7529, 0.8784), bg = mix(b1, b2, vec3((1.0 - abs(fragCoord.x - iResolution.x / 2.) / iResolution.y * 2.) * fragCoord.y / iResolution.x));   
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
                //col = calccolor(col, col, -rd, light, light2, nor);
                col = lightingv3(col, -rd, light, nor);


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