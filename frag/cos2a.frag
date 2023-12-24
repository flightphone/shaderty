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
vec3 col1 = vec3(.2);
//vec3 col1 = vec3(0.3137254901960784, 0.7843137254901961, 0.47058823529411764);
//vec3 col2 = vec3(0.7686274509803922, 0.8235294117647058, 0.8745098039215686);
vec3 col2 = vec3(0.4745098039215686, 0.26666666666666666, 0.23137254901960785);
vec3 col3 = vec3(1., 0.8431, 0.);
float resReflect = 0.5;
float sdfReflect;

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
          float f = acos(L / a) / n;
          d = min(2.0 * abs(sin((fi - f) / 2.0)) * L, d);
          d = min(2.0 * abs(sin((fi + f) / 2.0)) * L, d);
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


        float f = acos(p.x)/n, y = p.x * tan(f);
        if (p.y > 0.)
            d = min(abs(p.y - y), d);
        else    
            d = min(abs(-p.y - y), d);


    return d;

}

float knot(vec3 p)
{
    float z = clamp(p.z, -1., 1.);
    float d = dist_infin;
    float f = asin(-z)/3.0;
    for (float i = 0.; i < 3.; i++)
    {
        float fi = f + TAU/3.*i;
        vec3 v = vec3(sin(fi) + 2.*sin(2.*fi), cos(fi) - 2.*cos(2.*fi), z);
        d = min(length(p - vec3(v.xy/1., v.z)), d);

        fi = PI - f + TAU/3.*i;
        v = vec3(sin(fi) + 2.*sin(2.*fi), cos(fi) - 2.*cos(2.*fi), z);
        d = min(length(p - vec3(v.xy/1., v.z)), d);
    }
    d =d*0.3  - 0.05;
    return d;
}

int quadratic(float A, float B, float C, out vec2 x) {
   float D = B*B - 4.0*A*C;
   if (D < 0.0) return 0;
   D = sqrt(D);
   if (B < 0.0) D = -D;
   x[0] = (-B-D)/(2.0*A);
   x[1] = C/(A*x[0]);
   return 2;
}

float knot2(vec3 p, float shift)
{
    float dz = 1.; //node thickness
    p.xy *= rot(TAU/3.0*shift);

    float d  = dist_infin;
    float y = clamp(p.y, -3., 2. + 1./16.);

    vec2 x = vec2(dist_infin);
    int n = quadratic(-4., 1., 2. - y, x);
    
    if (n == 0)
    {
        x[0] = 0.;
        x[1] = 0.;
    }
    for (int i = 0; i <2; i++)
    {
        
        float f = acos(x[i]);
        float tx = sin(f) + 2.*sin(2.*f);
        d = min(d, length(vec3(p.x - tx, p.y - y, p.z + dz*sin(3.*f))));
        d = min(d, length(vec3(p.x + tx, p.y - y, p.z - dz*sin(3.*f))));
    }
    return d;
}

float knot3(vec3 p)
{
    float d = dist_infin;
    for (float i = 0.; i < 3.; i++)
        d = min(d, knot2(p, i));
    return d*0.45 - 0.1;    
}
/*
cos(f) - 2(cos^2 - 1 + cos^2) = cos(f) - 4cos^2(f) + 2
sin(f) + 2sin(2f) = sin(f) + 4(sin(f))cos(f)
sin^2 + 8 sin^2fcosf + 16sin^2cos^2f
*/


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

    //float d = knot(p);
    float d = sdCosNp(p.xy, 3., 2.);
    d = length(vec2(d, p.z)) * 0.45 - 0.02;

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
//https://mathcurve.com/courbes2d/larme/larme.shtml
float larme2(vec2 p, float a)
{
    float k = 0.3;
    p.x += a*(k - 0.95)/2.0;
    float zoom = 1.;
    float x = clamp(p.x, -a*0.95, a*k);
    float f = acos(x/a);
    float y = a*sin(f)*pow(sin(f/2.), 2.)*zoom;
    float d = length(p - vec2(x, y));
    //d = min(length(p - vec2(x, -y)), d);
    if (abs(p.y) < abs(y))
    {
        sdfReflect = 0.;
        sdfColor = col2;
    }
    else
    {
        sdfReflect = 0.1;
        sdfColor = col1;
    }
    
    f = acos(k);
    y = a*sin(f)*pow(sin(f/2.), 2.)*zoom;
    float y2 = clamp(p.y, -y, y);
    float d2 = length(p - vec2(k*a, y2));
    if (d2 < d)
    {
        d = d2;
        if (p.x > x)
        {
            sdfReflect = 0.1;
            sdfColor = col1;
        }
        else
        {
            sdfReflect = 0.;
            sdfColor = col2;
        }
    }
    
    
    return d;
}

float larme(vec3 p, float a)
{
    //float d = larme2(p.xy, a);
    //d = length(vec2(p.z, d));
    float l = length(p.xy);
    float d = larme2(vec2(p.z, l), a);
    return d*0.3 - 0.02;
}

//https://mathcurve.com/courbes2d.gb/bouche/bouche.shtml
float kiss2(vec2 p, float a)
{
    
    float x = clamp(p.x, -a, a);
    float f = acos(x/a);
    float y = a*pow(sin(f),3.)/3.;
    float d = length(p - vec2(x, y));
    d = min(length(p - vec2(x, -y)), d);
    return d;
}
float kiss(vec3 p, float a)
{
    //float d = kiss2(p.xy, a);
    //d = length(vec2(p.z, d));
    float l = length(p.xy);
    float d = kiss2(vec2(l, p.z), a);
    return d*0.5 - 0.03;
}
//https://mathcurve.com/courbes2d/cissoiddroite/cissoiddroite.shtml
float ciss(vec3 p, float a)
{
    float f = atan(p.y, p.x);
    float r = a*(1./cos(f) - cos(f));
    float x = r*cos(f);
    float y = r*sin(f);
    float d = length(p - vec3(x, y, 0.));
    return d*0.5 - 0.03;
}

float map(vec3 p) {

    float d = larme(p, 3.5);//knot3(p);//kiss(p, 2.);//ciss(p, 2.);//
    resColor = sdfColor;
    resReflect = sdfReflect;
    return d;
}


vec3 csky(vec3 p) {
    float n = 6., m = 6., dlat = PI / n, dlon = TAU / m;
    float lon = mod(atan(p.y, p.x), TAU), lat = atan(length(p.xy), p.z);
    float fo = fract(lon / dlon), fa = fract(lat / dlat);

    float pst = fo * fa * (1. - fo) * (1. - fa);
    pst=pow(pst, 2.);

    pst = smoothstep(0.0, pow(0.0625, 2.), pst);
    pst = clamp(pst, 0.1, 1.0);
    return vec3(pst);
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

void mainImage(out vec4 fragColor, in vec2 fragCoord) {
    vec3 light = normalize(vec3(0.0, 1.0, -2.5)); //light
    vec3 light2 = normalize(vec3(0.0, -1.0, 2.5)); //light
    vec2 mo = vec2(0.0, 0.0);
    //if  (iMouse.z > 0.0)
    {
        mo = (-iResolution.xy + 2.0 * (iMouse.xy)) / iResolution.y;
    }
    vec3 ro = vec3(0.0, 0.0, 6.0); // camera
    //camera rotation
    ro.yz *= rot(mo.y * PI);
    ro.xz *= rot(-mo.x * TAU);

    const float fl = 1.5; // focal length
    float dist = dist_infin;

    vec3 b1 = vec3(0.23529411764705882, 0.4235294117647059, 0.7725490196078432), b2 = vec3(0.3686274509803922, 0.5725490196078431, 0.8941176470588236);
    vec3 bg = mix(b2, b1*b1, fragCoord.y / iResolution.y);   
    //vec3 bg = vec3(0.);
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
    tot = tot / float(AA);
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