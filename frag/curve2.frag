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

const float dist_infin = 10.0;
#define nn 128
const float eps = 0.001;


float sdCosN(vec3 p, float a, float n) {

    float fi = atan(p.y, p.x);
    float L = length(p.xy);
    float d1 = dist_infin;
    float d2 = dist_infin;

    for(float i = 0.; i < 10.; i++) {
        if(fract(n * i) == 0. && i > 0.)
            break;
        float r = a * cos(n * fi + n * i * TAU);
        d1 = min(abs(L - r), d1);
        //d1 = min(L-r, d1);
    }
    float f = acos(L / a) / n;
    for(float i = 0.; i < 10.; i++) {
        
        d2 = min(2.0 * abs(sin((fi - (f + i * TAU / n)) / 2.0)) * L, d2);
        d2 = min(2.0 * abs(sin((fi + (f + i * TAU / n)) / 2.0)) * L, d2);
        d2 = min(2.0 * abs(sin((fi - (f - i * TAU / n)) / 2.0)) * L, d2);
        d2 = min(2.0 * abs(sin((fi + (f - i * TAU / n)) / 2.0)) * L, d2);
        
        
    }
    //float d = smin2(d1, d2);
    float e = 1.;
    float d = min(d1*e, d2);
    //float d = d2;
    if (d < 0.)
        d = 0.;
    d = length(vec2(p.z*e, d));
    d *= .6;
    d -= 0.03;
    return d;
}

mat3 rotateX(float f) {
    return mat3(vec3(1.0, 0.0, 0.0), vec3(0.0, cos(f), -sin(f)), vec3(.0, sin(f), cos(f)));
}

mat3 rotateZ(float f) {
    return mat3(vec3(cos(f), -sin(f), 0.0), vec3(sin(f), cos(f), 0.0), vec3(0.0, 0.0, 1.0));

}

mat3 rotateY(float f) {
    return mat3(vec3(cos(f), 0.0, sin(f)), vec3(0.0, 1.0, 0.0), vec3(-sin(f), 0.0, cos(f)));
}


float sdRound(vec3 p, float r, float f) {
    float d = abs(length(p.xy) - r * cos(f));
    d = length(vec2(p.z - r * sin(f), d));
    return d;
}
float f1(float x)
{
    return x*x*2.;
}

float sdff1(vec2 p, float a, float b)
{
    float x = clamp(p.x, a, b);
    float v = f1(x);
    return length(p - vec2(x, v));
}

float f2(float x)
{
    return sqrt(x/2.);
}

float sdff2(vec2 p, float a, float b)
{
    float x = clamp(p.x, a, b);
    float v = f2(x);
    return length(p - vec2(x, v));
}


float sdf_fun(vec3 p)
{
    float d1 = sdff1(p.xy, 0., .8);
    float d2 = sdff2(p.yx, 0., 1.28);
    float d =  min(d1, d2);
    d = length(vec2(p.z, d));
    d *= 0.5;
    d -= 0.06;
    return d;
}

// FIXME need better hair SDF
float dCyl(vec3 p, float r, float h)
{
	r = max(.0, length(p.xy) - r);
    if (p.z < 0.)
        h = p.z;
    else
        h = max(0., p.z - h);
	return length(vec2(r, h));
}

float sdCapsule( vec3 p, vec3 a, vec3 b, float r )
{
  vec3 pa = p - a, ba = b - a;
  float h = clamp( dot(pa,ba)/dot(ba,ba), 0.0, 1.0 );
  return length( pa - ba*h ) - r;
}

float heir(vec3 p)
{
    float n = 20., 
          m = 20.,
          h = .1, //height pimple
          r = 0.01; // radius pimple  

    vec3 ell = vec3(.6, 1., 1.2); //ellipsoid
    //vec3 ell = vec3(1., 1., 1.); //sphere
    float k0 = length(p/ell),
          k1 = length(p/(ell*ell)),
          dlon = TAU/n, dlat = PI/m, l = length(p.xy), 
          //https://iquilezles.org/articles/distfunctions/ 
          dz = k0*(k0-1.0)/k1, //distance to ellipsoid
          lon = mod(atan(p.y,p.x), TAU),
          lat = atan(l, p.z), //longitude and latitude
          lon1 = floor(lon/dlon)*dlon + 0.5*dlon,
          lat1 = floor(lat/dlat)*dlat + 0.5*dlat; //longitude and latitude nearest pimple
    
    float x1 = ell.x * sin(lat1)*cos(lon1), 
          y1 = ell.y*sin(lat1)*sin(lon1), 
          z1 = ell.z * cos(lat1), //coordinate nearest pimple
          x = ell.x * sin(lat)*cos(lon), 
          y = ell.y *sin(lat)*sin(lon), 
          z = ell.z * cos(lat), //coordinate point on ellipsoid
          dxy = length(vec3(x, y, z) - vec3(x1, y1, z1)),
          dp = (length(vec2(dxy, dz - clamp(dz, 0., h)))*0.5 - r); //dist to pimple  
   
    return min(dz, dp) ;
}


float map(in vec3 pos) {
    
    
    //return sdCosN(pos, 1.0, .3);
    //return sdf_fun(pos);
    //return dCyl(pos, 0.2, .5);
    return heir(pos);


    
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
    vec3 f = normalize(l - p), r = normalize(cross(vec3(0, 1, 0), f)), u = cross(f, r), c = f * z, i = c + uv.x * r + uv.y * u;
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
void mainImage(out vec4 fragColor, in vec2 fragCoord) {
    vec3 light = normalize(vec3(0.0, 1.0, 1.0)); //light
    vec3 light2 = normalize(vec3(0.0, 1.0, -1.0)); //light
    float t = iTime / 3.;
    vec2 m = vec2(0.0, 0.0);
    //if  (iMouse.z > 0.0)
    {
        m = (-iResolution.xy + 2.0 * (iMouse.xy)) / iResolution.y;
        t = 0.;
    }
    vec3 ro = vec3(0.0, 0.0, 2.5); // camera
    ro = rotateY(-m.x * TAU) * rotateX(-m.y * PI) * ro; //camera rotation

    const float fl = 1.5; // focal length
    float dist = dist_infin;
    mat3 rota = rotateZ(t) * rotateX(-t);
    mat3 rota_1 = rotateX(t) * rotateZ(-t);

    vec3 tot = vec3(0.0);

    //antiblick
    for(int m = 0; m < AA; m++) for(int n = 0; n < AA; n++) {
            vec2 o = vec2(float(m), float(n)) / float(AA) - 0.5;
            vec2 p = (-iResolution.xy + 2.0 * (fragCoord + o)) / iResolution.y;
            vec3 rd = GetRayDir(p, ro, vec3(0, 0., 0), fl); //ray direction
            vec3 col = vec3(0.0);

            HIT giper = giper3D(rota * ro, rota * rd);
            if(giper.dist < dist) {
                vec3 nor = rota_1 * giper.nor;
                float dif = clamp(dot(nor, light), 0.2, 1.0);
                float amb = 0.5 + 0.5 * dot(nor, light2);
                col = vec3(0.2, 0.3, 0.4) * amb + vec3(0.85, 0.75, 0.65) * dif;
            }
        // gamma        
            col = sqrt(col);
            tot += col;
        }
    //antiblick
    tot /= float(AA * AA);
    fragColor = vec4(tot, 1.0);
}
/////=====================================================================================
void main() {
    vec4 fragColor = vec4(0);
    mainImage(fragColor, gl_FragCoord.xy);
    gl_FragColor = fragColor;
}