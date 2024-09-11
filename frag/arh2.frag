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
golden tower
SDF,raymatch,noise,arch,architecture
SDF for simple arch. Noise for simple starry sky
*/


#define PI  3.14159265359
#define TAU 6.28318530718
#define rot(f) mat2(cos(f), -sin(f), sin(f), cos(f))

const float dist_infin = 10.0;
#define nn 128

const float eps = 0.001;
vec3 bg = vec3(0.08, 0.42, 0.87);
vec3 col1 = vec3(0.73, 0.7, 0.4);
float npp =15.;
float lev = 0.995;

float hash (vec3 p) {
    return fract(sin(dot(p, vec3(127.1,311.7, 74.7))) * 43758.5453123);
}
float arch(vec3 p, float R, float h, float l) {
    float res = 0.;
    if(p.z >= 0.)
        res = length(vec2(p.x, max(p.z - h, 0.))) - R;
    else
        res = length(vec2(p.z, abs(p.x) - R));
    res = length(vec2(max(abs(p.y) - l, 0.), res));
    return res;
}

float level(vec3 p, float R, float h)
{
    return max((arch(p, R, h, 100.)), (arch(p.yxz, R, h, 100.)));
}

float sdBox( vec3 p, vec3 b )
{
  vec3 q = vec3(0.);  
  abs(-p.z), (p.z - b.z);
  if (p.z > 0.)  
    q = vec3(p.z - b.z, abs(p.xy)-b.xy);
  else
    q = vec3(-p.z, abs(p.xy)-b.xy);  

  return length(max(q,0.0)) + min(max(q.x,max(q.y,q.z)),0.0);
}

float arch2(vec3 p, float R, float h)
{
    float res = 0.;
    res = min(length(vec2(p.x, max(p.z - h, 0.))), length(vec2(p.y, max(p.z - h, 0.)))) - R;
    return res;    
}

float level2 (vec3 p, float R, float h, float w, float H)
{
    float t = sdBox(p, vec3(w, w, H));
    float t2 = arch2(p, R, h);
    return max (t, -t2);
}


float map(vec3 p) {
    p.yz *= rot(PI/2.);
    p.z -= 1.;
    float t0 = level(p, 0.2, 0.3) - 0.1;
    p.z += 0.95;
    float t1 = level(p, 0.4, 0.5) - 0.15;
    p.z += 1.55;
    float t2 = level2(p, 0.4, 0.9, 0.65, 1.5)-0.03;
    return min(t0, min(t1, t2));
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
void mainImage(out vec4 fragColor, in vec2 fragCoord) {
    vec3 light = normalize(vec3(0.0, 1.0, -2.5)); 
    vec2 mo = vec2( -0.2 * iTime, 0.);
    //if  (iMouse.z > 0.0)
    {
        mo = (-iResolution.xy + 2.0 * (iMouse.xy)) / iResolution.y;
    }
    vec3 ro = vec3(0., 0., 3.4); // camera
    //camera rotation
    ro.yz *= rot(mo.y*2.);
    ro.xz *= rot(-mo.x*2.);
   
    const float fl = 1.5; // focal length

    //antialiasing
    vec3 tot = vec3(0.0);
    for(int m = 0; m < AA; m++) for(int n = 0; n < AA; n++) {
            vec2 o = vec2(float(m), float(n)) / float(AA) - 0.5;
            vec2 p = (-iResolution.xy + 2.0 * (fragCoord + o)) / iResolution.y;
            vec3 rd = GetRayDir(p, ro, vec3(0, 0., 0), fl); //ray direction
            vec3 col = bg*bg;
            //================================sky color========================
            float d = length(cross(ro, rd));
            float td = abs(dot(ro, rd));
            d = sqrt(dist_infin * dist_infin - d * d);
            vec3 pos = ro + rd * (td + d);
            vec3 pp = floor(pos*npp)/npp;;
            float fil = hash(pp);
            if(fil > lev) {
                if ((length(pos - (pp + vec3(0.5/npp, 0.5/npp, 0.5/npp)))) < 0.5/npp)
                    col = vec3(1.);
            }
            //==========================raymatch=============================
            td = 0.;
            pos = vec3(0.);
            for(int i = 0; i < nn; i++) {
                pos = ro + rd * td;
                float h = map(pos);
                if(h < eps || td >= dist_infin)
                    break;
                td += h;
            }
            //======================color====================================
            if(td < dist_infin) {
                col = col1*col1;
                vec3 nor = calcNormal(pos);
                vec3 R = reflect(light, nor);
                float specular = pow(max(abs(dot(R, rd)), 0.), 16.);
                float difu = abs(dot(nor, light));
                col = col * (clamp(difu, 0., 1.0) + 0.5) + vec3(1., .7, .4) * specular;
                float fre = pow(clamp(dot(nor, rd) + 1., .0, 1.), 3.); // Fresnel, for some mild glow.
                col += vec3(.4, .7, 1) * fre; //?
                col = sqrt(col);
            }
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