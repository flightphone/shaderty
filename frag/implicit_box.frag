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
Exact Schwarz P surface

raytracing,implicit,schwarz,bisect

In this shader, an exact surface is built, using the technology of displaying an implicit surface, not SDF. Just like in here:
[url]https://www.shadertoy.com/view/DlXXR2[/url]
[url]https://www.shadertoy.com/view/lfSBzt[/url]
[url]https://www.shadertoy.com/view/Mf2fWV[/url]
*/
#define PI  3.14159265359
#define TAU 6.28318530718
#define rot(f) mat2(cos(f), -sin(f), sin(f), cos(f))
#define nn 64.
#define newton 5

float dist_inf = 6.;
const float eps = 0.001;

float sch(vec3 p) {
    p *= 5. * PI;
    return cos(p.x) + cos(p.y) + cos(p.z);
}

float lines(vec3 p) {
    //return dist_inf;
    if(abs(p.x) > 0.5 || abs(p.y) > .5)
        return dist_inf;
    float h = 1.1, r = 0.125, n = 5.;
    p.xy = (p.xy + 1.) / 2.;
    p.xy = fract(p.xy * n) - 0.5;
    p.z -= clamp(p.z, -h, h);
    return length(p) - r;
}

vec3 colx = vec3(0.7, 0. , 0.);
vec3 coly = vec3(0., 0.7 , 0.);
vec3 colz = vec3(0., 0., 0.7);
vec3 col_axis = vec3(0.);
//==================implicit map=====================
float map(vec3 p, int model) {
    if(model == 0)
        return sch(p);
    if(model == 1) {
        
        float t2 = lines(p);
        float t = t2;
        col_axis = colx;
        
        float t3 = lines(p.yzx);
        if (t3 < t)
        {
            t = t3;
            col_axis = coly;
        }
        float t4 = lines(p.xzy);
        if (t4 < t)
        {
            t = t4;
            col_axis = colz;
        }
        return t;
    }
}

vec3 calcNormal(in vec3 p, int model) {
    const float eps = 0.0001;
    vec2 q = vec2(0.0, eps);
    vec3 res = vec3(map(p + q.yxx, model) - map(p - q.yxx, model), map(p + q.xyx, model) - map(p - q.xyx, model), map(p + q.xxy, model) - map(p - q.xxy, model));
    return normalize(res);
}

//https://iquilezles.org/articles/intersectors/
// axis aligned box centered at the origin, with size boxSize
vec2 boxIntersection(in vec3 ro, in vec3 rd, vec3 boxSize) {
    vec3 m = 1.0 / rd; // can precompute if traversing a set of aligned boxes
    vec3 n = m * ro;   // can precompute if traversing a set of aligned boxes
    vec3 k = abs(m) * boxSize;
    vec3 t1 = -n - k;
    vec3 t2 = -n + k;
    float tN = max(max(t1.x, t1.y), t1.z);
    float tF = min(min(t2.x, t2.y), t2.z);
    if(tN > tF || tF < 0.0)
        return vec2(dist_inf); // no intersection

    return vec2(tN, tF);
}

vec3 getPoint(vec3 a, vec3 b, float v0, float v1, int model) {
    vec3 m;
    //binary search with  n iterations, n = newton
    for(int i = 0; i < newton; i++) {
        m = (a + b) * 0.5;
        float v = map(m, model);
        if(v == 0.)
            break;

        if(sign(v) * sign(v0) <= 0.) {
            v1 = v;
            b = m;
        } else {
            v0 = v;
            a = m;
        }
    }
    return m;
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
    float hh = 3.;

    vec3 light = normalize(vec3(0.0, 1.0, -2.5)); //light
    vec2 mo = 1.5 * cos(0.5 * iTime + vec2(0, 11));
    //if  (iMouse.z > 0.0)
    {
        mo = (-iResolution.xy + 2.0 * (iMouse.xy)) / iResolution.y;
    }
    vec3 ro = vec3(0.0, 0.0, hh); // camera
    //camera rotation
    ro.yz *= rot(mo.y * 2.);
    ro.xz *= rot(-mo.x - 1.57);

    const float fl = 1.5; // focal length

    vec3 bg = vec3(0.7, 0.7, 0.9) * 0.6; 

    //antialiasing
    vec3 tot = vec3(0.0);
    for(int m = 0; m < AA; m++) for(int n = 0; n < AA; n++) {
            vec2 o = vec2(float(m), float(n)) / float(AA) - 0.5;
            vec2 p = (-iResolution.xy + 2.0 * (fragCoord + o)) / iResolution.y;
            vec3 rd = GetRayDir(p, ro, vec3(0, 0., 0), fl); //ray direction
            vec3 col = bg; // background  

            vec3 pos_gl = vec3(1000.);
            //===========================Render implicit=====================
            //STEP 1. Calculating bounding box
            int flInt = -1;
            float dist = 2.0;

            for(int model = 0; model < 2; model++) {
                vec2 d2 = boxIntersection(ro, rd, vec3(dist * (float(model)/4. + 1.)));
                if(d2[0] < dist_inf) {

            /*
            STEP 2.
            ray tracing inside the bounding box, 
            searching for a segment with different signs of the function value 
            at the ends of the segment
            */

                    vec3 pos0 = ro + rd * d2[0];
                    vec3 pos1 = ro + rd * d2[1];
                    vec3 rd0 = pos1 - pos0;
                    vec3 pos = pos0;
                    float val0 = map(pos0, model);
                    for(float i = 1.; i < nn; i++) {
                        pos1 = pos0 + rd0 * i / (nn - 1.);
                        float val1 = map(pos1, model);
                        if(sign(val0) * sign(val1) <= 0.) {
                    //different signs of the function value  at the ends of the segment
                    //STEP 3. binary search to clarify the intersection of a ray with a surface.
                            pos = getPoint(pos, pos1, val0, val1, model);
                            if(length(pos - ro) < length(pos_gl - ro))
                            {
                                flInt = model;
                                pos_gl = pos;
                            }
                            break;
                        }
                //if (sign(val1) < 0.) col = col2;
                        val0 = val1;
                        pos = pos1;
                    }
                }
            }

            //===================color====================================
            if(flInt > -1) {
                vec3 col1 = vec3(139. / 255., 0., 1.);
                vec3 col2 = vec3(0.8, 0.5, 0.3);
                vec3 col3 = vec3(1., 0., 0.);
                col1 *= col1;
                col2 *= col2;


                vec3 nor = calcNormal(pos_gl, flInt);
                if(flInt == 0)
                    col = col2;
                else
                    col = col_axis;

                if (dot(rd, nor) < 0.0 && flInt == 0)
                    col *= col;

            //  float ao = calculateAO(pos, nor); // Ambient occlusion, for self shadowing. 

                vec3 R = reflect(light, nor);
                float specular = pow(max(abs(dot(R, rd)), 0.), 16.);
                float difu = abs(dot(nor, light));
                col = col * (clamp(difu, 0., 1.0) + 0.5) + vec3(1., .7, .4) * specular;
                float fre = pow(clamp(dot(nor, rd) + 1., .0, 1.), 3.); // Fresnel, for some mild glow.
                col += vec3(.1, .1, 0.1) * fre; //?
                col = sqrt(col);
            }
            //===================color====================================
            tot += col;
        }
    tot = tot / float(AA) / float(AA);
    //antialiasing
    fragColor = vec4(tot, 1.0);
}
/////=====================================================================================
void main() {
    vec4 fragColor = vec4(0);
    mainImage(fragColor, gl_FragCoord.xy);
    gl_FragColor = fragColor;
}