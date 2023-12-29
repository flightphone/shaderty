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


float map(vec3 p) {
    p.xy *=rot(iTime);
    float pst = smoothstep(0.1, 0., abs(atan(p.y, p.x) - PI/4.)); 
    resColor = mix(col1, col2, pst);
    return (length(p) - 1.);
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

vec3 calccolor(vec3 col, vec3 backcol, vec3 rd, vec3 light1, vec3 light2, vec3 nor) {
    
    float difu1 = dot(nor, light1);
    float difu2 = dot(nor, light2);
    float difu = max(difu1, difu2);
    

    vec3 R1 = reflect (light1, nor);
    vec3 R2 = reflect (light2, nor);
    float shininess=10.0;
    float specular1    =  pow(max(dot(R1, rd), 0.), shininess);
    float specular2    =  pow(max(dot(R2, rd), 0.), shininess);
    float specular = max(specular1, specular2);

    
    //col = col*clamp(difu, 0.3, 1.0) + vec3(.5)*specular*specular;
    col = col*(col*clamp(difu, 0., 1.0) + 0.3) + vec3(.5)*specular*specular;

    return col;
}

vec3 ccolor(vec3 col, vec3 rd, vec3 light, vec3 nor) {
    
    //jorge2017a2
    
    vec3 R = reflect (light, nor);
    float specular    =  pow(max(dot(R, rd), 0.), 16.);
    
    //iq
    /*
    vec3  hal = normalize( light-rd );
    float specular = pow( clamp( dot( nor, hal ), 0.0, 1.0 ),16.0);
    */

    float difu = dot(nor, light);
    //col = col*clamp(difu, 0.3, 1.0) + vec3(.5)*specular*specular;
    col = col*(col*clamp(difu, 0., 1.0) + 0.3) + vec3(.5)*specular*specular;
    return col;
}

void mainImage(out vec4 fragColor, in vec2 fragCoord) {
    vec3 light = normalize(vec3(0.0, .0, 1.)); //light
    vec3 light2 = normalize(vec3(0.0, 0.0, -1.)); //light
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
                //col = calccolor(col, col, rd, light, light2, nor);
                col = ccolor(col, rd, light, nor);
            }
            //==========================raymatch=============================
            tot += col;
        }
    tot = sqrt(tot) / float(AA);
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