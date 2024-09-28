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
pumpkin
sdf,raymarching,noise,pumpkin,halloween
pumpkin of halloween
*/
float hash(float n) {
    return fract(sin(n) * 437558.5453123);
}

float hash(vec2 p) {
	return fract(1e4 * sin(17.0 * p.x + p.y * 0.1) * (0.1 + abs(sin(p.y * 13.0 + p.x))));
} 
//	<https://www.shadertoy.com/view/4dS3Wd>
//	By Morgan McGuire @morgan3d, http://graphicscodex.com



// This one has non-ideal tiling properties that I'm still tuning
float noise(vec3 x) {
    const vec3 step = vec3(110, 241, 171);

    vec3 i = floor(x);
    vec3 f = fract(x);

	// For performance, compute the base input to a 1D hash from the integer part of the argument and the 
	// incremental change to the 1D based on the 3D -> 1D wrapping
    float n = dot(i, step);

    vec3 u = f * f * (3.0 - 2.0 * f);
    return mix(mix(mix(hash(n + dot(step, vec3(0, 0, 0))), hash(n + dot(step, vec3(1, 0, 0))), u.x), mix(hash(n + dot(step, vec3(0, 1, 0))), hash(n + dot(step, vec3(1, 1, 0))), u.x), u.y), mix(mix(hash(n + dot(step, vec3(0, 0, 1))), hash(n + dot(step, vec3(1, 0, 1))), u.x), mix(hash(n + dot(step, vec3(0, 1, 1))), hash(n + dot(step, vec3(1, 1, 1))), u.x), u.y), u.z);
}

float hash(vec3 p) {
    return fract(sin(dot(p, vec3(127.1, 311.7, 74.7))) * 43758.5453123);
}

#define PI  3.14159265359
#define TAU 6.28318530718
#define rot(f) mat2(cos(f), -sin(f), sin(f), cos(f))

const float dist_infin = 10.0;
#define nn 128
const float eps = 0.001;

vec3 getSg(vec3 p, float nseg) {
    float fi = mod(atan(p.y, p.x), TAU);
    fi = mod(fi + PI / nseg, TAU);
    float n = floor(fi / TAU * nseg);
    p.xy *= rot(-n * TAU / nseg);
    return p;
}

float sdBox(vec3 p, vec3 b) {
    vec3 q = abs(p) - b;
    return length(max(q, 0.0)) + min(max(q.x, max(q.y, q.z)), 0.0);
}

float smin_out(float d1, float d2, float k) {
    float h = clamp(0.5 - 0.5 * (d2 + d1) / k, 0.0, 1.0);
    return mix(d2, -d1, h) + k * h * (1.0 - h);
}


float halloween(vec3 p) {
    p.yz *= rot(PI / 2.);
    //p.xy *= rot(iTime*0.3);
    p.xz *= rot(PI / 8.);
    float n1 = 14., R = 1.5;
    float d0 = length(p) - 1.45;

    float fi1 = TAU / n1 + PI / n1, li = PI / 4., r = R * cos(li);//PI/n1
    vec3 e1p = vec3(r * cos(fi1), r * sin(fi1), R * sin(li));
    float e1 = length(p - e1p) - .45;
    //fi1+= 3.*TAU/n1;
    fi1 = -fi1;
    e1p = vec3(r * cos(fi1), r * sin(fi1), R * sin(li));
    float e2 = length(p - e1p) - .45;

    fi1 = 0., li = -PI / 10.;
    vec3 lip = vec3(r * cos(fi1), r * sin(fi1), R * sin(li));
    vec3 lv = p - lip;

    float lips = sdBox(lv, vec3(0.9, 1., 0.2));

    p = getSg(p, n1);
    p.xy *= rot(-PI / n1 * sign(p.y));
    float d = abs(length(p.xz) - R);
    d = length(vec2(d, p.y)) - 0.35;

    d = max(d, -d0);
    d = smin_out(e1, d, 0.1);
    d = smin_out(e2, d, 0.1);
    d = smin_out(lips, d, 0.2);
    
    return d;
}

float map(vec3 p) {
    return halloween(p);
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
//https://www.shadertoy.com/view/4td3zj
// Compact, self-contained version of IQ's 3D value noise function.
float n3D(vec3 p){
    
	const vec3 s = vec3(7, 157, 113);
	vec3 ip = floor(p); p -= ip; 
    vec4 h = vec4(0., s.yz, s.y + s.z) + dot(ip, s);
    p = p*p*(3. - 2.*p); //p *= p*p*(p*(p * 6. - 15.) + 10.);
    h = mix(fract(sin(mod(h, 6.2831589))*43758.5453), 
            fract(sin(mod(h + s.x, 6.2831589))*43758.5453), p.x);
    h.xy = mix(h.xz, h.yw, p.y);
    return mix(h.x, h.y, p.z); // Range: [0, 1].
}

// Simple environment mapping. Pass the reflected vector in and create some
// colored noise with it. The normal is redundant here, but it can be used
// to pass into a 3D texture mapping function to produce some interesting
// environmental reflections.
vec3 envMap(vec3 rd, vec3 sn){
    
    vec3 sRd = rd; // Save rd, just for some mixing at the end.
    
    // Add a time component, scale, then pass into the noise function.
    rd.xy -= iTime*.25;
    rd *= 3.;
    
    float c = n3D(rd)*.57 + n3D(rd*2.)*.28 + n3D(rd*4.)*.15; // Noise value.
    c = smoothstep(.4, 1., c); // Darken and add contast for more of a spotlight look.
    
    vec3 col = vec3(c, c*c, c*c*c*c); // Simple, warm coloring.
    //vec3 col = vec3(min(c*1.5, 1.), pow(c, 2.5), pow(c, 12.)); // More color.
    
    // Mix in some more red to tone it down and return.
    return mix(col, col.yzx, sRd*.25+.25); 
    
}
/*
#if HW_PERFORMANCE==0
#define AA 1
#else
#define AA 2
#endif
*/
#define AA 2
vec3 col0 = vec3(0.776, 0.525, 0.301);
vec3 col2 = vec3(0.72, 0.01, 0.01);
float npp = 75.;
float lev = 0.995;

void mainImage(out vec4 fragColor, in vec2 fragCoord) {
    vec3 light = normalize(vec3(0.0, .0, 1.)); //light
    vec3 light2 = normalize(vec3(0.0, 0.0, -1.)); //light
    //vec2 mo = 1.5 * cos(0.3 * iTime + vec2(0, 11));
    vec2 mo = vec2( -0.3 * iTime, 0.);
    //if  (iMouse.z > 0.0)
    {
        mo = (-iResolution.xy + 2.0 * (iMouse.xy)) / iResolution.y;
        mo *= 1.6;
    }
    vec3 ro = vec3(0.0, 0.0, 5.); // camera
    //camera rotation
    ro.yz *= rot(mo.y);
    ro.xz *= rot(-mo.x - 1.7);

    const float fl = 1.5; // focal length
    float dist = dist_infin;

    vec3 b1 = vec3(0.23529411764705882, 0.4235294117647059, 0.7725490196078432), b2 = vec3(0.3686274509803922, 0.5725490196078431, 0.8941176470588236);
    vec3 bg = mix(b2, b1, fragCoord.y / iResolution.y);  
    //antialiasing
    vec3 tot = vec3(0.0);
    for(int m = 0; m < AA; m++) for(int n = 0; n < AA; n++) {
            vec2 o = vec2(float(m), float(n)) / float(AA) - 0.5;
            vec2 p = (-iResolution.xy + 2.0 * (fragCoord + o)) / iResolution.y;
            vec3 rd = GetRayDir(p, ro, vec3(0, 0., 0), fl); //ray direction
            vec3 col = bg * bg; // background  

            //================================sky color========================
            
            vec2 pp = floor(p * npp) / npp;
            float fil = hash(pp);
            if(fil > lev) {
                if((length(p - (pp + vec2(0.5 / npp, 0.5 / npp)))) < 0.5 / npp)
                    col = vec3(1.);
            }

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

                vec3 nor = calcNormal(pos);
                if(dot(pos, nor) < 0.0) {
                    col = col2;
                    float t = noise(pos * 4.);
	                t = fract(t * 2.);
                    t = smoothstep(0., 0.01, t - 0.3) * smoothstep(0., 0.01, 0.5 - t);
                    col *= 1.-t;
                    col*=envMap(rd, nor)*6.;
                } else {
                    col = col0 * col0;
                    float t = noise(pos);
                    t = smoothstep(0., 0.01, t - 0.3);
                    col *= t;
                }

                vec3 R = reflect(light, nor);
                float specular = pow(max(abs(dot(R, rd)), 0.), 16.);
                float difu = abs(dot(nor, light));
                col = col * (clamp(difu, 0., 1.0) + 0.5) + vec3(1., .7, .4) * specular;
                float fre = pow(clamp(dot(nor, rd) + 1., .0, 1.), 3.); // Fresnel, for some mild glow.
                col += vec3(.1, .1, 0.1) * fre; //?
                col = sqrt(col);
            }
            //==========================raymatch=============================
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