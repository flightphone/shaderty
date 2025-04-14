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
amulet for good luck
SDF,raymarching,knot,endless,tibet,buddhism,amulet 
Fork [url]https://www.shadertoy.com/view/W3SGDK[/url], [url]https://www.shadertoy.com/view/w3SGRm[/url]
*/
#define PI  3.14159265359
#define TAU 6.28318530718
#define rot(f) mat2(cos(f), -sin(f), sin(f), cos(f))

const float dist_infin = 10.0;
#define nn 128
const float eps = 0.001;
vec3 colf = vec3(222. / 255., 208. / 255., 159. / 255.);

float sdSegment(in vec2 p, in vec2 a, in vec2 b) {
    vec2 pa = p - a, ba = b - a;
    float h = clamp(dot(pa, ba) / dot(ba, ba), 0.0, 1.0);
    return length(pa - ba * h);
}

vec3 curColorN(float d1, vec3 col, vec3 col1, float dd) {
    float s1 = smoothstep(0., -5. / iResolution.y, d1);
    float cs = dd * dd - (dd - abs(d1)) * (dd - abs(d1));
    if(cs < 0.) {
        cs = 1.;
        d1 = 0.;
    } else {
        cs = sqrt(cs) / dd;
        d1 = d1 / dd;
    }
    vec3 norm = vec3(d1, 0., cs);
    vec3 light = vec3(0., 0., 1.);
    vec3 rd = vec3(0., 0., -1.);
    float difu = dot(norm, light);
    vec3 R1 = reflect(light, norm);
    float shininess = 15.0;
    float specular = pow(max(dot(R1, rd), 0.), shininess);
    vec3 colccs = col1 * difu + 0.8 * specular;
    return mix(col, colccs, vec3(s1));
}

vec3 draw_cyrcle(vec2 p, float r, float dd, vec3 col, vec3 col1) {
    float d1 = abs(length(p) - r) - dd;
    return curColorN(d1, col, col1, dd);

}

float bridge(vec2 p0, float n, float dd) {
    vec2 c = vec2(0.5 / n, (1. + dd) / n);
    float d = abs(length(p0 - c) - 0.5 / n) * n;
    return d;
}

float bridge2(vec2 p0, float n, float dd) {
    float d = sdSegment(p0, vec2(-1. / n, 1. / n), vec2(-1. / n, 1.5 / n)) * n;
    d = min(sdSegment(p0, vec2(-1. / n, 1. / n), vec2(-1.5 / n, 1. / n)) * n, d);

    vec2 c = p0 - vec2(-1.5 / n, 1.5 / n);
    float a = mod(atan(c.y, c.x), TAU);
    if(a <= 1.5 * PI)
        d = min(d, abs(length(c) - 0.5 / n) * n);

    return d;

}

vec3 endless(vec2 p0) {
    float n = 3.1, d = 10.;
    vec3 col1 = vec3(0.3, 0.3, 1.);
    ;
    vec3 col = colf;
    p0.xy *= rot(-PI / 4.);
    vec2 p = p0;
    p *= n;
    float numx = floor(p.x), numy = floor(p.y), num = floor(p.x) + floor(p.y);
    p = fract(p);
    float dd = 0.2, d2 = 10., d1 = 10.;
    
    col = draw_cyrcle(p0, 0.98, 0.02, col, col1);
    

    if(abs(p0.x) < (1. + dd) / n && abs(p0.y) < (1. + dd) / n) {
        d1 = min(p.x, 1. - p.x);
        d2 = min(p.y, 1. - p.y);
        d = min(d1, d2);
        if(d1 < dd && d2 < dd) {
            if(p.x < dd && p.y < dd) {
                d2 = p.y;
                d1 = p.x;
            }
            if(p.x < dd && 1. - p.y < dd) {
                d2 = 1. - p.y;
                d1 = p.x;
                num += 1.;
            }
            if(1. - p.x < dd && 1. - p.y < dd) {
                d2 = 1. - p.y;
                d1 = 1. - p.x;
                num += 2.;
            }
            if(1. - p.x < dd && p.y < dd) {
                d2 = p.y;
                d1 = 1. - p.x;
                num += 1.;
            }
            if(mod(num, 2.0) == 1.)
            {
                //d = d1;
                col = curColorN(d1 - dd, col, col1, dd);
                col = curColorN(d2 - dd, col, col1, dd);
            }
            else
            {
                //d = d2;
                col = curColorN(d2 - dd, col, col1, dd);
                col = curColorN(d1 - dd, col, col1, dd);
            }
        }
        else
        {
            col = curColorN(d - dd, col, col1, dd);
        }
    } else {
        if (p0.y > 0.)
            d = min(d, bridge(p0, n, dd));
        if (p0.y < 0.)    
            d = min(d, bridge(vec2(-p0.x, -p0.y), n, dd));
        if (p0.x > 0.)
            d = min(d, bridge(vec2(p0.y, p0.x), n, dd));    
        if (p0.x < 0.)
            d = min(d, bridge(vec2(-p0.y, -p0.x), n, dd));  

        
        if (p0.x < 0. && p0.y > 0.)          
        {
            d = min(bridge2(p0, n, dd), d);
        }
        if (p0.x > 0. && p0.y < 0.)          
        {
            d = min(bridge2(-p0, n, dd), d);
        }
        col = curColorN(d - dd, col, col1, dd);
        
    }
    
    
    
    return col;
}

vec3 vi5(vec2 p) {
    p.x *= TAU;
    vec3 col = colf;
    vec3 col0 = vec3(0.5, 0.5, 1.);
    vec3 col1 = vec3(1., 0.5, 0.5);
    vec3 col2 = vec3(0.5, 1., 0.5);
    float dd = 0.15, k = 0.25;
    float fa0 = 0., fa1 = -2. / 3. * PI, fa2 = -4. / 3. * PI;
    vec3 facol0 = col0, facol1 = col1, facol2 = col2;
    if(p.x >= PI / 3.) {
        fa0 = -4. / 3. * PI, fa1 = 0., fa2 = -2. / 3. * PI;
        facol0 = col2, facol1 = col0, facol2 = col1;
    }
    if(p.x >= 2. * PI / 3.) {
        fa1 = -4. / 3. * PI, fa2 = 0., fa0 = -2. / 3. * PI;
        facol1 = col2, facol2 = col0, facol0 = col1;
    }
    if(p.x >= PI) {
        fa2 = -4. / 3. * PI, fa0 = 0., fa1 = -2. / 3. * PI;
        facol2 = col2, facol0 = col0, facol1 = col1;
    }
    if(p.x >= 4. * PI / 3.) {
        fa0 = -4. / 3. * PI, fa1 = 0., fa2 = -2. / 3. * PI;
        facol0 = col2, facol1 = col0, facol2 = col1;
    }
    if(p.x >= 5. * PI / 3.) {
        fa1 = -4. / 3. * PI, fa2 = 0., fa0 = -2. / 3. * PI;
        facol1 = col2, facol2 = col0, facol0 = col1;
    }
    float y = k * cos(p.x + fa0), alf = PI / 2. - (atan(-sin(p.x + fa0), 1.));
    float d1 = abs(p.y - y) - dd / sin(alf);
    col = curColorN(d1, col, facol0, dd / sin(alf));
    y = k * cos(p.x + fa1), alf = PI / 2. - (atan(-sin(p.x + fa1), 1.));
    d1 = abs(p.y - y) - dd / sin(alf);
    col = curColorN(d1, col, facol1, dd / sin(alf));
    y = k * cos(p.x + fa2), alf = PI / 2. - (atan(-sin(p.x + fa2), 1.));
    d1 = abs(p.y - y) - dd / sin(alf);
    col = curColorN(d1, col, facol2, dd / sin(alf));
    return col;
}

vec3 vi6(vec3 p) {
    float h = 0.2, n = 15., d = TAU / n;
    p.xy *= rot(iTime * 0.1);

    float x = fract(mod(atan(p.y, p.x), TAU) / d);
    float y = (p.z) / h;
    vec3 col = vi5(vec2(x, y));
    return col;

}

float sdRoundedCylinder(vec3 p, float ra, float rb, float h) {
    vec2 d = vec2(length(p.xz) - 2.0 * ra + rb, abs(p.y) - h);
    return min(max(d.x, d.y), 0.0) + length(max(d, 0.0)) - rb;
}

float map(vec3 p) {
    return sdRoundedCylinder(p.xzy, 0.5, 0.03, 0.1);
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

void mainImage(out vec4 fragColor, in vec2 fragCoord) {
    vec3 light = normalize(vec3(0.0, 1.0, 1.)); //light
    vec3 light2 = normalize(vec3(0.0, -1.0, -1.)); //light
    //vec2 mo = 1.5 * cos(0.5 * iTime + vec2(0, 11));
    vec2 mo = vec2(-iTime * 0.2, 0.);
    
    //if  (iMouse.z > 0.0)
   
    {
        mo = (-iResolution.xy + 2.0 * (iMouse.xy)) / iResolution.y;
        mo*=1.7;
    }
   

    vec3 ro = vec3(0.0, 0.0, 1.85); // camera
    //camera rotation
    ro.yz *= rot(mo.y);
    ro.xz *= rot(-mo.x - 1.57);

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

                float d = length(pos.xy);
                if(d < 0.98)
                    col = endless(pos.xy);
                else if(d > 1.0)
                    col = vi6(pos);
                else {
                    vec3 col1 = endless(pos.xy);
                    vec3 col2 = vi6(pos);
                    float s = smoothstep(0.98, 1., d);
                    col = mix(col1, col2, vec3(s));
                }
                

                vec3 nor = calcNormal(pos);
                vec3 R = reflect(light, nor);
                float specular = pow(max(abs(dot(R, rd)), 0.), 55.);
                float difu = abs(dot(nor, light));
                col = col * (clamp(difu, 0., 1.0) + 0.5) + vec3(1., .7, .4) * specular;
                float fre = pow(clamp(dot(nor, rd) + 1., .0, 1.), 3.); // Fresnel, for some mild glow.
                col += vec3(.1, .1, 0.1) * fre; //?
                //col = sqrt(col);
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