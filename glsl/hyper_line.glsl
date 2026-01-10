#ifdef GL_ES
precision mediump float;
#endif



/////=====================================================================================
/*
*/

#define PI  3.14159265359
#define TAU 6.28318530718
#define rot(f) mat2(cos(f), -sin(f), sin(f), cos(f))

const float dist_infin = 40.0;
#define nn 128
const float eps = 0.01;

vec3 getSg(vec3 p, float nseg) {
    float fi = mod(atan(p.y, p.x), TAU);
    fi = mod(fi + PI / nseg, TAU);
    float n = floor(fi / TAU * nseg);
    p.xy *= rot(-n * TAU / nseg);
    return p;
}

float sdSegment(vec3 p, vec3 a, vec3 b) {
    
    vec3 pa = p - a, ba = b - a;
    float h = clamp(dot(pa, ba) / dot(ba, ba), 0.0, 1.0);
    return length(pa - ba * h);
}

float map(vec3 p) {
    p.yz*=rot(PI/2.);
    //float n1 = 14., R = 2., h = 2.1;
    float k = 3., r = 8./k, //big radius
          w = 0.07/k,//# radius edge
          n = 22.,//#count lines
          n2 = 4.,//  #section
          h = 16./k;//  #height hyperboloid

    float alf0 = n2/n*TAU, x0 = r*cos(alf0), y0 = r*sin(alf0), y1 = p.z/h * 2. * y0, alf1 = atan(y1, x0),
    dalf = alf0-alf1;
    vec3 p1 = vec3(p);
    p.xy *= rot(dalf);
    p1.xy *= rot(-dalf);
    p = getSg(p, n);
    p1 = getSg(p1, n);
    
    vec3 a0 = vec3(r, 0., h/2.), a1 = vec3(a0); 
    vec3 b0 = vec3(r*cos(n2/n * 2.* TAU), r*sin(n2/n * 2.* TAU), -h/2.);
    vec3 b1 = vec3(r*cos(n2/n * 2.* TAU), -r*sin(n2/n * 2.* TAU), -h/2.);
    b1.xy *= rot(dalf);
    a1.xy *= rot(dalf);
   
    a0.xy *= rot(-dalf);
    b0.xy *= rot(-dalf);
    float d = (sdSegment(p, a0, b0)* 0.8 - w);
    float d1 = (sdSegment(p1, a1, b1)* .8 - w);
    
    return min(d, d1);

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
#define AA 4
vec3 col0 = vec3(0.73, 0.59, 0.3);
vec3 col1 = vec3(0.75, 0.75, 0.75);
vec3 col2 = vec3(0.7, 0.7, 1.);
vec3 col4 = vec3(0.9, 0.1, 0.1);


void mainImage(out vec4 fragColor, in vec2 fragCoord) {
    vec3 light = normalize(vec3(0.0, .0, 1.)); //light
    vec3 light2 = normalize(vec3(0.0, 0.0, -1.)); //light
    vec2 mo = 1.5 * cos(0.5 * iTime + vec2(0, 11));
    if  (iMouse.z > 0.0)
    {
        mo = (-iResolution.xy + 2.0 * (iMouse.xy)) / iResolution.y;
        mo*=1.7;
    }
    vec3 ro = vec3(0.0, 0.0, 7.); // camera
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
                col = col0;
                vec3 nor = calcNormal(pos);
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
    //antialiasing
    fragColor = vec4(tot, 1.0);
}
/////=====================================================================================

void main() {
    vec4 fragColor = vec4(0);
    mainImage(fragColor, gl_FragCoord.xy);
    gl_FragColor = fragColor;
}
