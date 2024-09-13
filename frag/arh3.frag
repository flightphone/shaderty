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
SDF,raymatch,arch,architecture
SDF for simple architecture. 
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
float sdBox( vec3 p, vec3 b )
{
  b.z = b.z/2.;
  p.z -= b.z;
  vec3 q = abs(p) - b;
  return length(max(q,0.0)) + min(max(q.x,max(q.y,q.z)),0.0);
}

/*
float arch2(vec3 p, float R, float h)
{
    float res = 0.;
    res = min(length(vec2(p.x, max(p.z - h, 0.))), length(vec2(p.y, max(p.z - h, 0.)))) - R;
    return res;    
}
*/

float level2 (vec3 p, float R, float h, float w, float H)
{
    float t = sdBox(p, vec3(w, w, H));
    float fi = mod(atan(p.y, p.x), TAU);
    fi = mod(fi+PI/4., TAU);
    float n = floor(fi/(PI/2.));
    p.xy*= rot(-n*PI/2.);
    float t2 = length(vec2(p.y, max(p.z - h, 0.))) - R;
    //float t2 = arch2(p, R, h);
    return max (t, -t2);
}

float sdOctogon( in vec2 p, in float r )
{
    const vec3 k = vec3(-0.9238795325, 0.3826834323, 0.4142135623 );
    p = abs(p);
    p -= 2.0*min(dot(vec2( k.x,k.y),p),0.0)*vec2( k.x,k.y);
    p -= 2.0*min(dot(vec2(-k.x,k.y),p),0.0)*vec2(-k.x,k.y);
    p -= vec2(clamp(p.x, -k.z*r, k.z*r), r);
    return length(p)*sign(p.y);
}

float sdOctogon3( in vec3 p, in float r, float h)
{
    p.z -= h/2.;
    float d = sdOctogon(p.xy, r);
    vec2 w = vec2( d, abs(p.z) - h/2. );
    return min(max(w.x,w.y),0.0) + length(max(w,0.0));
}

float level8(vec3 p, float r, float h, float R, float H)
{
    float t = sdOctogon3(p, R, H);
    float fi = mod(atan(p.y, p.x), TAU);
    fi = mod(fi+PI/8., TAU);
    float n = floor(fi/(PI/4.));
    p.xy*= rot(-n*PI/4.);
    float t2 = length(vec2(p.y, max(p.z - h, 0.))) - r;
    return max (t, -t2);
    /*
    float t0 = arch2(p, r, h);
    p.xy *= rot(PI/4.);
    float t1 = arch2(p, r, h);
    float t2 = min(t1, t0);
    return max (t, -t2);
    */
    
    
}



//https://iquilezles.org/articles/distfunctions/
float dot2( in vec3 v ) { return dot(v,v); }
float udTriangle( vec3 p, vec3 a, vec3 b, vec3 c )
{
  vec3 ba = b - a; vec3 pa = p - a;
  vec3 cb = c - b; vec3 pb = p - b;
  vec3 ac = a - c; vec3 pc = p - c;
  vec3 nor = cross( ba, ac );

  
  return sqrt(
    (sign(dot(cross(ba,nor),pa)) +
     sign(dot(cross(cb,nor),pb)) +
     sign(dot(cross(ac,nor),pc))<2.0)
     ?
     min( min(
     dot2(ba*clamp(dot(ba,pa)/dot2(ba),0.0,1.0)-pa),
     dot2(cb*clamp(dot(cb,pb)/dot2(cb),0.0,1.0)-pb) ),
     dot2(ac*clamp(dot(ac,pc)/dot2(ac),0.0,1.0)-pc) )
     :
     dot(nor,pa)*dot(nor,pa)/dot2(nor) );
}
//https://iquilezles.org/articles/distfunctions/
float udQuad( vec3 p, vec3 a, vec3 b, vec3 c, vec3 d )
{
  vec3 ba = b - a; vec3 pa = p - a;
  vec3 cb = c - b; vec3 pb = p - b;
  vec3 dc = d - c; vec3 pc = p - c;
  vec3 ad = a - d; vec3 pd = p - d;
  vec3 nor = cross( ba, ad );

  return sqrt(
    (sign(dot(cross(ba,nor),pa)) +
     sign(dot(cross(cb,nor),pb)) +
     sign(dot(cross(dc,nor),pc)) +
     sign(dot(cross(ad,nor),pd))<3.0)
     ?
     min( min( min(
     dot2(ba*clamp(dot(ba,pa)/dot2(ba),0.0,1.0)-pa),
     dot2(cb*clamp(dot(cb,pb)/dot2(cb),0.0,1.0)-pb) ),
     dot2(dc*clamp(dot(dc,pc)/dot2(dc),0.0,1.0)-pc) ),
     dot2(ad*clamp(dot(ad,pd)/dot2(ad),0.0,1.0)-pd) )
     :
     dot(nor,pa)*dot(nor,pa)/dot2(nor) );
}

float dome(vec3 p, float R, float h)
{
    float fi = mod(atan(p.y, p.x), TAU);
    fi = mod(fi + PI/8., TAU);
    float n = floor(fi/(PI/4.)), fi0 = -PI/8. + n*PI/4., fi1 = PI/8. + n*PI/4.;
    vec3 a = vec3(R*cos(fi0), R*sin(fi0), 0.);
    vec3 b = vec3(R*cos(fi1), R*sin(fi1), 0.);
    vec3 c = vec3(0., 0., h);
    return udTriangle(p, a, b, c);
}

float roofh = 0.6/(2.39/2.), roofw = (2.39/2. - 0.44)/(2.39/2.);
float roof(vec3 p, float R, float r)
{
    
    float h = roofh, w = roofw;
    float fi = mod(atan(p.y, p.x), TAU), n = floor(fi/(PI/4.)), turn = floor((n + 1.)/2.);
    p.xy *= rot(-turn*PI/2.0);
    vec3 a = vec3(0., 0., R*h), b = vec3 (R, 0., R*h), c = vec3(0), d = vec3(0);
    vec3 a1 = vec3(0.), b1 = vec3(0.), c1 = vec3(0.);
    vec3 a2 = vec3(R, -R*w, 0.), b2 = vec3(R, R*w, 0.), c2 = vec3(R, 0., R*h);
    
    if (mod(n, 2.0) == 0.)
    {
        c = vec3(R, R*w, 0.);
        d = vec3(0, R*w, 0.);
        a1 = vec3(R*w, R*w, 0.0); 
        b1 = vec3(R, R*w, 0.); 
        c1 = vec3(R, R, 0.); 
        
    }
    else
    {
        c = vec3(R, -R*w, 0.);
        d = vec3(0, -R*w, 0.);
        a1 = vec3(R*w, -R*w, 0.0); 
        b1 = vec3(R, -R*w, 0.); 
        c1 = vec3(R, -R, 0.); 
        

        
    }
    float t0 = udQuad(p, a, b, c, d),
          t1 = udTriangle(p, a1, b1, c1),
          t2 = udTriangle(p, a2, b2, c2);
    return min((min(t0, t1) - r), t2);
}

float dome1R = 2.51/2. - 0.35, dome0R = 2.51/2. - 0.4, t0R = 2.51/2. - 0.4, t1R = 2.51/2.-0.05;
float t2R = (2.51-0.6*2.)/2., t2w = 2.51/2.;

float map(vec3 p) {
    p.yz *= rot(PI/2.);
    p.z -= 2.1;
    float d1 =  dome(p, dome1R, 0.9)-0.03;
    p.z += 1.38;
    float d0 =  level8(p, 0.17, 0.84, dome0R, 1.39);
    p.z += 0.7;
    float t0 = sdBox(p, vec3(t0R, t0R, 0.7))-0.03;
    float t1 = roof(p, t1R, 0.03) - 0.05;
    p.z+=2.6;
    float t2 =  level2(p, t2R, 1.7, t2w, 2.6);
    return min(min(min(min(t1, t2), t0), d0), d1);
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

//converts a vector on a sphere to longitude and latitude
vec2 lonlat (vec3 p)
{
    float lon = mod(atan(p.y, p.x), TAU)/TAU;
    float lat = atan(p.z, length(p.xy))/PI;
    return vec2(1.0-lon, lat);
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
    vec3 ro = vec3(0., 0., 5.5); // camera
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
            vec3 col = bg;
            
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
            //======================color====================================
            if(td < dist_infin) {
                //col = col1*col1;
                
                vec2 tx = lonlat(pos);
                tx = fract(tx * vec2(10., 10.));
                col = texture(iChannel0, tx).rgb;
                col = col*col;
                //if (abs(pos.z) < 0.1)
                //    col = vec3(1., 0., 0.);
                
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