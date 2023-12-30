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
float sdfReflect = 0.;
vec3 col1 = vec3(0.3764, 0.8196, 0.3725);
vec3 col2 = vec3(0.8117, 0.1764, 0.8078);
vec3 col4 = vec3(0.85, 0.75, 0.65);

//https://mathcurve.com/courbes2d/larme/larme.shtml
float larme2(vec2 p, float a)
{
    float k = 0.3;
    p.x += a*(k - 0.95)/2.0;
    float zoom = 1.2;
    float x = clamp(p.x, -a*0.95, a*k);
    float f = acos(x/a);
    float y = a*sin(f)*pow(sin(f/2.), 2.)*zoom;
    float d = length(p - vec2(x, y));
    if (abs(p.y) < abs(y))
    {
        sdfReflect = 0.;
        sdfColor = col2;
    }
    else
    {
        sdfReflect = 0.05;
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
            sdfReflect = 0.05;
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
    float l = length(p.xy);
    float d = larme2(vec2(p.z, l), a);
    return d*0.3 - 0.02;
}


float sdBacket(vec3 p, float a, float b, float m, float n) {
    float fi = atan(p.y, p.x); //aafi(p.xy)
    float w = 20.;
    for(float i = 0.; i < 10.; i++) {
        if(mod(i, m) == 0.0 && i > 0.)
            break;
        float t = (fi + TAU * i) / m;
        float wt = abs(p.z - b * sin(n * t));
        w = min(w, wt);
    }
    float r = length(vec2(length(p.xy) - a, w)) / 2.0;
    return r - 0.03;
}

float glzoom = 1.0;

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
float sdLonLat(vec3 p, float r) {
        #define ll 20.
    float fi = atan(p.x, p.y);
        //fi += step(p.y, 0.0)*TAU;
    fi = mod(fi, TAU);
    float ln = floor(fi / TAU * ll);
    float l1 = ln * TAU / ll;
    float l2 = l1 + TAU / ll;
    float d = min(sdRound(rotateX(l1) * rotateY(PI / 2.) * p, r, 0.), sdRound(rotateX(l2) * rotateY(PI / 2.) * p, r, 0.));

    fi = atan(p.z, length(p.xy));
    float mm = ll / 4.0;
    ln = floor(abs(fi) / PI * 2.0 * mm);
    l1 = ln * PI / 2.0 / mm;
    l2 = l1 + PI / 2.0 / mm;
    float d2 = min(sdRound(p, r, l1 * sign(p.z)), sdRound(p, r, l2 * sign(p.z)));
    d = min(d2, d);
    return d - 0.03;
}

float ff3(float x, float y)
{
    float k = 15.;
    return 0.1*(sin(k*x) + sin(k*y));
}

float sdEggBox(in vec3 p, float a, float b)
{
    float x = clamp(p.x, a, b);
    float y = clamp(p.y, a, b);
    vec3 val = vec3(x, y, ff3(x, y));
    float d = length(p - val)/4.;
    return d;
}

//https://www.shadertoy.com/view/tt23RR
float sdHexagram( in vec2 p, in float r )
{
    const vec4 k = vec4(-0.5,0.86602540378,0.57735026919,1.73205080757);
    
    p = abs(p);
    p -= 2.0*min(dot(k.xy,p),0.0)*k.xy;
    p -= 2.0*min(dot(k.yx,p),0.0)*k.yx;
    p -= vec2(clamp(p.x,r*k.z,r*k.w),r);
    return length(p)*sign(p.y);
}

//Extrussion
//https://www.shadertoy.com/view/4lyfzw
float sdHexagram3( in vec3 p, in float h, in float r )
{
    float d = sdHexagram(p.xy, r);
    vec2 w = vec2( d, abs(p.z) - h );
    return min(max(w.x,w.y),0.0) + length(max(w,0.0)) - 0.05;
}

float sdf(vec3 pos) {
    return sdHexagram3(pos, 0.15, 0.35);

}

float sdBox( in vec2 p, in vec2 b )
{
    vec2 d = abs(p)-b;
    return length(max(d,0.0)) + min(max(d.x,d.y),0.0);
    
}

//https://www.shadertoy.com/view/mlcyDj
vec4 texChar(float char, vec2 uv) {
    vec2 pt = uv/16.0;
    pt.x += char/16.0;
    pt.y += 12.0/16.0;
    return texture(iChannel0, pt);
}

float sdTextBox ( vec2 p, vec2 b, float char ) {
    float l = sdBox(p,b);
    vec2 pn = (p.xy / b.xy) * .5 + .5;
    float lt = (texChar(char, pn).w - 0.5);
    return max(lt,l); 
}
//===========================================

float sdBox3( in vec3 p, in vec2 b, float h, float char)
{
    float d = sdTextBox(p.xy, b, char);
    vec2 w = vec2( d, abs(p.z) - h );
    return min(max(w.x,w.y),0.0) + length(max(w,0.0)) - 0.02;
}
float sdCosNp(vec2 p, float a) {
    
    float fi = atan(p.y, p.x);
    float L = length(p.xy);
    float d = dist_infin;
    float r = a * cos(2. * fi);
    if(p.x < 0.)
        r = -.5;
    d = min(abs(L - r), d);
    float f = acos(L / a) / 2.;
    d = min(2.0 * abs(sin((fi - f) / 2.0)) * L, d);
    d = min(2.0 * abs(sin((fi + f) / 2.0)) * L, d);
    return d;
}

float scale5(float x, float y)
{
    float r = (1.-asin(abs((x - 0.5)*2.))*2./PI)*0.6 + 0.4 -y;
    float f = smoothstep(-0.05, 0.05, r);
    f *= y*y*0.5;
    f *= (1. + 6.0*cos((x-0.5)*PI)*cos((y-0.5)*PI));
    float dl = f*0.01;
    return dl;   
}


float toorow(float x, float y, float step)
{
    float x1 = x;
    float y1 = y*0.5 + 0.5;
    float x2 = mod(x + step, 1.0);
    float y2 = y*0.5;
    float dl1 = scale5(x1, y1);
    float dl2 = scale5(x2, y2);
    return max(dl1, dl2);
}

float sdConePine(vec3 p, float a)
{
    float l = length(p.xy);
    float h = -p.z + a/2.;
    float d = sdCosNp(vec2(h, l), a);
    
    float n = 10.;
    float m = 8.;
    float step = 0.5;
    float y = 1.0 - h/a;
    float x = mod(atan(p.y, p.x), TAU)/TAU;
    float row = floor(y*n);
    y = y*n - row;
    float shift = mod(step*row, 1.0);
    x = mod(x - shift/m, 1.0);
    x = fract(x*m);
    //float dl = displ(x, y);
    float dl = toorow(x, y, step);
    d*=0.3;
    d -= dl;
    return d;
}




float map(vec3 p) {
    float d = dist_infin;
    float glzoom2 = 1.0 - glzoom;
    p.x += glzoom-glzoom2;
    if (glzoom > 0.01)
        d = min(d, sdLonLat(p, 0.8*glzoom)/*sdEggBox(p, -glzoom, glzoom)larme(p, 0.8*glzoom) sdBox3(p, vec2(glzoom, glzoom), 0.2*glzoom, 9.)larme(vec3(p.x+glzoom-glzoom2, p.y, p.z), glzoom)*/);
    
    if (glzoom2 > 0.01)
        d = min(d, sdEggBox(p, -glzoom2, glzoom2)/*sdConePine(vec3(p.x+glzoom-glzoom2, p.y, p.z), 2.*glzoom2) sdEggBox(p, -glzoom2, glzoom2)*//*sdHexagram3(p, 0.15*glzoom2, 0.55*glzoom2)*/);     
    return d;
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


// IQ's vec2 to float hash.
float hash21(vec3 p){  
    return fract(sin(mod(dot(normalize(p), vec3(27.609, 57.583, 11.2345)), 6.2831853))*43758.5453); 
}

float npp = 200.;
float level = 0.95;
vec3 point(vec3 p) {
    return floor(p*npp)/npp;
}

vec3 hsb2rgb( in vec3 c )
{
    vec3 rgb = clamp(abs(mod(c.x*6.0+vec3(0.0,4.0,2.0),
                             6.0)-3.0)-1.0,
                     0.0,
                     1.0 );
    rgb = rgb*rgb*(3.0-2.0*rgb);
    return (c.z * mix( vec3(1.0), rgb, c.y));
}
float glz()
{
    float t = iTime/2.;
    float st = mod(floor(t), 4.);
    float res;
    if (st == 0.)
        res = 1.;
    if (st == 1.)    
        res = cos(fract(t)*PI/2.);
    if (st == 2.)    
        res = 0.;
    if (st == 3.)
        res = sin(fract(t)*PI/2.);    
    return res;    
}

void mainImage(out vec4 fragColor, in vec2 fragCoord) {
    //https://www.shadertoy.com/view/lcfGWH
    //vec3 snowBgcol = snowBackground( fragCoord ).rgb;
    //vec3 bg = snowBgcol;
    vec3 b1 = vec3(0.23529411764705882, 0.4235294117647059, 0.7725490196078432), b2 = vec3(0.3686274509803922, 0.5725490196078431, 0.8941176470588236);
    vec3 bg = mix(b2, b1*b1, fragCoord.y / iResolution.y);   

    vec3 light = normalize(vec3(0.0, .0, 1.)); //light
    vec3 light2 = normalize(vec3(0.0, 0.0, -1.)); //light
    vec2 mo = vec2(.0, 1.0);
    //if  (iMouse.z > 0.0)
    {
        mo = (-iResolution.xy + 2.0 * (iMouse.xy)) / iResolution.y;
    }
    glzoom = glz();//clamp(1.0 - clamp(cos(iTime/2.), 0., 1.0), 0.01, 1.0);
    vec3 ro = vec3(0.0, 0.0, 2.); // camera
    //camera rotation
    ro.yz *= rot(mo.y * PI);
    ro.xz *= rot(-mo.x * TAU);

    const float fl = 1.5; // focal length
    float dist = dist_infin;

    
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
                
                if(h < eps)
                {
                    vec3 pp = point(pos);
                    float fil = hash21(pp);
                    if (fil > level /*&& length(pos-pp) < 1./npp*/)
                        break;
                    else
                        h = 0.05;    
                }
                td += h;
                if (td >= dist_infin)
                    break;
            }
            if(td < dist_infin) {
                vec3 pp = point(pos);
                float fil = hash21(pp);
                if (fil > level /*&& length(pos-pp) < 1./npp*/)
                {
                    //float blink=1.0;//-cos(5.0*2.0*iTime);
                    //col = hsb2rgb(vec3(fract(fil*1000.)*3., 1., 2.)); //blink+0.9
                    col = col4;
                }
            }
            //==========================raymatch=============================
            tot += col;
        }
    //tot = sqrt(tot) / float(AA);
    tot = tot / float(AA);
    //antialiasing
    fragColor = vec4(tot, 1.0);
}
/////=====================================================================================
void main() {
    vec4 fragColor = vec4(0);
    mainImage(fragColor, gl_FragCoord.xy);
    gl_FragColor = fragColor;
}