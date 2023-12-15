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

//https://iquilezles.org/articles/smin/
// polynomial smooth min 1 (k=0.1)
float smin(float a, float b, float k) {
    float h = clamp(0.5 + 0.5 * (b - a) / k, 0.0, 1.0);
    return mix(b, a, h) - k * h * (1.0 - h);
}

float aafi(vec2 p) {
    float fi = atan(p.y, p.x);
    fi += step(p.y, 0.0) * TAU;
    return fi;
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
        float a = 0.2;
        float b = 0.15;
        return a*x+ b;
}


float f2(float x, float x1, float r, float h)
{
    return r*cos((x-x1)/h*PI/2.0);
}



float sdf_scale(vec3 p, float x0, float x1, float h)
{
    /*
    float x0 = 0.;
    float x1 = 0.5;
    float h = 0.2;
    */
    float k = 0.15;
    float x = clamp(p.x, x0, x1);
    float y = clamp(p.y, -f1(x), f1(x));
    float z = k*cos(y/f1(x)*PI/2.0)*(x-x0)/(x1-x0);
    
    float d = length(p - vec3(x, y, 0.));

    float d2 = length(p - vec3(x, y, z));
    d2 -= smoothstep(0.01, 0., abs(p.y)) * 0.01;   
    
    d = min(d, d2);
    
    float r = f1(x1);
   
    x = clamp(p.x, x1, x1+h);
    y = clamp(p.y, -f2(x, x1, r, h), f2(x, x1, r, h));
    d2 = length(p - vec3(x, y, 0.));
    d = min(d, d2);

    z = k*cos(y/f2(x, x1, r, h)*PI/2.0)*(1. - (x-x1)/h);
    d2 = length(p - vec3(x, y, z));
    d2 -= smoothstep(0.01, 0., abs(p.y)) * 0.01;   
    d = min(d, d2);
    d *= 0.3;
    d -= 0.005;
    return d;
}

float def3(in vec3 p, float x0, float x1, float h)
{
    float y = (p.y < x1)?p.y/f1(p.x): p.y/f2(p.x, x1, f1(x1), h);
    float x = (p.x-x0)/(x1 + h - x0);
    float d = sdf_scale(rotateX(-y*PI/4.*x)*p, x0, x1, h);
    return d;

}

float def2(in vec3 p, float x0, float x1, float h)
    {
        float x = (p.x-x0)/(x1 + h - x0);
        //x = pow(abs(x), 3.);
        float d = def3(rotateY(-x*PI/4.)*p, x0, x1, h);
        return d;
        
    }

float def1(in vec3 p, float x0, float x1, float h)
{
    /*
    float x0 = 0.;
    float x1 = 0.5;
    float h = 0.2;
    */
    
    float x = (p.x-x0)/(x1 + h - x0);
    const float k = PI/10.0; // or some other amount
    float c = cos(x*k);
    float s = sin(k*x);
    mat2  m = mat2(c,-s,s,c);
    vec3  q = vec3(m*p.zx,p.y);
    //q = p;
    float d = def2(q, x0, x1, h);
    return d;
}




float scale1(float x, float y)
{
    float f = step(0.25, x)*step(x, 0.75)*smoothstep(0., 0.1, abs(y-1.));
    if (f == 0.)
        f += smoothstep(-0.1, 0.1, 4.*x - y)*smoothstep(-0.1, 0.1, 4. - 4.*x - y);
    f = clamp(f, 0., 1.)*y*y*(1. + .5*sin(x*PI));
    float dl = f*0.01;
    return dl;   
}

float scale2(float x, float y)
{
    float r = 1. - length(vec2(x*2., y) - vec2(1., 0.));
    float f = smoothstep(-0.05, 0.05, r);
    f = f*y*y*y*(1. + .5*sin(x*PI));
    
    float dl = f*0.005;
    return dl;   
}

float scale3(float x, float y)
{
    float r = (1.-asin(abs((x - 0.5)*2.))*2./PI)-y;
    float f = smoothstep(-0.05, 0.05, r);
    f *= y*y;
    float dl = f*0.01;
    return dl;   
}

float sdHexagon( in vec2 p, in float r )
{
    const vec3 k = vec3(-0.866025404,0.5,0.577350269);
    p = abs(p);
    p -= 2.0*min(dot(k.xy,p),0.0)*k.xy;
    p -= vec2(clamp(p.x, -k.z*r, k.z*r), r);
    return length(p)*sign(p.y);
}


float scale4(float x, float y)
{
    float r = sdHexagon(vec2(x,y) - 0.5, 0.3) - 0.15; 
    float f = smoothstep(-0.05, 0.05, -r);
    if (f > 0.)
        f += y;
    //f *= (1. + 1.0*cos((y-0.5)*PI));
    float dl = f*0.01;
    return dl;   
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


float toorow(float x, float y, float st)
{
    float x1 = x;
    float y1 = y*0.5 + 0.5;
    float x2 = mod(x - st, 1.0);
    float y2 = y*0.5;
    float dl1 = scale2(x1, y1);
    float dl2 = scale2(x2, y2);
    return max(dl1, dl2);
}

float sdConePine(vec3 p, float a)
{
    float l = length(p.xy);
    float h = -p.z + a/2.;
    float d = sdCosNp(vec2(h, l), a);
    
    float n = 10.;
    float m = 8.;
    float st = 0.25;
    float y = 1.0 - h/a;
    float x = mod(atan(p.y, p.x), TAU)/TAU;
    float row = floor(y*n);
    y = y*n - row;
    float shift = mod(st*row, 1.0);
    x = mod(x - shift/m, 1.0);
    x = fract(x*m);
    //float dl = displ(x, y);
    float dl = toorow(x, y, st);
    d*=0.3;
    d -= dl;
    return d;
}

/*
float toorow3D(float x, float y, float z, float st, float w, float h)
{
    float x1 = x;
    float y1 = y*0.5 + 0.5;
    float x2 = mod(x - st, 1.0);
    float y2 = y*0.5;
    float dl1 = petal(vec3(x1*w, y1*h, z), w, h);
    float dl2 = petal(vec3(x2*w, y2*h, z), w, h);
    return min(dl1, dl2);
}
*/



float petal(vec3 p, float w, float h)
{
    float y = clamp(p.y/h, 0., 1.);
    float lh = .6;
    float yh = clamp((y - 1. + lh)/lh, 0., 1.); 
    float d = cos(yh*PI/2.0)/2.0;
    float x = clamp(p.x/w, 0.5-d, 0.5+d);
    float z = (1.+cos((x-0.5)*PI)*cos((y-0.5)*PI))*.2*(w+h)*y*y;
    vec3 val = vec3(x*w, y*h, z);
    return length(p - val)*0.4 - 0.01*(h+w);
}

float petal4(vec3 p, float w, float h)
{
    float y = clamp(p.y/h, 0., 1.);
    float lh = .6;
    float yh = y;//clamp((y - 1. + lh)/lh, 0., 1.); 
    //float d = 1.0; 
    float d = sqrt(1. - y*y); //cos(yh*PI/2.0)/2.0;
    float x = clamp(p.x/w, 0.5-d, 0.5+d);
    float z = (1.+cos((x-0.5)*PI)*cos((y-0.5)*PI))*.2*(w+h)*y*y;
    vec3 val = vec3(x*w, y*h, z);
    return length(p - val)*0.4 - 0.005*(h+w);
}


float petal3(vec3 p, float w, float h)
{
    float x = p.x/w;
    float y = p.y/h;

    float d = sdHexagon(vec2(x,y) - 0.5, 0.4)*(h+w)*0.5; 
    if (d>0.)
    {
        d =  length(vec2(p.z, d));
    }
    else
    {
        float z = (1.+cos((x-0.5)*PI)*cos((y-0.5)*PI))*.2*(w+h)*y*y;
        d =  (p.z - z);
    }
    return d*0.5 - 0.02*(h+w);
}


vec3 sdfColor;
vec3 resColor;
vec3 col1 = vec3(0.5019607843137255, 0.6705882352941176, 0.34509803921568627);
vec3 col2 = vec3(0.4588235294117647, 0.16862745098039217, 0.21176470588235294);
vec3 col3 = vec3(0.13333333333333333, 0.3254901960784314, 0.08235294117647059);
vec3 col4 = vec3(0.4, 0.6509803921568628, 0.10588235294117647);
float fi = PI/7.; 
float fcolor(float x)
{
    return 0.35 + 0.1*sin(x*10.0); 
}

float petalArti(vec3 pos, float w, float h)
{
    //PI/7.
    vec3 p = rotateX(fi)*pos;
    
    float y = clamp(p.y/h, 0., 1.);
    float lh = .5;
    float yh = clamp((y - 1. + lh)/lh, 0., 1.); 
    float d = cos(yh*(PI/2.2))/2.0;
    float x = clamp(p.x/w, 0.5-d, 0.5+d);
    float z = sin((d - abs(x - 0.5))/d*PI/2.1)*(w+h)*0.05 * cos(abs(y-0.6)/0.6*PI/2.);//(1.+cos((x-0.5)*PI)*cos((y-0.5)*PI))*.1*(w+h)*y*y;
    z*=(1. + y);
    
    float cl = fcolor(x) + fcolor(x + PI); 
    sdfColor = mix(col2, col1, vec3(smoothstep(-0.3, 0.3, y - cl)));
    
    vec3 val = vec3(x*w, y*h, z);

    
    return length(p - val)*0.4 - 0.001*(h+w);
}

float petalPine(vec3 pos, float w, float h)
{
    float rs = 4.*cos(iTime) + 8.;
    vec3 p = rotateX(PI/rs)*pos;
    
    float y = clamp(p.y/h, 0., 1.);
    float lh = .5;
    float yh = clamp((y - 1. + lh)/lh, 0., 1.); 
    float d = cos(yh*PI/2.)/2.0;
    
    float x = clamp(p.x/w, 0.5-d, 0.5+d);
    float midh = sin((0.5 - abs(x - 0.5))/0.5*PI/2.2)*(w+h)*0.03;

    float z = 0.;
    if (y < 1.-lh)
        z = midh*y/(1.-lh);
    else
    {
        float yyh = acos(2.*abs(x - 0.5))*2./PI*lh;
        z = midh*((yyh + (1.-lh))-y)/(yyh + 0.001);
    }
    z*=(1. + y)*(1. + y);
    
    sdfColor = mix(col3, col1, vec3(y));
    float pst = smoothstep(0.1, -0.1, (length(vec2(x-0.5, y - 0.5)) - 0.15));
    z += pst*.01*(h+w);
    sdfColor = mix(sdfColor, pow(col4, vec3(0.7)), pst);
    
    return length(p - vec3(x*w, y*h, z))*0.4 - 0.005*(h+w);
}

float sdConePine3D(vec3 p, float a)
{
    float l = length(p.xy);
    float h = -p.z + a/2.;
    float d = sdCosNp(vec2(h, l), a);
    
    float n = 13.;
    float m = 8.;
    float st = 0.5;
    float y = 1.0 - h/a;
    float x = mod(atan(p.y, p.x), TAU)/TAU;
    float row = floor(y*n);
    y = y*n - row;
    float shift = mod(st*row, 1.0);
    x = mod(x - shift/m, 1.0);
    x = fract(x*m);
    
    float rh = 2.*a/n;
    float rw = (l-d)*TAU/m;

    float x1 = x;
    float y1 = y*0.5 + 0.5;
    float x2 = mod(x - st, 1.0);
    float y2 = y*0.5;
    
    float dl1 = petalPine(vec3(x1*rw, y1*rh, d), rw, rh);
        resColor = sdfColor;
    float dl2 = petalPine(vec3(x2*rw, y2*rh, d), rw, rh);
    if (dl2 < dl1)
        resColor = sdfColor;

    float dl = min(dl1, dl2);
    
    d*=0.3;
    if (d < dl)
        resColor = col3;
    d = min(d, dl);
    return d;
}





float sdArtichoke(vec3 p, float r)
{
    float l = length(p.xy);
    float d = length(p)-r;
    
    float n = 8.;
    float m = 12.;
    float st = 0.5;
    
    float y = 1. - (atan(l, p.z) + PI/2.)/PI;
    float x = mod(atan(p.y, p.x), TAU)/TAU;
    float row = floor(y*n);
    y = y*n - row;
    
    float shift = mod(st*row, 1.0);
    x = mod(x - shift/m, 1.0);
    x = fract(x*m);
    float rh = 1.5*PI*r/n;///cos(fi);
    float rw = l/length(p)*r*TAU/m;
    //float dl = toorow3D(x, y, d, step, rw, rh);
    //calc patal in two row
    float x1 = x;
    float y1 = y*0.5 + 0.5;
    float x2 = mod(x - st, 1.0);
    float y2 = y*0.5;
    float dl1 = petalArti(vec3(x1*rw, y1*rh, d), rw, rh);
    resColor = sdfColor; 
    float dl2 = petalArti(vec3(x2*rw, y2*rh, d), rw, rh);
    if (dl2 < dl1)
        resColor = sdfColor;
    float dl = min(dl1, dl2);
    //calc patal
    d*=0.8;
    if (d < dl)
        resColor = col2;
    d = min(d, dl);
    return d;
    
}



float map(in vec3 pos) {
   
   //return def1(pos, 0., 1., 1.);
   //return sdConePine(pos, 2.0);
   //return petalArti(pos, 1., 1.);
   //return sdArtichoke(pos, .8);
   //return petal(pos, .6, .6);
   return sdConePine3D(pos, 2.);
   //return petalPine(pos, 1., 1.);
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
    vec3 f = normalize(l - p), r = normalize(vec3(f.z,0,-f.x)), u = cross(f, r), c = f * z, i = c + uv.x * r + uv.y * u;
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

vec3 calccolor(vec3 col_in, vec3 backcol, vec3 rd, vec3 light1, vec3 light2, vec3 nor)
{
    vec3 col = col_in;
    float d = dot(rd, nor);
    if (d < 0.0)
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
    float t = iTime / 3.;
    vec2 m = vec2(0.0, 0.0);
    //if  (iMouse.z > 0.0)
    {
        m = (-iResolution.xy + 2.0 * (iMouse.xy)) / iResolution.y;
        t = 0.;
    }
    vec3 ro = vec3(0.0, 0.0, 2.); // camera
    ro = rotateY(-m.x * TAU) * rotateX(-m.y * PI) * ro; //camera rotation

    const float fl = 1.5; // focal length
    float dist = dist_infin;
    mat3 rota  = rotateX(-PI/2.0)*rotateY(-t);
    mat3 rota_1  = rotateY(t)*rotateX(PI/2.0);
    

    vec3 tot = vec3(0.0);

    //antialiasing
    for(int m = 0; m < AA; m++) for(int n = 0; n < AA; n++) {
            vec2 o = vec2(float(m), float(n)) / float(AA) - 0.5;
            vec2 p = (-iResolution.xy + 2.0 * (fragCoord + o)) / iResolution.y;
            vec3 rd = GetRayDir(p, ro, vec3(0, 0., 0), fl); //ray direction
            vec3 col = vec3(1.); // background  

            HIT giper = giper3D(rota * ro, rota * rd);
            if(giper.dist < dist) {
                vec3 nor = rota_1 * giper.nor;
                col = resColor;
                col = calccolor(col, col, -rd, light, light2, nor);
            }
            tot += col;
        }
    
    tot = pow(tot, vec3(0.7))/float(AA);
    //antialiasing
    fragColor = vec4(tot, 1.0);
}
/////=====================================================================================
void main() {
    vec4 fragColor = vec4(0);
    mainImage(fragColor, gl_FragCoord.xy);
    gl_FragColor = fragColor;
}