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

#define PI 3.14159265359
#define TAU 6.28318530718
#define nn 20.0
const float  eps = 0.05;
const float dist_infin = 200.0;
#define AA 1


mat3 rotateX(float f)
{
    return mat3(vec3(1.0,    0.0,      0.0), vec3(0.0,	 cos(f),  -sin(f)), 	vec3(.0, sin(f), cos(f)));
}

mat3 rotateY(float f)
{
    return mat3(vec3(cos(f), 0.0,  sin(f)),vec3(0.0,	 1.0,  0.0),vec3(-sin(f), 0.0, cos(f)));
}

mat3 rotateZ(float f)
{
    return mat3(vec3(cos(f),    -sin(f),  0.0),vec3(sin(f),	 cos(f),  0.0), 	vec3(0.0, 0.0, 1.0));
}

vec3 plane(vec3 ro, vec3 rd, vec3 po, vec3 nor)
{
    float t = dot(nor, (po - ro)) / dot(nor, rd);
    if (t < 0.0)
        return(vec3(dist_infin, dist_infin, dist_infin));
    return ro + rd*t;
}

float rand(float t, float n)
{
    float res =  fract(sin(t)*n);
    return res;
}
float getfi(float t, float n)
{
    if (t == 0.0)   
        return PI/3.0;
    float r = rand(t, n);
    return  r*PI/3.0 + PI/20.0;
}

float yBezier(vec2 v0, vec2 v1, vec2 v2, float x)
{
    //t^2(P0-2P1+P2) + t(2P1-2P0) + P0
    vec2 a2 = v0 - 2.0*v1 + v2;
    vec2 a1 = (2.0*v1 - 2.0*v0);
    vec2 a0 = v0;
    float t = x;
    float val = t*t*a2.y + t*a1.y + a0.y;
    return val;
}

float tileCurve(float p_x, float n)
{
    float t0 = floor(p_x);
    float t1 = t0 + 1.0;
    float px = fract(p_x);
    float a = getfi(t1, n) ;
    float b = getfi(t0, n) ;
    float x = tan(b)/(tan(a) + tan(b));
    float h = x*tan(a);
    x = 1.0 - x;
    float y = h;
    if (mod(t1, 2.0) == 0.0)
        y = - h;

    vec2 v0 = vec2(0.,0.);
    vec2 v1 = vec2(x,y);
    vec2 v2 = vec2(1.0,0.);

    float val = yBezier(v0, v1, v2, px);
    return val;
}


float sdfBezier(vec3 p)
{
    float u = tileCurve(p.x, 287.5668);
    float v = tileCurve(p.y, 117.3456);
    return abs(p.z - 15.0*u*v);
}


float sdf(vec3 p)
{
    float r = 4.0;
    if ((p.x*p.x + p.y*p.y) > r*r)
        return abs(p.z);

    float h =  r*r - p.x*p.x - p.y*p.y  ;   
    h = sqrt(h);
    return abs(p.z - h);
    //return abs(p.z - sin(p.x) - sin(p.y));
    //return sdfBezier(p);
}
float sdf2(vec3 p, float w, float h)
{
    return sdf(p);
    float x = (p.x + w)/2.0/w;
    float y = (p.y + h)/2.0/h;
    float z = texture(iChannel0, vec2(x,y)).r;
    return abs(p.z - 2.*z);
}

vec3 grid(vec3 ro, vec3 rd, float w, float h)
{
    vec3 col = vec3(0.0);
    //x axis
    vec3 nor = vec3(1.0, 0., 0.);
    for (float i = 0.; i < nn; i++)
    {
        
        vec3 po = vec3(-w + 2.0*w*i/nn, 0.0, 0.0);
        vec3 p = plane(ro, rd, po, nor);
        if (abs(p.y) > h)
            continue;
        
        float pst = smoothstep(eps, 0.0, sdf2(p, w, h));
        //float pst = step(sdf2(p, w, h), eps);
        col = mix(col, vec3(1.0), pst);
    }

    //y axis
    nor = vec3(0.0, 1., 0.);
    for (float i = 0.; i < nn; i++)
    {
        
        vec3 po = vec3(0.0, -h+ 2.0*h*i/nn, 0.0);
        vec3 p = plane(ro, rd, po, nor);
        if (abs(p.x) > w)
            continue;
        
        float pst = smoothstep(eps, 0.0, sdf2(p, w, h));
        //float pst = step(sdf2(p, w, h), eps);
        col = mix(col, vec3(1.0), pst);
    }
    return col;
}

vec3 GetRayDir(vec2 uv, vec3 p, vec3 l, float z) {
    vec3 
        f = normalize(l-p),
        r = normalize(cross(vec3(0,1,0), f)),
        u = cross(f,r),
        c = f*z,
        i = c + uv.x*r + uv.y*u;
    return normalize(i);
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    //vec2 p = (2.0*fragCoord-iResolution.xy)/iResolution.y;
    float t = iTime/4.0;
    vec2 m = vec2(0.0, 0.0);
    //if  (iMouse.z > 0.0)
    {
        m = (-iResolution.xy + 2.0*(iMouse.xy))/iResolution.y;
        t = 0.;
    }
    vec3 ro = vec3(0.0, -10.0, 5.); // camera
    const float fl = 1.5; // focal length
    ro = rotateY(-m.x*TAU)*rotateX(-m.y*PI)*ro; //camera rotation
    mat3 rota  = rotateZ(t + PI/4.0);
    //mat3 rota_1  = rotateY(t)*rotateZ(-t);
    float w = 8.;
    float h = 8.;
    vec3 tot = vec3(0.0);
    
    //antiblick
    for( int m=0; m<AA; m++ )
    for( int n=0; n<AA; n++ )
    {
            // pixel coordinates
        vec2 o = vec2(float(m),float(n)) / float(AA) - 0.5;
        vec2 p = (-iResolution.xy + 2.0*(fragCoord+o))/iResolution.y;
        vec3 rd = GetRayDir(p, ro, vec3(0,0.,0), fl); //ray direction
        vec3 col = grid(rota*ro, rota*rd, w, h);
        tot += col;
    }

    //antiblick
    tot /= float(AA*AA);
    fragColor = vec4(tot,1.0);
    
}

/////=====================================================================================
void main()
{
    vec4 fragColor = vec4(0);
    mainImage(fragColor,  gl_FragCoord.xy);
    gl_FragColor = fragColor;
}