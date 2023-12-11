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

const float dist_infin = 6.0;
#define nn 128
const float eps = 0.001;


float sdPolar(vec3 p, float a) {
    float fi = atan(p.y, p.x);//aafi(p.xy);
    float r = a * cos(fi) * cos(fi);
    
    float L = length(p.xy);
    float d = abs(L - r);
    float fi2 = acos(sqrt(L / a));
    float fi3 = PI - fi2;
    d = min(2.0 * abs(sin((fi - fi2) / 2.0)) * L, d);
    d = min(2.0 * abs(sin((fi - fi3) / 2.0)) * L, d);
    d = min(2.0 * abs(sin((fi + fi2) / 2.0)) * L, d);
    d = min(2.0 * abs(sin((fi + fi3) / 2.0)) * L, d);
    
    d = sqrt(p.z * p.z  + d * d);
    d *= .5;
    d -= 0.03;
    return d;
}



float sdBacket(vec3 p, float a, float b, float m, float n)
{
    float fi = atan(p.y, p.x); //aafi(p.xy)
    float w = 20.;
    for (float i = 0.; i < 10.; i++ )
    {
        if (mod(i,m) == 0.0 && i > 0.)
            break;
        float t = (fi + TAU*i)/m;    
        float wt = abs(p.z - b*sin(n*t));
        w = min(w, wt);
    }
    float r = length(vec2(length(p.xy) - a, w))/2.0;
    return r - 0.03;
}

float sdCosN(vec3 p, float a, float n) {
    
    float fi = atan(p.y, p.x);
    float L = length(p.xy);
    float d = dist_infin;
    fi += step(p.y, 0.0)*TAU;
   
    for (float i = 0.; i < 10.; i++)
    {   
        if (fract(n*i) == 0. && i > 0.)
            break;  
        float r = a * cos(n*fi + n*i*TAU);
        d = min(abs(L - r), d);
    }
    float f = acos(L/a)/n;
    for (float i = 0.; i < 10.; i++)
    {
        d = min(2.0 * abs(sin((fi - (f + i*TAU/n)) / 2.0)) * L, d);    
        d = min(2.0 * abs(sin((fi + (f + i*TAU/n)) / 2.0)) * L, d);    
        d = min(2.0 * abs(sin((fi - (f - i*TAU/n)) / 2.0)) * L, d);    
        d = min(2.0 * abs(sin((fi + (f - i*TAU/n)) / 2.0)) * L, d);    
    }

    d = sqrt(p.z * p.z  + d * d);
    d *= .7;
    d -= 0.02;
    return d;
}


float sdCosN2(vec2 p, float a, float n) {
    
    float fi = atan(p.y, p.x);
    float L = length(p.xy);
    float d = dist_infin;
    fi += step(p.y, 0.0)*TAU;
   
    for (float i = 0.; i < 5.; i++)
    {   
        if (fract(n*i) == 0. && i > 0.)
            break;  
        float r = a * cos(n*fi + n*i*TAU);
        d = min(abs(L - r), d);
    }
    float f = acos(L/a)/n;
    for (float i = 0.; i < 2.; i++)
    {
        d = min(2.0 * abs(sin((fi - (f + i*TAU/n)) / 2.0)) * L, d);    
        d = min(2.0 * abs(sin((fi + (f + i*TAU/n)) / 2.0)) * L, d);    
        d = min(2.0 * abs(sin((fi - (f - i*TAU/n)) / 2.0)) * L, d);    
        d = min(2.0 * abs(sin((fi + (f - i*TAU/n)) / 2.0)) * L, d);    
    }
    return d;
}

float opExtrusion( in vec3 p,  float a, float n, float h)
{
    float d = sdCosN2(p.xy, a, n);
    vec2 w = vec2( d, abs(p.z) - h );
    d =  min(max(w.x,w.y),0.0) + length(max(w,0.0));
    d *= .5;
    d -= 0.01;
    return d;
}

float map( in vec3 pos )
{
    //return sdPolar(pos, 1.0);
    //return sdCosN(pos, 1.0, 0.25);
    return opExtrusion(pos, 1.0, 0.25, 0.15);

    //return sdBacket(pos, .9, .4, 3.5, 9.);
    //return sdEggBox(pos, -1., 1.);
}

// https://iquilezles.org/articles/normalsSDF
vec3 calcNormal( in vec3 pos )
{
    const float h = 0.0001; // replace by an appropriate value
    const vec2 k = vec2(1,-1);
    return normalize( k.xyy*map(pos + k.xyy*h ) + 
                      k.yyx*map(pos + k.yyx*h ) + 
                      k.yxy*map(pos + k.yxy*h ) + 
                      k.xxx*map(pos + k.xxx*h ) );
}

struct HIT
{
    float dist;
    vec3 nor;
    vec3 pos;
};
const HIT hit_inf = HIT(dist_infin, vec3(0.0), vec3(0.0));

HIT giper3D(vec3 ro, vec3 rd)
{
    float t  = 0.;
    for (int i = 0; i < nn; i++)
    {
        vec3 pos = ro + rd*t;
        float h = map(pos);
        if (h < eps || t >= dist_infin)
            break;
        t += h;  
    }    

    if (t >= dist_infin)
        return hit_inf;
      
    vec3 pos = ro + t*rd;
    vec3 nor = calcNormal(pos);
    return HIT(t, nor, pos);
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


mat3 rotateX(float f)
{
    return mat3(
    vec3(1.0,    0.0,      0.0),
    vec3(0.0,	 cos(f),  -sin(f)), 	
	vec3(.0, sin(f), cos(f))
    );
}


mat3 rotateZ(float f)
{
    return mat3(
    vec3(cos(f),    -sin(f),  0.0),
    vec3(sin(f),	 cos(f),  0.0), 	
	vec3(0.0, 0.0, 1.0)
    );
    
}


mat3 rotateY(float f)
{
    return mat3(
    vec3(cos(f), 0.0,  sin(f)),
    vec3(0.0,	 1.0,  0.0), 	
	vec3(-sin(f), 0.0, cos(f))
    );
}

/*
#if HW_PERFORMANCE==0
#define AA 1
#else
#define AA 2
#endif
*/
#define AA 1
void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    vec3 light = normalize(vec3(0.0, 1.0, 1.0)); //light
    vec3 light2 = normalize(vec3(0.0, 1.0, -1.0)); //light
    float t = iTime;
    vec2 m = vec2(0.0, 0.0);
    //if  (iMouse.z > 0.0)
    {
        m = (-iResolution.xy + 2.0*(iMouse.xy))/iResolution.y;
        t = 0.;
    }
    vec3 ro = vec3(0.0, 0.0, 2.0); // camera
    ro = rotateY(-m.x*TAU)*rotateX(-m.y*PI)*ro; //camera rotation
    
    
    const float fl = 1.5; // focal length
    float dist = dist_infin;
    mat3 rota  = rotateX(PI/2.0)*rotateZ(t)*rotateX(-t);
    mat3 rota_1  = rotateX(t)*rotateZ(-t)*rotateX(-PI/2.0);
    
    vec3 tot = vec3(0.0);
    
    
    //antiblick
    for( int m=0; m<AA; m++ )
    for( int n=0; n<AA; n++ )
    {
        vec2 o = vec2(float(m),float(n)) / float(AA) - 0.5;
        vec2 p = (-iResolution.xy + 2.0*(fragCoord+o))/iResolution.y;
        vec3 rd = GetRayDir(p, ro, vec3(0,0.,0), fl); //ray direction
        vec3 col = vec3(0.0);
        
        HIT giper = giper3D(rota*ro, rota*rd);
        if (giper.dist < dist)
        {
            vec3 nor = rota_1*giper.nor;
            float dif = clamp( dot(nor,light), 0.2, 1.0 );
            float amb = 0.5 + 0.5*dot(nor,light2);
            col = vec3(0.2,0.3,0.4)*amb + vec3(0.85,0.75,0.65)*dif;
        }
        // gamma        
        col = sqrt( col );
	    tot += col;
    }
    //antiblick
    tot /= float(AA*AA);
    fragColor = vec4(tot,1.0);
}
/////=====================================================================================
void main() {
    vec4 fragColor = vec4(0);
    mainImage(fragColor, gl_FragCoord.xy);
    gl_FragColor = fragColor;
}