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
#define TAU 6.283185

const float dist_infin = 4.0;
#define nn 128
const float eps = 0.001;

float dot2( in vec2 v ) { return dot(v,v); }
float dot2( in vec3 v ) { return dot(v,v); }
float ndot( in vec2 a, in vec2 b ) { return a.x*b.x - a.y*b.y; }

//https://iquilezles.org/articles/distfunctions/
float sdSphere( vec3 p, float s )
{
  return length(p)-s;
}

float sdEllipsoid( vec3 p, vec3 r )
{
  float k0 = length(p/r);
  float k1 = length(p/(r*r));
  return k0*(k0-1.0)/k1;
}

float sdCappedCone( vec3 p, float h, float r1, float r2 )
{
  vec2 q = vec2( length(p.xy), p.z );
  vec2 k1 = vec2(r2,h);
  vec2 k2 = vec2(r2-r1,2.0*h);
  vec2 ca = vec2(q.x-min(q.x,(q.y<0.0)?r1:r2), abs(q.y)-h);
  vec2 cb = q - k1 + k2*clamp( dot(k1-q,k2)/dot2(k2), 0.0, 1.0 );
  float s = (cb.x<0.0 && ca.y<0.0) ? -1.0 : 1.0;
  return s*sqrt( min(dot2(ca),dot2(cb)) );
}

float sdRoundedCylinder( vec3 p, float ra, float rb, float h )
{
  vec2 d = vec2( length(p.xy)-2.0*ra+rb, abs(p.z) - h );
  return min(max(d.x,d.y),0.0) + length(max(d,0.0)) - rb;
}

//https://iquilezles.org/articles/smin/
// polynomial smooth min 1 (k=0.1)
float smin( float a, float b, float k )
{
    float h = clamp( 0.5+0.5*(b-a)/k, 0.0, 1.0 );
    return mix( b, a, h ) - k*h*(1.0-h);
}

float map( in vec3 pos )
{
    vec3 pos1 = pos - vec3(0., 0., 0.3);
    float d1 = sdSphere(pos1 - vec3(0., 0., -.85), 0.2);
    float d2 = sdEllipsoid(pos1 - vec3(0., 0.0, -0.65), vec3(0.3, 0.3, 0.05));
    d1 = smin(d1, d2, 0.02);
    d2 = sdCappedCone(pos1-vec3(0.0, 0.0, -0.4), 0.2, 0.15, 0.22);
    //d1 = smin(d1, d2, 0.01);
    d1 = smin(d1, d2, 0.02);
    d2 = sdRoundedCylinder(pos1 - vec3(0.0, 0.0, 0.0), 0.17, 0.05, 0.15);
    d1 = min(d1, d2);
    d2 = sdEllipsoid(pos1 - vec3(0., 0.0, 0.), vec3(0.5, 0.5, 0.2));
    d1 = smin(d1, d2, 0.05);
    d2 = sdRoundedCylinder(pos1 - vec3(0.0, 0.0, 0.17), 0.25, 0.01, 0.02);
    d1 = smin(d1, d2, 0.02);
    return d1;
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

float aafi(vec2 p) {
    float fi = atan(p.y, p.x);
    fi += step(p.y, 0.0)*TAU;
    return fi;
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    vec3 light = normalize(vec3(0.0, 1.0, 1.0)); //light
    vec3 light2 = normalize(vec3(0.0, 1.0, -1.0)); //light
    float t = iTime;
    vec2 m = vec2(0.0, 0.0);
    //if  (iMouse.z > 0.0)
    {
        m = (-iResolution.xy + 2.0*(iMouse.xy))/iResolution.y;
        //t = 0.;
    }
    vec3 ro = vec3(0.0, 0.0, 1.8); // camera
    ro = rotateY(-m.x*TAU)*rotateX(-m.y*PI)*ro; //camera rotation
    
    
    const float fl = 1.5; // focal length
    float dist = dist_infin;
    mat3 rota  = rotateX(PI/2.0)*rotateZ(t)*rotateX(-t);
    mat3 rota_1  = rotateX(t)*rotateZ(-t)*rotateX(-PI/2.0);
    
    vec3 tot = vec3(0.0);
    
    #define AA 2
    //antiblick
    for( int m=0; m<AA; m++ )
    for( int n=0; n<AA; n++ )
    {
        vec2 o = vec2(float(m),float(n)) / float(AA) - 0.5;
        vec2 p = (-iResolution.xy + 2.0*(fragCoord+o))/iResolution.y;
        vec3 rd = GetRayDir(p, ro, vec3(0,0.,0), fl); //ray direction
        vec3 col = vec3(0.5, 0.5, 0.7);
        
        HIT giper = giper3D(rota*ro, rota*rd);
        if (giper.dist < dist)
        {
            vec3 nor = rota_1*giper.nor;
            float dif = clamp( dot(nor,light), 0.2, 1.0 );
            float amb = 0.5 + 0.5*dot(nor,light2);
            float yy = aafi(vec2(giper.pos.x, giper.pos.y)) /TAU;
            float xx = fract(giper.pos.z);
            vec3 tx = texture(iChannel0, vec2(xx,yy)).rgb;
            
            col = tx*amb + tx*dif;
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
void main()
{
    vec4 fragColor = vec4(0);
    mainImage(fragColor,  gl_FragCoord.xy);
    gl_FragColor = fragColor;
}