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

const float dist_infin = 3.0;
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

float sdCappedCylinder( vec3 p, float h, float r )
{
  vec2 d = abs(vec2(length(p.xy),p.z)) - vec2(r,h);
  return min(max(d.x,d.y),0.0) + length(max(d,0.0));
}

float sdTorus( vec3 p, vec2 t )
{
  vec2 q = vec2(length(p.xy)-t.x,p.z);
  return length(q)-t.y;
}

//https://iquilezles.org/articles/smin/
// polynomial smooth min 1 (k=0.1)
float smin( float a, float b, float k )
{
    float h = clamp( 0.5+0.5*(b-a)/k, 0.0, 1.0 );
    return mix( b, a, h ) - k*h*(1.0-h);
}

float sdOctahedron( vec3 p, float s)
{
  p = abs(p);
  return (p.x+p.y+p.z-s)*0.57735027;
}

float map( in vec3 pos )
{

    float d = sdOctahedron(pos, 0.7);
    return d;
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

vec2 lonlat (vec3 p)
{
    float lon = aafi(p.xy)/TAU;
    float lat = aafi(vec2(p.z, length(p.xy)))/PI;
    return vec2(1.0-lon, lat);
}

HIT plane(vec3 ro, vec3 rd, vec3 po, vec3 nor)
{
    float t = dot(nor, (po - ro)) / dot(nor, rd);
    if (t < 0.0)
        return hit_inf;
    vec3 pos = ro + t*rd;
    return HIT(t, nor, pos);    
    
}
vec3 skycol(vec3 rd)
{
    vec2 fon = lonlat(rd); //get longitude and latitude
    return  texture(iChannel0, fon).rgb;
    //return  texture(iChannel0, rd).rgb;
}
float sgn(float x) {
  return x < 0.0? -1.0: 1.0; // Return 1 for x == 0
}

int quadratic(float A, float B, float C, out vec2 res) {
  float x1,x2;
  float b = -0.5*B;
  float q = b*b - A*C;
  if (q < 0.0) return 0;
  float r = b + sgn(b)*sqrt(q);
  if (r == 0.0) {
    x1 = C/A; x2 = -x1;
  } else {
    x1 = C/r; x2 = r/A;
  }
  res = vec2(x1,x2);
  return 2;
}

HIT cyl3D(vec3 ro, vec3 rd, float t)
{
    float a = ro.x;
    float b = rd.x;
    float c = ro.y;
    float d = rd.y;
    float e = ro.z;
    float f = rd.z;
    
    //https://github.com/flightphone/shaderty/blob/master/staples_polynomial.py
    //for generate this expression used python script staples_polynomial.py

    float a0 = 1.*a*a + 1.*c*c + 1.*e*e-1.*t*t;
    float a1 = 2.*a*b + 2.*c*d + 2.*e*f;
    float a2 = 1.*b*b + 1.*d*d + 1.*f*f;

    //https://github.com/flightphone/shaderty/blob/master/staples_polynomial.py

    vec2 roots = vec2(dist_infin);
    int nroots = quadratic(a2, a1, a0, roots); //cubic(a3, a2, a1, a0, roots);
        
    float dist = dist_infin;
    vec3 pos = vec3(0.0);
    vec3 nor = vec3(0.0);
    for (int i = 0; i < 2; i++)
    {
        if (i >= nroots)
            break;
        if (roots[i] < 0.0)
            continue;
        vec3 p = ro + roots[i]*rd;
        
        if (roots[i] < dist)    
        {
            dist = roots[i];
            pos = p;
        }
    }
    
    if (dist < dist_infin)
    {
        nor = vec3(0.+2.*pos.x, 0.+2.*pos.y, 0.+2.*pos.z);
        nor = normalize(nor);
    }
    
    return HIT(dist, nor, pos);
    
}

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
    //return HIT(t + giper.dist, nor, pos);
    return HIT(t, nor, pos);
}


void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    vec3 light = normalize(vec3(0.0, 1.0, 1.0)); //light
    vec3 light2 = normalize(vec3(0.0, 1.0, -1.0)); //light
    float t = -iTime*0.3;
    vec2 m = vec2(0.0, 0.0);
    //if  (iMouse.z > 0.0)
    {
        m = (-iResolution.xy + 2.0*(iMouse.xy))/iResolution.y;
        //t = 0.;
    }
    vec3 ro = vec3(0.0, 0.0, 2.1); // camera
    //ro = rotateY(-m.x*TAU)*rotateX(-m.y*PI)*ro; //camera rotation
    
    
    const float fl = 1.5; // focal length
    
    mat3 rota  = rotateX(PI/2.0)*rotateZ(m.x*TAU)*rotateX(-m.y*PI);
    mat3 rota_1  = rotateX(m.y*PI)*rotateZ(-m.x*TAU)*rotateX(-PI/2.0);
    //mat3 sky = rotateY(t);
    mat3 sky = rotateZ(0.)*rotateX(PI/2.0);
    
    vec3 tot = vec3(0.0);
    
    #define AA 1
    
        //antiblick
        for( int m=0; m<AA; m++ )
        for( int n=0; n<AA; n++ )
        {
            float dist = dist_infin;
            vec2 o = vec2(float(m),float(n)) / float(AA) - 0.5;
            vec2 p = (-iResolution.xy + 2.0*(fragCoord+o))/iResolution.y;
            vec3 rd = GetRayDir(p, ro, vec3(0,0.,0), fl); //ray direction
            vec3 col = skycol(sky*rd); 
            
            
            HIT giper = cyl3D(rota*ro, rota*rd, 0.7);
            
            if (giper.dist >= dist_infin)
            {
                tot += col;
                continue;
            }
            
                        
            giper = giper3D(rota*ro, rota*rd);
            if (giper.dist < dist)
            {
                
                vec3 nor = rota_1*giper.nor;
                float dif = clamp( dot(nor,light), 0.2, 1.0 );
                //float amb = 0.5 + 0.5*dot(nor,light2);
                
                //vec3 tx = vec3(80./255.,200.0/255.,120.0/255.);
                vec3 tx = vec3(0.698,0.098, 0.176);
                tx = tx*dif;
                
                //refract
                float n12 = 1.5;
                vec3 rd1 = refract(rd, nor, n12);
                vec3 ro1 = ro + giper.dist * rd + 0.7 * dist_infin * rd1;
                rd1 = -rd1;
                giper = giper3D(rota*ro1, rota*rd1);
                if (giper.dist < dist_infin)
                {
                    nor = rota_1*giper.nor;
                    ro1 = ro1 + rd1 * giper.dist;
                    rd1 = refract(-rd1, -nor, 1.0/n12);
                    
                    vec3 tx2 = skycol(sky*rd1); 
                    col = mix(tx, tx2, 0.6);
                }
                else
                    col = tx;

            }
            
            // gamma        
            //col = sqrt( col );
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