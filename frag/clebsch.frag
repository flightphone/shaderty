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


struct HIT
{
    float dist;
    vec3 nor;
    vec3 pos;
};

float aafi(vec2 p)
{
    float l = length(p);
    float fi = asin(abs(p.y)/l);
    float pst = step(0.0, p.y)*step(p.x, 0.0);
    fi = fi + pst*(PI - 2.0*fi);
    pst = step(p.y, 0.0)*step(p.x, 0.0);
    fi = fi + pst*PI;
    pst = step(p.y, 0.0)*step(0.0, p.x);
    fi = fi + pst*(2.0*PI - 2.0*fi);
    return fi;    
}

//converts a vector on a sphere to longitude and latitude
vec2 lonlat (vec3 p)
{
    float lon = aafi(p.xy)/2.0/PI;
    float lat = aafi(vec2(p.z, length(p.xy)))/PI;
    return vec2(1.0-lon, lat);
}

const float dist_infin = 100000.0;
const HIT hit_inf = HIT(100000.0, vec3(0.0), vec3(0.0));

vec3 calcSkyReflect(vec3 rd, vec3 nor, mat3 sky)
{
    vec3 n = nor;
    float d = dot(rd, nor);
    n = nor*sign(d);
    vec3 r = reflect(rd, n);
    vec2 fon = lonlat(sky*r); //get longitude and latitude
    vec3 col = texture(iChannel0, fon).rgb;
    return col;

}

vec3 culccolor(vec3 col_in, vec3 backcol, vec3 rd, vec3 light1, vec3 light2, vec3 nor)
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

/// https://www.shadertoy.com/view/4lsfzj
int quadratic(float A, float B, float C, out vec3 x) {
   float D = B*B - 4.0*A*C;
   if (D < 0.0) return 0;
   D = sqrt(D);
   if (B < 0.0) D = -D;
   x[0] = (-B-D)/(2.0*A);
   x[1] = C/(A*x[0]);
   return 2;
}

// Numerical Recipes algorithm for solving cubic equation
int cubic0(float a, float b, float c, float d, out vec3 x) {
  if (a == 0.0) return quadratic(b,c,d,x);
  //if (d == 0.0) return quadratic(a,b,c,x); // Need 0 too.
  float tmp = a; a = b/tmp; b = c/tmp; c = d/tmp;
  // solve x^3 + ax^2 + bx + c = 0
  float Q = (a*a-3.0*b)/9.0;
  float R = (2.0*a*a*a - 9.0*a*b + 27.0*c)/54.0;
  float R2 = R*R, Q3 = Q*Q*Q;
  if (R2 < Q3) {
    float X = clamp(R/sqrt(Q3),-1.0,1.0);
    float theta = acos(X);
    float S = sqrt(Q); // Q must be positive since 0 <= R2 < Q3
    x[0] = -2.0*S*cos(theta/3.0) - a/3.0;
    x[1] = -2.0*S*cos((theta+2.0*PI)/3.0) - a/3.0;
    x[2] = -2.0*S*cos((theta+4.0*PI)/3.0) - a/3.0;
    return 3;
  } else {
    float A = -sign(R)*pow(abs(R)+sqrt(R2-Q3),0.3333);
    float B = A == 0.0 ? 0.0 : Q/A;
    x[0] = (A+B) - a/3.0;
    return 1;
  }
}

int cubic(float A, float B, float C, float D, out vec3 x) {
  int nroots;
  // Some ill-conditioned coeffs can cause problems
  // The worst is fixed by solving for reciprocal
  if (abs(A) > abs(D)) {
    nroots = cubic0(A,B,C,D,x);
  } else {
    nroots = cubic0(D,C,B,A,x);
    for (int i = 0; i < 3; i++) {
      x[i] = 1.0/x[i];
    }
  }
  return nroots;
}
/// https://www.shadertoy.com/view/4lsfzj


HIT giper3D(vec3 ro, vec3 rd, float t, float r)
{
    

    float a = ro.x;
    float b = rd.x;
    float c = ro.y;
    float d = rd.y;
    float e = ro.z;
    float f = rd.z;

    
    
    float k = 3.5;
    float a0 = 1.*a*a*a*k + 1.*c*c*c*k + 1.*e*e*e*k + 1.*k*t*t*t-1.*a*a*a-3.*a*a*c-3.*a*a*e-3.*a*a*t-3.*a*c*c-6.*a*c*e-6.*a*c*t-3.*a*e*e-6.*a*e*t-3.*a*t*t-1.*c*c*c-3.*c*c*e-3.*c*c*t-3.*c*e*e-6.*c*e*t-3.*c*t*t-1.*e*e*e-3.*e*e*t-3.*e*t*t-1.*t*t*t;
    float a1 = 3.*a*a*b*k + 3.*c*c*d*k + 3.*e*e*f*k-3.*a*a*b-3.*a*a*d-3.*a*a*f-6.*a*b*c-6.*a*c*d-6.*a*c*f-6.*a*b*e-6.*a*d*e-6.*a*e*f-6.*a*b*t-6.*a*d*t-6.*a*f*t-3.*b*c*c-3.*c*c*d-3.*c*c*f-6.*b*c*e-6.*c*d*e-6.*c*e*f-6.*b*c*t-6.*c*d*t-6.*c*f*t-3.*b*e*e-3.*d*e*e-3.*e*e*f-6.*b*e*t-6.*d*e*t-6.*e*f*t-3.*b*t*t-3.*d*t*t-3.*f*t*t;
    float a2 = 3.*a*b*b*k + 3.*c*d*d*k + 3.*e*f*f*k-3.*a*b*b-6.*a*b*d-6.*a*b*f-3.*a*d*d-6.*a*d*f-3.*a*f*f-3.*b*b*c-6.*b*c*d-6.*b*c*f-3.*c*d*d-6.*c*d*f-3.*c*f*f-3.*b*b*e-6.*b*d*e-6.*b*e*f-3.*d*d*e-6.*d*e*f-3.*e*f*f-3.*b*b*t-6.*b*d*t-6.*b*f*t-3.*d*d*t-6.*d*f*t-3.*f*f*t;      
    float a3 = 1.*b*b*b*k + 1.*d*d*d*k + 1.*f*f*f*k-1.*b*b*b-3.*b*b*d-3.*b*b*f-3.*b*d*d-6.*b*d*f-3.*b*f*f-1.*d*d*d-3.*d*d*f-3.*d*f*f-1.*f*f*f;
    
    vec3 roots = vec3(dist_infin);
    int nroots = cubic(a3, a2, a1, a0, roots);
    
    
    float dist = dist_infin;
    vec3 pos = vec3(0.0);
    vec3 nor = vec3(0.0);
    for (int i = 0; i < 3; i++)
    {
        if (i >= nroots)
            break;
        if (roots[i] < 0.0)
            continue;
        vec3 p = ro + roots[i]*rd;
        if (length(p) > r)    
            continue;
        if (roots[i] < dist)    
        {
            dist = roots[i];
            pos = p;
        }

    }
    if (dist < dist_infin)
    {
        float s = pos.x + pos.y + pos.z + t;
        s = s*s*3.0;
        nor = k*3.0*pos*pos - s;
        nor = normalize(nor);
    }
    return HIT(dist, nor, pos);
    
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
    vec3 light = normalize(vec3(0.0, 0.0, -1.0)); //light
    vec3 light2 = normalize(vec3(0.0, 0.0, 1.0)); //light

    float ra = 12.0;
    float g = 2.0;

    float t = iTime/2.0;
    vec2 m = vec2(0.0, 0.0);
    //if  (iMouse.z > 0.0)
    {
    m = (-iResolution.xy + 2.0*(iMouse.xy))/iResolution.y;
    t = 0.;
    }
    vec3 ro = vec3(0.0, 0.0, 30.); // camera
    ro = rotateY(-m.x*TAU)*rotateX(-m.y*PI)*ro; //camera rotation
    
    
    const float fl = 1.5; // focal length
    float dist = dist_infin;
    mat3 rota  = rotateZ(t)*rotateY(-t);
    mat3 rota_1  = rotateY(t)*rotateZ(-t);
    mat3 sky = rotateZ(0.0)*rotateX(PI/2.0);
    
    vec2 torus = vec2(1.0,0.3);
    vec3 tot = vec3(0.0);
    
    #define AA 2
    //antiblick
    for( int m=0; m<AA; m++ )
    for( int n=0; n<AA; n++ )
    {
        vec3 col = vec3(0.7, 0.7, 0.9); // background        
            // pixel coordinates
        vec2 o = vec2(float(m),float(n)) / float(AA) - 0.5;
        vec2 p = (-iResolution.xy + 2.0*(fragCoord+o))/iResolution.y;
        //vec2 p = (2.0*fragCoord-iResolution.xy)/iResolution.y;

        //vec3 rd = normalize( vec3(p,fl) ); // ray direction
        vec3 rd = GetRayDir(p, ro, vec3(0,0.,0), fl); //ray direction
        HIT giper = giper3D(rota*ro, rota*rd, g, ra);
        if (giper.dist < dist)
        {
            col = vec3(0.5, 0.5, 1.0);
            vec3 backcol = vec3(1.0, 0.2, 0.2);
            vec3 nor = rota_1*giper.nor;
            col = culccolor(col, backcol, -rd, light, light2, nor);
            //col = culccolor(col, backcol, -rd, light2, nor);
            // gamma
            //col = pow( col, vec3(0.4545) ); 
            //reflect
            //col = calcSkyReflect(-rd, nor, sky);
        }
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