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

/*
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
*/
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

    float k = 2.;
    float a0 = 1.*a*a*c + 3.*a*c*e + 1.*a*a*e + 1.*a*c*c + 1.*c*c*e + 1.*c*e*e + 1.*a*e*e-1.*a*c*t-1.*c*e*t-1.*a*e*t-1.*a*c*e*k;
    float a1 = 1.*a*a*d + 2.*a*b*c + 3.*a*c*f + 3.*a*d*e + 2.*a*b*e + 1.*a*a*f + 2.*a*c*d + 1.*b*c*c + 1.*c*c*f + 2.*c*d*e + 3.*b*c*e + 2.*c*e*f + 1.*d*e*e + 1.*b*e*e + 2.*a*e*f-1.*a*d*t-1.*b*c*t-1.*c*f*t-1.*d*e*t-1.*b*e*t-1.*a*f*t-1.*a*c*f*k-1.*a*d*e*k-1.*b*c*e*k;
    float a2 = 2.*a*b*d + 3.*a*d*f + 2.*a*b*f + 2.*b*c*d + 2.*c*d*f + 3.*b*c*f + 3.*b*d*e + 2.*d*e*f + 2.*b*e*f-1.*b*d*t-1.*d*f*t-1.*b*f*t + 1.*b*b*c + 1.*b*b*e + 1.*a*d*d + 1.*d*d*e + 1.*c*f*f + 1.*a*f*f-1.*a*d*f*k-1.*b*c*f*k-1.*b*d*e*k;
    float a3 = 1.*b*b*d + 3.*b*d*f + 1.*b*b*f + 1.*b*d*d + 1.*d*d*f + 1.*d*f*f + 1.*b*f*f-1.*b*d*f*k;

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
        // -3x^2 + 3y^2, 6yx, 3z^2 - t^2
        //normal = (xy+yz+zx) + (y+z)(x+y+z-a) - kyz, (xy+yz+zx) + (x+z)(x+y+z-a) - kxz,  (xy+yz+zx) + (y+x)(x+y+z-a) - kxy
        float s = pos.x*pos.y + pos.y*pos.z + pos.x*pos.z;
        float s1 = pos.x + pos.y + pos.z - t;
        nor = vec3(s1*(pos.y + pos.z) - k*pos.y*pos.z, s1*(pos.x + pos.z) - k*pos.x*pos.z, s1*(pos.x + pos.y) - k*pos.x*pos.y);
        nor += s;
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

    float ra = 3.0;
    float g = 1.0;

    float t = iTime/2.0;
    vec2 m = vec2(0.0, 0.0);
    //if  (iMouse.z > 0.0)
    {
        //m = (-iResolution.xy + 2.0*(iMouse.xy))/iResolution.y;
        //t = 0.;
    }
    vec3 ro = vec3(0.0, 0.0, 6.); // camera
    ro = rotateY(-m.x*TAU)*rotateX(-m.y*PI)*ro; //camera rotation
    
    
    const float fl = 1.5; // focal length
    float dist = dist_infin;
    float fi = PI/4.5;
    
    mat3 rota  = rotateX(fi)*rotateY(-fi)*rotateX(-PI/2.)*rotateY(t);
    mat3 rota_1  = rotateY(-t)*rotateX(PI/2.)*rotateY(fi)*rotateX(-fi);
    //mat3 sky = rotateZ(0.0)*rotateY(PI/2.0);
    
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