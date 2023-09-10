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
//roman  surface https://mathcurve.com/surfaces.gb/romaine/romaine.shtml
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

//===================https://www.shadertoy.com/view/wsXGWS======================
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

int quadratic(vec3 coeffs, out vec2 res) {
  return quadratic(coeffs[0],coeffs[1],coeffs[2],res);
}

void eval(float X, float A, float B, float C, float D,
          out float Q, out float Q1, out float B1,out float C2) {
  float q0 = A*X;
  B1 = q0+B;
  C2 = B1*X+C;
  Q1 = (q0+B1)*X + C2;
  Q = C2*X + D;
}

// Solve: Ax^3 + Bx^2 + Cx + D == 0
// Find one real root, then reduce to quadratic.
int cubic(float A, float B, float C, float D, out vec3 res) {
  float X,b1,c2;
  if (A == 0.0) {
    X = 1e8; A = B; b1 = C; c2 = D;
  } else if (D == 0.0) {
    X = 0.0; b1 = B; c2 = C;
  } else {
    X = -(B/A)/3.0;
    float t,r,s,q,dq,x0;
    eval(X,A,B,C,D,q,dq,b1,c2);
    t = q/A; r = pow(abs(t),1.0/3.0); s = sgn(t);
    t = -dq/A; if (t > 0.0) r = 1.324718*max(r,sqrt(t));
    x0 = X - s*r;
    if (x0 != X) {
      for (int i = 0; i < 6; i++) {
        X = x0;
        eval(X,A,B,C,D,q,dq,b1,c2);
        if (dq == 0.0) break;
        x0 -= (q/dq);
      }
      if (abs(A)*X*X > abs(D/X)) {
        c2 = -D/X; b1 = (c2 - C)/X;
      }
    }
  }
  res.x = X;
  return 1 + quadratic(A,b1,c2,res.yz);
}

int cubic(vec4 coeffs, out vec3 res) {
  float A = coeffs[0], B = coeffs[1], C = coeffs[2], D = coeffs[3];
  return cubic(A,B,C,D,res);
}

// Special wrapper for cubic function for solving quartic.
// Find largest real root of x**3 + a*x**2 + b*x + c
// Assume c < 0
float qcubic(float a, float b, float c) {
  // c is always <= 0, but may be very
  // small, in which case we return an
  // approximation. Never return < 0.
  //assert(c <= 0.0);
  if (c == 0.0) return 0.0;
  
  vec3 res;
  int nroots = cubic(1.0,a,b,c,res);
  if (nroots == 1) return res.x;
  else return max(res.x,max(res.y,res.z));
}

int quartic(vec4 coeffs, out vec4 res) {
  float c1 = coeffs[0];
  float c2 = coeffs[1];
  float c3 = coeffs[2];
  float c4 = coeffs[3];
  float alpha = 0.5*c1;
  float A = c2-alpha*alpha;
  float B = c3-alpha*A;
  float a,b,beta,psi;
  psi = qcubic(2.0*A-alpha*alpha, A*A+2.0*B*alpha-4.0*c4, -B*B);
  //assert(!isnan(psi));
  //assert(!isinf(psi));
  //assert(psi >= 0.0);
  a = sqrt(psi);
  beta = 0.5*(A + psi);
  if (psi <= 0.0) {
    b = sqrt(max(beta*beta-c4,0.0));
  } else {
    b = 0.5*a*(alpha-B/psi);
  }
  int resn = quadratic(1.0,alpha+a,beta+b,res.xy);
  vec2 tmp;
  if (quadratic(1.0,alpha-a,beta-b,tmp) != 0) { 
    res.zw = res.xy;
    res.xy = tmp;
    resn += 2;
  }
  return resn;
}

int quartic(float A, float B, float C, float D, float E, out vec4 roots) {
  int nroots;
  vec4 coeffs = vec4(B,C,D,E)/A;
  nroots = quartic(coeffs,roots);
  return nroots;
}

//https://www.shadertoy.com/view/wsXGWS

HIT giper3D(vec3 ro, vec3 rd, float t, float r)
{
    float a = ro.x;
    float b = rd.x;
    float c = ro.y;
    float d = rd.y;
    float e = ro.z;
    float f = rd.z;
    
    //https://github.com/flightphone/shaderty/blob/master/staples_polynomial.py
    //for generate this expression used python script staples_polynomial.py

    float a0 = 2.*a*a*c*c + 2.*c*c*e*e + 2.*a*a*e*e-1.*a*a*t*t-1.*c*c*t*t-1.*e*e*t*t-1.*t*t*t*t;
    float a1 = 4.*a*a*c*d + 4.*a*b*c*c + 4.*c*c*e*f + 4.*c*d*e*e + 4.*a*a*e*f + 4.*a*b*e*e-2.*a*b*t*t-2.*c*d*t*t-2.*e*f*t*t;
    float a2 = 2.*a*a*d*d + 8.*a*b*c*d + 2.*b*b*c*c + 2.*c*c*f*f + 8.*c*d*e*f + 2.*d*d*e*e + 2.*a*a*f*f + 8.*a*b*e*f + 2.*b*b*e*e-1.*b*b*t*t-1.*d*d*t*t-1.*f*f*t*t;
    float a3 = 4.*a*b*d*d + 4.*b*b*c*d + 4.*c*d*f*f + 4.*d*d*e*f + 4.*a*b*f*f + 4.*b*b*e*f;
    float a4 = 2.*b*b*d*d + 2.*d*d*f*f + 2.*b*b*f*f;

    //https://github.com/flightphone/shaderty/blob/master/staples_polynomial.py

    vec4 roots = vec4(dist_infin);
    int nroots = quartic(a4, a3, a2, a1, a0, roots); //cubic(a3, a2, a1, a0, roots);
    
    
    float dist = dist_infin;
    vec3 pos = vec3(0.0);
    vec3 nor = vec3(0.0);
    for (int i = 0; i < 4; i++)
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
        nor = vec3(0.+4.*pos.x*pos.y*pos.y+4.*pos.x*pos.z*pos.z-2.*t*t*pos.x, 0.+4.*pos.x*pos.x*pos.y+4.*pos.y*pos.z*pos.z-2.*t*t*pos.y, 0.+4.*pos.y*pos.y*pos.z+4.*pos.x*pos.x*pos.z-2.*t*t*pos.z);
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
    vec3 ro = vec3(0.0, 0.0, 6.5); // camera
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