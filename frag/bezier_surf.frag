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
#define TAU 6.28318530718
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

struct BEZIER2
{
    vec2 v0;
    vec2 v1;
    vec2 v2;
};

struct BEZIER_K
{
    vec2 a0;
    vec2 a1;
    vec2 a2;
};



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
BEZIER_K ck(BEZIER2 b)
{
    //t^2(P0-2P1+P2) + t(2P1-2P0) + P0
    BEZIER_K res;
    res.a2 = b.v0 - 2.0*b.v1 + b.v2;
    res.a1 = (2.0*b.v1 - 2.0*b.v0);
    res.a0 = b.v0;
    return res;
}    
//https://www.shadertoy.com/view/wsXGWS

float rand(float t, float n)
{
    float res =  fract(sin(t)*n);
    return res;
}
float getfi(float t, float n)
{
    float r = rand(t, n);
    return  r*PI/3.0 + PI/20.0;
}

BEZIER2 tileCurve(float t0, float n)
{
    float t1 = t0 + 1.0;
    float a = getfi(t1, n) ;
    float b = getfi(t0, n) ;
    float x = tan(b)/(tan(a) + tan(b));
    float h = x*tan(a);
    x = 1.0 - x;
    float y = 0.5 + h;
    if (mod(t1, 2.0) == 0.0)
        y = 0.5 - h;

    BEZIER2 res;
    res.v0 = vec2(0.,0.5);
    res.v1 = vec2(x,y);
    res.v2 = vec2(1.0,0.5);
    return res;
    
}


HIT giper3D(vec3 ro, vec3 rd, BEZIER2 xb, BEZIER2 yb, float xstart, float ystart)
{
    float r = ro.x;
    float R = rd.x;
    float s = ro.y;
    float S = rd.y;
    float t = ro.z;
    float T = rd.z;
    //t^2(P0-2P1+P2) + t(2P1-2P0) + P0
    BEZIER_K ko = ck(xb);
    float a = ko.a2.x;
    float b = ko.a1.x;
    float c = ko.a0.x;
    float g = ko.a2.y;
    float h = ko.a1.y;
    float l = ko.a0.y;

    ko = ck(yb);
    float d = ko.a2.x;
    float e = ko.a1.x;
    float f = ko.a0.x;
    float i = ko.a2.y;
    float k = ko.a1.y;
    l += ko.a0.y;
    float A = g/a;
    float B = i/d;
    float dd = (h - A*b);
    float F = 1.0 / dd;
    float C = -A/dd;
    float D = -B/dd;
    float G = -(k - B*e)/dd;
    float H = -(l - A*c - B*f)/dd;
    float K = a*G*G/d;

    /*
    x = au2 + bu + c
    y = dv2 + ev + f
    z = gu2 + hu + iv2 + kv + l

    A = g/a
    B = i/d
    */
    
    //https://github.com/flightphone/shaderty/blob/master/staples_polynomial.py
    //for generate this expression used python script staples_polynomial.py
    float a0 = 1.*d*r*r;
a0 = a0 -2.*K*d*r*s;
a0 = a0 -2.*C*C*a*d*r*r*r;
a0 = a0 -4.*C*D*a*d*r*r*s;
a0 = a0 -4.*C*F*a*d*r*r*t;
a0 = a0 -4.*C*H*a*d*r*r;
a0 = a0 -2.*D*D*a*d*r*s*s;
a0 = a0 -4.*D*F*a*d*r*s*t;
a0 = a0 -4.*D*H*a*d*r*s;
a0 = a0 -2.*F*F*a*d*r*t*t;
a0 = a0 -4.*F*H*a*d*r*t;
a0 = a0 -2.*H*H*a*d*r;
a0 = a0 -2.*C*b*d*r*r;
a0 = a0 -2.*D*b*d*r*s;
a0 = a0 -2.*F*b*d*r*t;
a0 = a0 -2.*H*b*d*r;
a0 = a0  + 2.*K*d*f*r;
a0 = a0  + 1.*K*K*d*s*s;
a0 = a0  + 2.*C*C*K*a*d*r*r*s;
a0 = a0  + 4.*C*D*K*a*d*r*s*s;
a0 = a0  + 4.*C*F*K*a*d*r*s*t;
a0 = a0  + 4.*C*H*K*a*d*r*s;
a0 = a0  + 2.*D*D*K*a*d*s*s*s;
a0 = a0  + 4.*D*F*K*a*d*s*s*t;
a0 = a0  + 4.*D*H*K*a*d*s*s;
a0 = a0  + 2.*F*F*K*a*d*s*t*t;
a0 = a0  + 4.*F*H*K*a*d*s*t;
a0 = a0  + 2.*H*H*K*a*d*s;
a0 = a0  + 2.*C*K*b*d*r*s;
a0 = a0  + 2.*D*K*b*d*s*s;
a0 = a0  + 2.*F*K*b*d*s*t;
a0 = a0  + 2.*H*K*b*d*s;
a0 = a0 -2.*K*K*d*f*s;
a0 = a0  + 1.*C*C*C*C*a*a*d*r*r*r*r;
a0 = a0  + 4.*C*C*C*D*a*a*d*r*r*r*s;
a0 = a0  + 4.*C*C*C*F*a*a*d*r*r*r*t;
a0 = a0  + 4.*C*C*C*H*a*a*d*r*r*r;
a0 = a0  + 6.*C*C*D*D*a*a*d*r*r*s*s;
a0 = a0  + 12.*C*C*D*F*a*a*d*r*r*s*t;
a0 = a0  + 12.*C*C*D*H*a*a*d*r*r*s;
a0 = a0  + 6.*C*C*F*F*a*a*d*r*r*t*t;
a0 = a0  + 12.*C*C*F*H*a*a*d*r*r*t;
a0 = a0  + 6.*C*C*H*H*a*a*d*r*r;
a0 = a0  + 2.*C*C*C*a*b*d*r*r*r;
a0 = a0  + 6.*C*C*D*a*b*d*r*r*s;
a0 = a0  + 6.*C*C*F*a*b*d*r*r*t;
a0 = a0  + 6.*C*C*H*a*b*d*r*r;
a0 = a0 -2.*C*C*K*a*d*f*r*r;
a0 = a0  + 4.*C*D*D*D*a*a*d*r*s*s*s;
a0 = a0  + 12.*C*D*D*F*a*a*d*r*s*s*t;
a0 = a0  + 12.*C*D*D*H*a*a*d*r*s*s;
a0 = a0  + 12.*C*D*F*F*a*a*d*r*s*t*t;
a0 = a0  + 24.*C*D*F*H*a*a*d*r*s*t;
a0 = a0  + 12.*C*D*H*H*a*a*d*r*s;
a0 = a0  + 6.*C*D*D*a*b*d*r*s*s;
a0 = a0  + 12.*C*D*F*a*b*d*r*s*t;
a0 = a0  + 12.*C*D*H*a*b*d*r*s;
a0 = a0 -4.*C*D*K*a*d*f*r*s;
a0 = a0  + 4.*C*F*F*F*a*a*d*r*t*t*t;
a0 = a0  + 12.*C*F*F*H*a*a*d*r*t*t;
a0 = a0  + 12.*C*F*H*H*a*a*d*r*t;
a0 = a0  + 6.*C*F*F*a*b*d*r*t*t;
a0 = a0  + 12.*C*F*H*a*b*d*r*t;
a0 = a0 -4.*C*F*K*a*d*f*r*t;
a0 = a0  + 4.*C*H*H*H*a*a*d*r;
a0 = a0  + 6.*C*H*H*a*b*d*r;
a0 = a0 -4.*C*H*K*a*d*f*r;
a0 = a0  + 1.*D*D*D*D*a*a*d*s*s*s*s;
a0 = a0  + 4.*D*D*D*F*a*a*d*s*s*s*t;
a0 = a0  + 4.*D*D*D*H*a*a*d*s*s*s;
a0 = a0  + 6.*D*D*F*F*a*a*d*s*s*t*t;
a0 = a0  + 12.*D*D*F*H*a*a*d*s*s*t;
a0 = a0  + 6.*D*D*H*H*a*a*d*s*s;
a0 = a0  + 2.*D*D*D*a*b*d*s*s*s;
a0 = a0  + 6.*D*D*F*a*b*d*s*s*t;
a0 = a0  + 6.*D*D*H*a*b*d*s*s;
a0 = a0 -2.*D*D*K*a*d*f*s*s;
a0 = a0  + 4.*D*F*F*F*a*a*d*s*t*t*t;
a0 = a0  + 12.*D*F*F*H*a*a*d*s*t*t;
a0 = a0  + 12.*D*F*H*H*a*a*d*s*t;
a0 = a0  + 6.*D*F*F*a*b*d*s*t*t;
a0 = a0  + 12.*D*F*H*a*b*d*s*t;
a0 = a0 -4.*D*F*K*a*d*f*s*t;
a0 = a0  + 4.*D*H*H*H*a*a*d*s;
a0 = a0  + 6.*D*H*H*a*b*d*s;
a0 = a0 -4.*D*H*K*a*d*f*s;
a0 = a0  + 1.*F*F*F*F*a*a*d*t*t*t*t;
a0 = a0  + 4.*F*F*F*H*a*a*d*t*t*t;
a0 = a0  + 6.*F*F*H*H*a*a*d*t*t;
a0 = a0  + 2.*F*F*F*a*b*d*t*t*t;
a0 = a0  + 6.*F*F*H*a*b*d*t*t;
a0 = a0 -2.*F*F*K*a*d*f*t*t;
a0 = a0  + 4.*F*H*H*H*a*a*d*t;
a0 = a0  + 6.*F*H*H*a*b*d*t;
a0 = a0 -4.*F*H*K*a*d*f*t;
a0 = a0  + 1.*H*H*H*H*a*a*d;
a0 = a0  + 2.*H*H*H*a*b*d;
a0 = a0 -2.*H*H*K*a*d*f;
a0 = a0  + 1.*C*C*b*b*d*r*r;
a0 = a0  + 2.*C*D*b*b*d*r*s;
a0 = a0  + 2.*C*F*b*b*d*r*t;
a0 = a0  + 2.*C*H*b*b*d*r;
a0 = a0 -2.*C*K*b*d*f*r;
a0 = a0  + 1.*D*D*b*b*d*s*s;
a0 = a0  + 2.*D*F*b*b*d*s*t;
a0 = a0  + 2.*D*H*b*b*d*s;
a0 = a0 -2.*D*K*b*d*f*s;
a0 = a0  + 1.*F*F*b*b*d*t*t;
a0 = a0  + 2.*F*H*b*b*d*t;
a0 = a0 -2.*F*K*b*d*f*t;
a0 = a0  + 1.*H*H*b*b*d;
a0 = a0 -2.*H*K*b*d*f;
a0 = a0  + 1.*K*K*d*f*f;
a0 = a0  + 2.*C*G*a*e*r*r;
a0 = a0  + 2.*D*G*a*e*r*s;
a0 = a0  + 2.*F*G*a*e*r*t;
a0 = a0  + 2.*G*H*a*e*r;
a0 = a0  + 1.*G*b*e*r;
a0 = a0 -1.*K*e*e*r;
a0 = a0  + 2.*C*G*K*a*e*r*s;
a0 = a0  + 2.*D*G*K*a*e*s*s;
a0 = a0  + 2.*F*G*K*a*e*s*t;
a0 = a0  + 2.*G*H*K*a*e*s;
a0 = a0  + 1.*G*K*b*e*s;
a0 = a0 -2.*C*C*C*G*a*a*e*r*r*r;
a0 = a0 -6.*C*C*D*G*a*a*e*r*r*s;
a0 = a0 -6.*C*C*F*G*a*a*e*r*r*t;
a0 = a0 -6.*C*C*G*H*a*a*e*r*r;
a0 = a0 -3.*C*C*G*a*b*e*r*r;
a0 = a0  + 1.*C*C*K*a*e*e*r*r;
a0 = a0 -6.*C*D*D*G*a*a*e*r*s*s;
a0 = a0 -12.*C*D*F*G*a*a*e*r*s*t;
a0 = a0 -12.*C*D*G*H*a*a*e*r*s;
a0 = a0 -6.*C*D*G*a*b*e*r*s;
a0 = a0  + 2.*C*D*K*a*e*e*r*s;
a0 = a0 -6.*C*F*F*G*a*a*e*r*t*t;
a0 = a0 -12.*C*F*G*H*a*a*e*r*t;
a0 = a0 -6.*C*F*G*a*b*e*r*t;
a0 = a0  + 2.*C*F*K*a*e*e*r*t;
a0 = a0 -6.*C*G*H*H*a*a*e*r;
a0 = a0 -6.*C*G*H*a*b*e*r;
a0 = a0  + 2.*C*H*K*a*e*e*r;
a0 = a0 -2.*D*D*D*G*a*a*e*s*s*s;
a0 = a0 -6.*D*D*F*G*a*a*e*s*s*t;
a0 = a0 -6.*D*D*G*H*a*a*e*s*s;
a0 = a0 -3.*D*D*G*a*b*e*s*s;
a0 = a0  + 1.*D*D*K*a*e*e*s*s;
a0 = a0 -6.*D*F*F*G*a*a*e*s*t*t;
a0 = a0 -12.*D*F*G*H*a*a*e*s*t;
a0 = a0 -6.*D*F*G*a*b*e*s*t;
a0 = a0  + 2.*D*F*K*a*e*e*s*t;
a0 = a0 -6.*D*G*H*H*a*a*e*s;
a0 = a0 -6.*D*G*H*a*b*e*s;
a0 = a0  + 2.*D*H*K*a*e*e*s;
a0 = a0 -2.*F*F*F*G*a*a*e*t*t*t;
a0 = a0 -6.*F*F*G*H*a*a*e*t*t;
a0 = a0 -3.*F*F*G*a*b*e*t*t;
a0 = a0  + 1.*F*F*K*a*e*e*t*t;
a0 = a0 -6.*F*G*H*H*a*a*e*t;
a0 = a0 -6.*F*G*H*a*b*e*t;
a0 = a0  + 2.*F*H*K*a*e*e*t;
a0 = a0 -2.*G*H*H*H*a*a*e;
a0 = a0 -3.*G*H*H*a*b*e;
a0 = a0  + 1.*H*H*K*a*e*e;
a0 = a0 -1.*C*G*b*b*e*r;
a0 = a0  + 1.*C*K*b*e*e*r;
a0 = a0 -1.*D*G*b*b*e*s;
a0 = a0  + 1.*D*K*b*e*e*s;
a0 = a0 -1.*F*G*b*b*e*t;
a0 = a0  + 1.*F*K*b*e*e*t;
a0 = a0 -1.*G*H*b*b*e;
a0 = a0  + 1.*H*K*b*e*e;
a0 = a0 -2.*C*G*K*a*e*f*r;
a0 = a0 -2.*D*G*K*a*e*f*s;
a0 = a0 -2.*F*G*K*a*e*f*t;
a0 = a0 -2.*G*H*K*a*e*f;
a0 = a0 -1.*G*K*b*e*f;
a0 = a0  + 4.*C*C*G*G*a*a*f*r*r;
a0 = a0  + 8.*C*D*G*G*a*a*f*r*s;
a0 = a0  + 8.*C*F*G*G*a*a*f*r*t;
a0 = a0  + 8.*C*G*G*H*a*a*f*r;
a0 = a0  + 4.*C*G*G*a*b*f*r;
a0 = a0  + 4.*D*D*G*G*a*a*f*s*s;
a0 = a0  + 8.*D*F*G*G*a*a*f*s*t;
a0 = a0  + 8.*D*G*G*H*a*a*f*s;
a0 = a0  + 4.*D*G*G*a*b*f*s;
a0 = a0  + 4.*F*F*G*G*a*a*f*t*t;
a0 = a0  + 8.*F*G*G*H*a*a*f*t;
a0 = a0  + 4.*F*G*G*a*b*f*t;
a0 = a0  + 4.*G*G*H*H*a*a*f;
a0 = a0  + 4.*G*G*H*a*b*f;
a0 = a0  + 1.*G*G*b*b*f;
a0 = a0 -4.*C*C*G*G*a*a*r*r*s;
a0 = a0 -8.*C*D*G*G*a*a*r*s*s;
a0 = a0 -8.*C*F*G*G*a*a*r*s*t;
a0 = a0 -8.*C*G*G*H*a*a*r*s;
a0 = a0 -4.*C*G*G*a*b*r*s;
a0 = a0 -4.*D*D*G*G*a*a*s*s*s;
a0 = a0 -8.*D*F*G*G*a*a*s*s*t;
a0 = a0 -8.*D*G*G*H*a*a*s*s;
a0 = a0 -4.*D*G*G*a*b*s*s;
a0 = a0 -4.*F*F*G*G*a*a*s*t*t;
a0 = a0 -8.*F*G*G*H*a*a*s*t;
a0 = a0 -4.*F*G*G*a*b*s*t;
a0 = a0 -4.*G*G*H*H*a*a*s;
a0 = a0 -4.*G*G*H*a*b*s;
a0 = a0 -1.*G*G*b*b*s;
a0 = a0 ;
float a1 = 2.*R*d*r;
a1 = a1 -2.*K*S*d*r;
a1 = a1 -6.*C*C*R*a*d*r*r;
a1 = a1 -4.*C*D*S*a*d*r*r;
a1 = a1 -8.*C*D*R*a*d*r*s;
a1 = a1 -4.*C*F*T*a*d*r*r;
a1 = a1 -8.*C*F*R*a*d*r*t;
a1 = a1 -8.*C*H*R*a*d*r;
a1 = a1 -4.*D*D*S*a*d*r*s;
a1 = a1 -4.*D*F*T*a*d*r*s;
a1 = a1 -4.*D*F*S*a*d*r*t;
a1 = a1 -4.*D*H*S*a*d*r;
a1 = a1 -4.*F*F*T*a*d*r*t;
a1 = a1 -4.*F*H*T*a*d*r;
a1 = a1 -4.*C*R*b*d*r;
a1 = a1 -2.*D*S*b*d*r;
a1 = a1 -2.*F*T*b*d*r;
a1 = a1 -2.*K*R*d*s;
a1 = a1  + 2.*K*K*S*d*s;
a1 = a1  + 4.*C*C*K*R*a*d*r*s;
a1 = a1  + 8.*C*D*K*S*a*d*r*s;
a1 = a1  + 4.*C*D*K*R*a*d*s*s;
a1 = a1  + 4.*C*F*K*T*a*d*r*s;
a1 = a1  + 4.*C*F*K*R*a*d*s*t;
a1 = a1  + 4.*C*H*K*R*a*d*s;
a1 = a1  + 6.*D*D*K*S*a*d*s*s;
a1 = a1  + 4.*D*F*K*T*a*d*s*s;
a1 = a1  + 8.*D*F*K*S*a*d*s*t;
a1 = a1  + 8.*D*H*K*S*a*d*s;
a1 = a1  + 4.*F*F*K*T*a*d*s*t;
a1 = a1  + 4.*F*H*K*T*a*d*s;
a1 = a1  + 2.*C*K*R*b*d*s;
a1 = a1  + 4.*D*K*S*b*d*s;
a1 = a1  + 2.*F*K*T*b*d*s;
a1 = a1  + 2.*C*C*K*S*a*d*r*r;
a1 = a1  + 4.*C*C*C*C*R*a*a*d*r*r*r;
a1 = a1  + 4.*C*C*C*D*S*a*a*d*r*r*r;
a1 = a1  + 12.*C*C*C*D*R*a*a*d*r*r*s;
a1 = a1  + 4.*C*C*C*F*T*a*a*d*r*r*r;
a1 = a1  + 12.*C*C*C*F*R*a*a*d*r*r*t;
a1 = a1  + 12.*C*C*C*H*R*a*a*d*r*r;
a1 = a1  + 12.*C*C*D*D*S*a*a*d*r*r*s;
a1 = a1  + 12.*C*C*D*F*T*a*a*d*r*r*s;
a1 = a1  + 12.*C*C*D*F*S*a*a*d*r*r*t;
a1 = a1  + 12.*C*C*D*H*S*a*a*d*r*r;
a1 = a1  + 12.*C*C*F*F*T*a*a*d*r*r*t;
a1 = a1  + 12.*C*C*F*H*T*a*a*d*r*r;
a1 = a1  + 6.*C*C*C*R*a*b*d*r*r;
a1 = a1  + 6.*C*C*D*S*a*b*d*r*r;
a1 = a1  + 6.*C*C*F*T*a*b*d*r*r;
a1 = a1  + 12.*C*C*D*D*R*a*a*d*r*s*s;
a1 = a1  + 24.*C*C*D*F*R*a*a*d*r*s*t;
a1 = a1  + 24.*C*C*D*H*R*a*a*d*r*s;
a1 = a1  + 12.*C*D*D*D*S*a*a*d*r*s*s;
a1 = a1  + 12.*C*D*D*F*T*a*a*d*r*s*s;
a1 = a1  + 24.*C*D*D*F*S*a*a*d*r*s*t;
a1 = a1  + 24.*C*D*D*H*S*a*a*d*r*s;
a1 = a1  + 24.*C*D*F*F*T*a*a*d*r*s*t;
a1 = a1  + 24.*C*D*F*H*T*a*a*d*r*s;
a1 = a1  + 12.*C*C*D*R*a*b*d*r*s;
a1 = a1  + 12.*C*D*D*S*a*b*d*r*s;
a1 = a1  + 12.*C*D*F*T*a*b*d*r*s;
a1 = a1  + 4.*C*F*K*S*a*d*r*t;
a1 = a1  + 12.*C*C*F*F*R*a*a*d*r*t*t;
a1 = a1  + 24.*C*C*F*H*R*a*a*d*r*t;
a1 = a1  + 12.*C*D*F*F*S*a*a*d*r*t*t;
a1 = a1  + 24.*C*D*F*H*S*a*a*d*r*t;
a1 = a1  + 12.*C*F*F*F*T*a*a*d*r*t*t;
a1 = a1  + 24.*C*F*F*H*T*a*a*d*r*t;
a1 = a1  + 12.*C*C*F*R*a*b*d*r*t;
a1 = a1  + 12.*C*D*F*S*a*b*d*r*t;
a1 = a1  + 12.*C*F*F*T*a*b*d*r*t;
a1 = a1  + 4.*C*H*K*S*a*d*r;
a1 = a1  + 12.*C*C*H*H*R*a*a*d*r;
a1 = a1  + 12.*C*D*H*H*S*a*a*d*r;
a1 = a1  + 12.*C*F*H*H*T*a*a*d*r;
a1 = a1  + 12.*C*C*H*R*a*b*d*r;
a1 = a1  + 12.*C*D*H*S*a*b*d*r;
a1 = a1  + 12.*C*F*H*T*a*b*d*r;
a1 = a1 -2.*D*D*R*a*d*s*s;
a1 = a1  + 4.*C*D*D*D*R*a*a*d*s*s*s;
a1 = a1  + 12.*C*D*D*F*R*a*a*d*s*s*t;
a1 = a1  + 12.*C*D*D*H*R*a*a*d*s*s;
a1 = a1  + 4.*D*D*D*D*S*a*a*d*s*s*s;
a1 = a1  + 4.*D*D*D*F*T*a*a*d*s*s*s;
a1 = a1  + 12.*D*D*D*F*S*a*a*d*s*s*t;
a1 = a1  + 12.*D*D*D*H*S*a*a*d*s*s;
a1 = a1  + 12.*D*D*F*F*T*a*a*d*s*s*t;
a1 = a1  + 12.*D*D*F*H*T*a*a*d*s*s;
a1 = a1  + 6.*C*D*D*R*a*b*d*s*s;
a1 = a1  + 6.*D*D*D*S*a*b*d*s*s;
a1 = a1  + 6.*D*D*F*T*a*b*d*s*s;
a1 = a1 -4.*D*F*R*a*d*s*t;
a1 = a1  + 12.*C*D*F*F*R*a*a*d*s*t*t;
a1 = a1  + 24.*C*D*F*H*R*a*a*d*s*t;
a1 = a1  + 12.*D*D*F*F*S*a*a*d*s*t*t;
a1 = a1  + 24.*D*D*F*H*S*a*a*d*s*t;
a1 = a1  + 12.*D*F*F*F*T*a*a*d*s*t*t;
a1 = a1  + 24.*D*F*F*H*T*a*a*d*s*t;
a1 = a1  + 12.*C*D*F*R*a*b*d*s*t;
a1 = a1  + 12.*D*D*F*S*a*b*d*s*t;
a1 = a1  + 12.*D*F*F*T*a*b*d*s*t;
a1 = a1 -4.*D*H*R*a*d*s;
a1 = a1  + 12.*C*D*H*H*R*a*a*d*s;
a1 = a1  + 12.*D*D*H*H*S*a*a*d*s;
a1 = a1  + 12.*D*F*H*H*T*a*a*d*s;
a1 = a1  + 12.*C*D*H*R*a*b*d*s;
a1 = a1  + 12.*D*D*H*S*a*b*d*s;
a1 = a1  + 12.*D*F*H*T*a*b*d*s;
a1 = a1 -2.*F*F*R*a*d*t*t;
a1 = a1  + 2.*F*F*K*S*a*d*t*t;
a1 = a1  + 4.*C*F*F*F*R*a*a*d*t*t*t;
a1 = a1  + 12.*C*F*F*H*R*a*a*d*t*t;
a1 = a1  + 4.*D*F*F*F*S*a*a*d*t*t*t;
a1 = a1  + 12.*D*F*F*H*S*a*a*d*t*t;
a1 = a1  + 4.*F*F*F*F*T*a*a*d*t*t*t;
a1 = a1  + 12.*F*F*F*H*T*a*a*d*t*t;
a1 = a1  + 6.*C*F*F*R*a*b*d*t*t;
a1 = a1  + 6.*D*F*F*S*a*b*d*t*t;
a1 = a1  + 6.*F*F*F*T*a*b*d*t*t;
a1 = a1 -4.*F*H*R*a*d*t;
a1 = a1  + 4.*F*H*K*S*a*d*t;
a1 = a1  + 12.*C*F*H*H*R*a*a*d*t;
a1 = a1  + 12.*D*F*H*H*S*a*a*d*t;
a1 = a1  + 12.*F*F*H*H*T*a*a*d*t;
a1 = a1  + 12.*C*F*H*R*a*b*d*t;
a1 = a1  + 12.*D*F*H*S*a*b*d*t;
a1 = a1  + 12.*F*F*H*T*a*b*d*t;
a1 = a1 -2.*H*H*R*a*d;
a1 = a1  + 2.*H*H*K*S*a*d;
a1 = a1  + 4.*C*H*H*H*R*a*a*d;
a1 = a1  + 4.*D*H*H*H*S*a*a*d;
a1 = a1  + 4.*F*H*H*H*T*a*a*d;
a1 = a1  + 6.*C*H*H*R*a*b*d;
a1 = a1  + 6.*D*H*H*S*a*b*d;
a1 = a1  + 6.*F*H*H*T*a*b*d;
a1 = a1  + 2.*C*K*S*b*d*r;
a1 = a1  + 2.*C*C*R*b*b*d*r;
a1 = a1  + 2.*C*D*S*b*b*d*r;
a1 = a1  + 2.*C*F*T*b*b*d*r;
a1 = a1 -2.*D*R*b*d*s;
a1 = a1  + 2.*C*D*R*b*b*d*s;
a1 = a1  + 2.*D*D*S*b*b*d*s;
a1 = a1  + 2.*D*F*T*b*b*d*s;
a1 = a1 -2.*F*R*b*d*t;
a1 = a1  + 2.*F*K*S*b*d*t;
a1 = a1  + 2.*C*F*R*b*b*d*t;
a1 = a1  + 2.*D*F*S*b*b*d*t;
a1 = a1  + 2.*F*F*T*b*b*d*t;
a1 = a1 -2.*H*R*b*d;
a1 = a1  + 2.*H*K*S*b*d;
a1 = a1  + 2.*C*H*R*b*b*d;
a1 = a1  + 2.*D*H*S*b*b*d;
a1 = a1  + 2.*F*H*T*b*b*d;
a1 = a1  + 2.*K*R*d*f;
a1 = a1 -2.*K*K*S*d*f;
a1 = a1 -4.*C*C*K*R*a*d*f*r;
a1 = a1 -4.*C*D*K*S*a*d*f*r;
a1 = a1 -4.*C*D*K*R*a*d*f*s;
a1 = a1 -4.*C*F*K*T*a*d*f*r;
a1 = a1 -4.*C*F*K*R*a*d*f*t;
a1 = a1 -4.*C*H*K*R*a*d*f;
a1 = a1 -4.*D*D*K*S*a*d*f*s;
a1 = a1 -4.*D*F*K*T*a*d*f*s;
a1 = a1 -4.*D*F*K*S*a*d*f*t;
a1 = a1 -4.*D*H*K*S*a*d*f;
a1 = a1 -4.*F*F*K*T*a*d*f*t;
a1 = a1 -4.*F*H*K*T*a*d*f;
a1 = a1 -2.*C*K*R*b*d*f;
a1 = a1 -2.*D*K*S*b*d*f;
a1 = a1 -2.*F*K*T*b*d*f;
a1 = a1  + 4.*C*G*R*a*e*r;
a1 = a1  + 2.*D*G*S*a*e*r;
a1 = a1  + 2.*F*G*T*a*e*r;
a1 = a1  + 2.*C*G*K*R*a*e*s;
a1 = a1  + 4.*D*G*K*S*a*e*s;
a1 = a1  + 2.*F*G*K*T*a*e*s;
a1 = a1 -6.*C*C*C*G*R*a*a*e*r*r;
a1 = a1 -6.*C*C*D*G*S*a*a*e*r*r;
a1 = a1 -6.*C*C*F*G*T*a*a*e*r*r;
a1 = a1 -12.*C*C*D*G*R*a*a*e*r*s;
a1 = a1 -12.*C*D*D*G*S*a*a*e*r*s;
a1 = a1 -12.*C*D*F*G*T*a*a*e*r*s;
a1 = a1 -12.*C*C*F*G*R*a*a*e*r*t;
a1 = a1 -12.*C*D*F*G*S*a*a*e*r*t;
a1 = a1 -12.*C*F*F*G*T*a*a*e*r*t;
a1 = a1 -12.*C*C*G*H*R*a*a*e*r;
a1 = a1 -12.*C*D*G*H*S*a*a*e*r;
a1 = a1 -12.*C*F*G*H*T*a*a*e*r;
a1 = a1 -6.*C*D*D*G*R*a*a*e*s*s;
a1 = a1 -6.*D*D*D*G*S*a*a*e*s*s;
a1 = a1 -6.*D*D*F*G*T*a*a*e*s*s;
a1 = a1 -12.*C*D*F*G*R*a*a*e*s*t;
a1 = a1 -12.*D*D*F*G*S*a*a*e*s*t;
a1 = a1 -12.*D*F*F*G*T*a*a*e*s*t;
a1 = a1 -12.*C*D*G*H*R*a*a*e*s;
a1 = a1 -12.*D*D*G*H*S*a*a*e*s;
a1 = a1 -12.*D*F*G*H*T*a*a*e*s;
a1 = a1 -6.*C*F*F*G*R*a*a*e*t*t;
a1 = a1 -6.*D*F*F*G*S*a*a*e*t*t;
a1 = a1 -6.*F*F*F*G*T*a*a*e*t*t;
a1 = a1 -12.*C*F*G*H*R*a*a*e*t;
a1 = a1 -12.*D*F*G*H*S*a*a*e*t;
a1 = a1 -12.*F*F*G*H*T*a*a*e*t;
a1 = a1 -6.*C*G*H*H*R*a*a*e;
a1 = a1 -6.*D*G*H*H*S*a*a*e;
a1 = a1 -6.*F*G*H*H*T*a*a*e;
a1 = a1 -6.*C*C*G*R*a*b*e*r;
a1 = a1 -6.*C*D*G*S*a*b*e*r;
a1 = a1 -6.*C*F*G*T*a*b*e*r;
a1 = a1 -6.*C*D*G*R*a*b*e*s;
a1 = a1 -6.*D*D*G*S*a*b*e*s;
a1 = a1 -6.*D*F*G*T*a*b*e*s;
a1 = a1 -6.*C*F*G*R*a*b*e*t;
a1 = a1 -6.*D*F*G*S*a*b*e*t;
a1 = a1 -6.*F*F*G*T*a*b*e*t;
a1 = a1 -6.*C*G*H*R*a*b*e;
a1 = a1 -6.*D*G*H*S*a*b*e;
a1 = a1 -6.*F*G*H*T*a*b*e;
a1 = a1 -2.*C*G*K*R*a*e*f;
a1 = a1 -2.*D*G*K*S*a*e*f;
a1 = a1 -2.*F*G*K*T*a*e*f;
a1 = a1  + 2.*D*G*R*a*e*s;
a1 = a1  + 2.*F*G*R*a*e*t;
a1 = a1  + 2.*G*H*R*a*e;
a1 = a1  + 1.*G*R*b*e;
a1 = a1 -1.*K*R*e*e;
a1 = a1  + 2.*C*G*K*S*a*e*r;
a1 = a1  + 2.*F*G*K*S*a*e*t;
a1 = a1  + 2.*G*H*K*S*a*e;
a1 = a1  + 1.*G*K*S*b*e;
a1 = a1  + 2.*C*C*K*R*a*e*e*r;
a1 = a1  + 2.*C*D*K*S*a*e*e*r;
a1 = a1  + 2.*C*D*K*R*a*e*e*s;
a1 = a1  + 2.*C*F*K*T*a*e*e*r;
a1 = a1  + 2.*C*F*K*R*a*e*e*t;
a1 = a1  + 2.*C*H*K*R*a*e*e;
a1 = a1  + 2.*D*D*K*S*a*e*e*s;
a1 = a1  + 2.*D*F*K*T*a*e*e*s;
a1 = a1  + 2.*D*F*K*S*a*e*e*t;
a1 = a1  + 2.*D*H*K*S*a*e*e;
a1 = a1  + 2.*F*F*K*T*a*e*e*t;
a1 = a1  + 2.*F*H*K*T*a*e*e;
a1 = a1 -1.*C*G*R*b*b*e;
a1 = a1  + 1.*C*K*R*b*e*e;
a1 = a1 -1.*D*G*S*b*b*e;
a1 = a1  + 1.*D*K*S*b*e*e;
a1 = a1 -1.*F*G*T*b*b*e;
a1 = a1  + 1.*F*K*T*b*e*e;
a1 = a1  + 8.*C*C*G*G*R*a*a*f*r;
a1 = a1  + 8.*C*D*G*G*S*a*a*f*r;
a1 = a1  + 8.*C*F*G*G*T*a*a*f*r;
a1 = a1  + 8.*C*D*G*G*R*a*a*f*s;
a1 = a1  + 8.*D*D*G*G*S*a*a*f*s;
a1 = a1  + 8.*D*F*G*G*T*a*a*f*s;
a1 = a1  + 8.*C*F*G*G*R*a*a*f*t;
a1 = a1  + 8.*D*F*G*G*S*a*a*f*t;
a1 = a1  + 8.*F*F*G*G*T*a*a*f*t;
a1 = a1  + 8.*C*G*G*H*R*a*a*f;
a1 = a1  + 8.*D*G*G*H*S*a*a*f;
a1 = a1  + 8.*F*G*G*H*T*a*a*f;
a1 = a1  + 4.*C*G*G*R*a*b*f;
a1 = a1  + 4.*D*G*G*S*a*b*f;
a1 = a1  + 4.*F*G*G*T*a*b*f;
a1 = a1 -8.*C*C*G*G*R*a*a*r*s;
a1 = a1 -16.*C*D*G*G*S*a*a*r*s;
a1 = a1 -8.*C*F*G*G*T*a*a*r*s;
a1 = a1 -8.*C*D*G*G*R*a*a*s*s;
a1 = a1 -12.*D*D*G*G*S*a*a*s*s;
a1 = a1 -8.*D*F*G*G*T*a*a*s*s;
a1 = a1 -8.*C*F*G*G*R*a*a*s*t;
a1 = a1 -16.*D*F*G*G*S*a*a*s*t;
a1 = a1 -8.*F*F*G*G*T*a*a*s*t;
a1 = a1 -8.*C*G*G*H*R*a*a*s;
a1 = a1 -16.*D*G*G*H*S*a*a*s;
a1 = a1 -8.*F*G*G*H*T*a*a*s;
a1 = a1 -4.*C*G*G*R*a*b*s;
a1 = a1 -8.*D*G*G*S*a*b*s;
a1 = a1 -4.*F*G*G*T*a*b*s;
a1 = a1 -4.*C*C*G*G*S*a*a*r*r;
a1 = a1 -8.*C*F*G*G*S*a*a*r*t;
a1 = a1 -8.*C*G*G*H*S*a*a*r;
a1 = a1 -4.*C*G*G*S*a*b*r;
a1 = a1 -4.*F*F*G*G*S*a*a*t*t;
a1 = a1 -8.*F*G*G*H*S*a*a*t;
a1 = a1 -4.*F*G*G*S*a*b*t;
a1 = a1 -4.*G*G*H*H*S*a*a;
a1 = a1 -4.*G*G*H*S*a*b;
a1 = a1 -1.*G*G*S*b*b;
a1 = a1 ;
float a2 = -6.*C*C*R*R*a*d*r;
a2 = a2 -8.*C*D*R*S*a*d*r;
a2 = a2 -8.*C*F*R*T*a*d*r;
a2 = a2 -2.*D*D*S*S*a*d*r;
a2 = a2 -4.*D*F*S*T*a*d*r;
a2 = a2 -2.*F*F*T*T*a*d*r;
a2 = a2  + 2.*C*C*K*R*R*a*d*s;
a2 = a2  + 8.*C*D*K*R*S*a*d*s;
a2 = a2  + 4.*C*F*K*R*T*a*d*s;
a2 = a2  + 6.*D*D*K*S*S*a*d*s;
a2 = a2  + 8.*D*F*K*S*T*a*d*s;
a2 = a2  + 2.*F*F*K*T*T*a*d*s;
a2 = a2  + 6.*C*C*C*C*R*R*a*a*d*r*r;
a2 = a2  + 12.*C*C*C*D*R*S*a*a*d*r*r;
a2 = a2  + 12.*C*C*C*F*R*T*a*a*d*r*r;
a2 = a2  + 6.*C*C*D*D*S*S*a*a*d*r*r;
a2 = a2  + 12.*C*C*D*F*S*T*a*a*d*r*r;
a2 = a2  + 6.*C*C*F*F*T*T*a*a*d*r*r;
a2 = a2  + 12.*C*C*C*D*R*R*a*a*d*r*s;
a2 = a2  + 24.*C*C*D*D*R*S*a*a*d*r*s;
a2 = a2  + 24.*C*C*D*F*R*T*a*a*d*r*s;
a2 = a2  + 12.*C*D*D*D*S*S*a*a*d*r*s;
a2 = a2  + 24.*C*D*D*F*S*T*a*a*d*r*s;
a2 = a2  + 12.*C*D*F*F*T*T*a*a*d*r*s;
a2 = a2  + 12.*C*C*C*F*R*R*a*a*d*r*t;
a2 = a2  + 24.*C*C*D*F*R*S*a*a*d*r*t;
a2 = a2  + 24.*C*C*F*F*R*T*a*a*d*r*t;
a2 = a2  + 12.*C*D*D*F*S*S*a*a*d*r*t;
a2 = a2  + 24.*C*D*F*F*S*T*a*a*d*r*t;
a2 = a2  + 12.*C*F*F*F*T*T*a*a*d*r*t;
a2 = a2  + 12.*C*C*C*H*R*R*a*a*d*r;
a2 = a2  + 24.*C*C*D*H*R*S*a*a*d*r;
a2 = a2  + 24.*C*C*F*H*R*T*a*a*d*r;
a2 = a2  + 12.*C*D*D*H*S*S*a*a*d*r;
a2 = a2  + 24.*C*D*F*H*S*T*a*a*d*r;
a2 = a2  + 12.*C*F*F*H*T*T*a*a*d*r;
a2 = a2  + 6.*C*C*D*D*R*R*a*a*d*s*s;
a2 = a2  + 12.*C*D*D*D*R*S*a*a*d*s*s;
a2 = a2  + 12.*C*D*D*F*R*T*a*a*d*s*s;
a2 = a2  + 6.*D*D*D*D*S*S*a*a*d*s*s;
a2 = a2  + 12.*D*D*D*F*S*T*a*a*d*s*s;
a2 = a2  + 6.*D*D*F*F*T*T*a*a*d*s*s;
a2 = a2  + 12.*C*C*D*F*R*R*a*a*d*s*t;
a2 = a2  + 24.*C*D*D*F*R*S*a*a*d*s*t;
a2 = a2  + 24.*C*D*F*F*R*T*a*a*d*s*t;
a2 = a2  + 12.*D*D*D*F*S*S*a*a*d*s*t;
a2 = a2  + 24.*D*D*F*F*S*T*a*a*d*s*t;
a2 = a2  + 12.*D*F*F*F*T*T*a*a*d*s*t;
a2 = a2  + 12.*C*C*D*H*R*R*a*a*d*s;
a2 = a2  + 24.*C*D*D*H*R*S*a*a*d*s;
a2 = a2  + 24.*C*D*F*H*R*T*a*a*d*s;
a2 = a2  + 12.*D*D*D*H*S*S*a*a*d*s;
a2 = a2  + 24.*D*D*F*H*S*T*a*a*d*s;
a2 = a2  + 12.*D*F*F*H*T*T*a*a*d*s;
a2 = a2  + 6.*C*C*F*F*R*R*a*a*d*t*t;
a2 = a2  + 12.*C*D*F*F*R*S*a*a*d*t*t;
a2 = a2  + 12.*C*F*F*F*R*T*a*a*d*t*t;
a2 = a2  + 6.*D*D*F*F*S*S*a*a*d*t*t;
a2 = a2  + 12.*D*F*F*F*S*T*a*a*d*t*t;
a2 = a2  + 6.*F*F*F*F*T*T*a*a*d*t*t;
a2 = a2  + 12.*C*C*F*H*R*R*a*a*d*t;
a2 = a2  + 24.*C*D*F*H*R*S*a*a*d*t;
a2 = a2  + 24.*C*F*F*H*R*T*a*a*d*t;
a2 = a2  + 12.*D*D*F*H*S*S*a*a*d*t;
a2 = a2  + 24.*D*F*F*H*S*T*a*a*d*t;
a2 = a2  + 12.*F*F*F*H*T*T*a*a*d*t;
a2 = a2  + 6.*C*C*H*H*R*R*a*a*d;
a2 = a2  + 12.*C*D*H*H*R*S*a*a*d;
a2 = a2  + 12.*C*F*H*H*R*T*a*a*d;
a2 = a2  + 6.*D*D*H*H*S*S*a*a*d;
a2 = a2  + 12.*D*F*H*H*S*T*a*a*d;
a2 = a2  + 6.*F*F*H*H*T*T*a*a*d;
a2 = a2  + 6.*C*C*C*R*R*a*b*d*r;
a2 = a2  + 12.*C*C*D*R*S*a*b*d*r;
a2 = a2  + 12.*C*C*F*R*T*a*b*d*r;
a2 = a2  + 6.*C*D*D*S*S*a*b*d*r;
a2 = a2  + 12.*C*D*F*S*T*a*b*d*r;
a2 = a2  + 6.*C*F*F*T*T*a*b*d*r;
a2 = a2  + 6.*C*C*D*R*R*a*b*d*s;
a2 = a2  + 12.*C*D*D*R*S*a*b*d*s;
a2 = a2  + 12.*C*D*F*R*T*a*b*d*s;
a2 = a2  + 6.*D*D*D*S*S*a*b*d*s;
a2 = a2  + 12.*D*D*F*S*T*a*b*d*s;
a2 = a2  + 6.*D*F*F*T*T*a*b*d*s;
a2 = a2  + 6.*C*C*F*R*R*a*b*d*t;
a2 = a2  + 12.*C*D*F*R*S*a*b*d*t;
a2 = a2  + 12.*C*F*F*R*T*a*b*d*t;
a2 = a2  + 6.*D*D*F*S*S*a*b*d*t;
a2 = a2  + 12.*D*F*F*S*T*a*b*d*t;
a2 = a2  + 6.*F*F*F*T*T*a*b*d*t;
a2 = a2  + 6.*C*C*H*R*R*a*b*d;
a2 = a2  + 12.*C*D*H*R*S*a*b*d;
a2 = a2  + 12.*C*F*H*R*T*a*b*d;
a2 = a2  + 6.*D*D*H*S*S*a*b*d;
a2 = a2  + 12.*D*F*H*S*T*a*b*d;
a2 = a2  + 6.*F*F*H*T*T*a*b*d;
a2 = a2 -2.*C*C*K*R*R*a*d*f;
a2 = a2 -4.*C*D*K*R*S*a*d*f;
a2 = a2 -4.*C*F*K*R*T*a*d*f;
a2 = a2 -2.*D*D*K*S*S*a*d*f;
a2 = a2 -4.*D*F*K*S*T*a*d*f;
a2 = a2 -2.*F*F*K*T*T*a*d*f;
a2 = a2  + 1.*R*R*d;
a2 = a2 -2.*K*R*S*d;
a2 = a2 -4.*C*D*R*R*a*d*s;
a2 = a2 -4.*C*F*R*R*a*d*t;
a2 = a2 -4.*C*H*R*R*a*d;
a2 = a2 -4.*D*D*R*S*a*d*s;
a2 = a2 -4.*D*F*R*T*a*d*s;
a2 = a2 -4.*D*F*R*S*a*d*t;
a2 = a2 -4.*D*H*R*S*a*d;
a2 = a2 -4.*F*F*R*T*a*d*t;
a2 = a2 -4.*F*H*R*T*a*d;
a2 = a2 -2.*C*R*R*b*d;
a2 = a2 -2.*D*R*S*b*d;
a2 = a2 -2.*F*R*T*b*d;
a2 = a2  + 1.*K*K*S*S*d;
a2 = a2  + 4.*C*C*K*R*S*a*d*r;
a2 = a2  + 4.*C*D*K*S*S*a*d*r;
a2 = a2  + 4.*C*F*K*S*T*a*d*r;
a2 = a2  + 4.*C*F*K*R*S*a*d*t;
a2 = a2  + 4.*C*H*K*R*S*a*d;
a2 = a2  + 4.*D*F*K*S*S*a*d*t;
a2 = a2  + 4.*D*H*K*S*S*a*d;
a2 = a2  + 4.*F*F*K*S*T*a*d*t;
a2 = a2  + 4.*F*H*K*S*T*a*d;
a2 = a2  + 2.*C*K*R*S*b*d;
a2 = a2  + 2.*D*K*S*S*b*d;
a2 = a2  + 2.*F*K*S*T*b*d;
a2 = a2  + 1.*C*C*R*R*b*b*d;
a2 = a2  + 2.*C*D*R*S*b*b*d;
a2 = a2  + 2.*C*F*R*T*b*b*d;
a2 = a2  + 1.*D*D*S*S*b*b*d;
a2 = a2  + 2.*D*F*S*T*b*b*d;
a2 = a2  + 1.*F*F*T*T*b*b*d;
a2 = a2  + 2.*C*G*R*R*a*e;
a2 = a2  + 2.*D*G*R*S*a*e;
a2 = a2  + 2.*F*G*R*T*a*e;
a2 = a2  + 2.*C*G*K*R*S*a*e;
a2 = a2  + 2.*D*G*K*S*S*a*e;
a2 = a2  + 2.*F*G*K*S*T*a*e;
a2 = a2 -6.*C*C*C*G*R*R*a*a*e*r;
a2 = a2 -12.*C*C*D*G*R*S*a*a*e*r;
a2 = a2 -12.*C*C*F*G*R*T*a*a*e*r;
a2 = a2 -6.*C*D*D*G*S*S*a*a*e*r;
a2 = a2 -12.*C*D*F*G*S*T*a*a*e*r;
a2 = a2 -6.*C*C*D*G*R*R*a*a*e*s;
a2 = a2 -12.*C*D*D*G*R*S*a*a*e*s;
a2 = a2 -12.*C*D*F*G*R*T*a*a*e*s;
a2 = a2 -6.*C*F*F*G*T*T*a*a*e*r;
a2 = a2 -6.*C*C*F*G*R*R*a*a*e*t;
a2 = a2 -12.*C*D*F*G*R*S*a*a*e*t;
a2 = a2 -12.*C*F*F*G*R*T*a*a*e*t;
a2 = a2 -6.*C*C*G*H*R*R*a*a*e;
a2 = a2 -12.*C*D*G*H*R*S*a*a*e;
a2 = a2 -12.*C*F*G*H*R*T*a*a*e;
a2 = a2 -6.*D*D*D*G*S*S*a*a*e*s;
a2 = a2 -12.*D*D*F*G*S*T*a*a*e*s;
a2 = a2 -6.*D*F*F*G*T*T*a*a*e*s;
a2 = a2 -6.*D*D*F*G*S*S*a*a*e*t;
a2 = a2 -12.*D*F*F*G*S*T*a*a*e*t;
a2 = a2 -6.*D*D*G*H*S*S*a*a*e;
a2 = a2 -12.*D*F*G*H*S*T*a*a*e;
a2 = a2 -6.*F*F*F*G*T*T*a*a*e*t;
a2 = a2 -6.*F*F*G*H*T*T*a*a*e;
a2 = a2 -3.*C*C*G*R*R*a*b*e;
a2 = a2 -6.*C*D*G*R*S*a*b*e;
a2 = a2 -6.*C*F*G*R*T*a*b*e;
a2 = a2 -3.*D*D*G*S*S*a*b*e;
a2 = a2 -6.*D*F*G*S*T*a*b*e;
a2 = a2 -3.*F*F*G*T*T*a*b*e;
a2 = a2  + 1.*C*C*K*R*R*a*e*e;
a2 = a2  + 2.*C*D*K*R*S*a*e*e;
a2 = a2  + 2.*C*F*K*R*T*a*e*e;
a2 = a2  + 1.*D*D*K*S*S*a*e*e;
a2 = a2  + 2.*D*F*K*S*T*a*e*e;
a2 = a2  + 1.*F*F*K*T*T*a*e*e;
a2 = a2  + 4.*C*C*G*G*R*R*a*a*f;
a2 = a2  + 8.*C*D*G*G*R*S*a*a*f;
a2 = a2  + 8.*C*F*G*G*R*T*a*a*f;
a2 = a2  + 4.*D*D*G*G*S*S*a*a*f;
a2 = a2  + 8.*D*F*G*G*S*T*a*a*f;
a2 = a2  + 4.*F*F*G*G*T*T*a*a*f;
a2 = a2 -4.*C*C*G*G*R*R*a*a*s;
a2 = a2 -16.*C*D*G*G*R*S*a*a*s;
a2 = a2 -8.*C*F*G*G*R*T*a*a*s;
a2 = a2 -12.*D*D*G*G*S*S*a*a*s;
a2 = a2 -16.*D*F*G*G*S*T*a*a*s;
a2 = a2 -4.*F*F*G*G*T*T*a*a*s;
a2 = a2 -8.*C*C*G*G*R*S*a*a*r;
a2 = a2 -8.*C*D*G*G*S*S*a*a*r;
a2 = a2 -8.*C*F*G*G*S*T*a*a*r;
a2 = a2 -8.*C*F*G*G*R*S*a*a*t;
a2 = a2 -8.*D*F*G*G*S*S*a*a*t;
a2 = a2 -8.*F*F*G*G*S*T*a*a*t;
a2 = a2 -8.*C*G*G*H*R*S*a*a;
a2 = a2 -8.*D*G*G*H*S*S*a*a;
a2 = a2 -8.*F*G*G*H*S*T*a*a;
a2 = a2 -4.*C*G*G*R*S*a*b;
a2 = a2 -4.*D*G*G*S*S*a*b;
a2 = a2 -4.*F*G*G*S*T*a*b;
a2 = a2 ;
float a3 = -2.*C*C*R*R*R*a*d;
a3 = a3 -4.*C*D*R*R*S*a*d;
a3 = a3 -4.*C*F*R*R*T*a*d;
a3 = a3 -2.*D*D*R*S*S*a*d;
a3 = a3 -4.*D*F*R*S*T*a*d;
a3 = a3 -2.*F*F*R*T*T*a*d;
a3 = a3  + 2.*C*C*K*R*R*S*a*d;
a3 = a3  + 4.*C*D*K*R*S*S*a*d;
a3 = a3  + 4.*C*F*K*R*S*T*a*d;
a3 = a3  + 2.*D*D*K*S*S*S*a*d;
a3 = a3  + 4.*D*F*K*S*S*T*a*d;
a3 = a3  + 2.*F*F*K*S*T*T*a*d;
a3 = a3  + 4.*C*C*C*C*R*R*R*a*a*d*r;
a3 = a3  + 12.*C*C*C*D*R*R*S*a*a*d*r;
a3 = a3  + 12.*C*C*C*F*R*R*T*a*a*d*r;
a3 = a3  + 12.*C*C*D*D*R*S*S*a*a*d*r;
a3 = a3  + 24.*C*C*D*F*R*S*T*a*a*d*r;
a3 = a3  + 12.*C*C*F*F*R*T*T*a*a*d*r;
a3 = a3  + 4.*C*D*D*D*S*S*S*a*a*d*r;
a3 = a3  + 12.*C*D*D*F*S*S*T*a*a*d*r;
a3 = a3  + 12.*C*D*F*F*S*T*T*a*a*d*r;
a3 = a3  + 4.*C*C*C*D*R*R*R*a*a*d*s;
a3 = a3  + 12.*C*C*D*D*R*R*S*a*a*d*s;
a3 = a3  + 12.*C*C*D*F*R*R*T*a*a*d*s;
a3 = a3  + 12.*C*D*D*D*R*S*S*a*a*d*s;
a3 = a3  + 24.*C*D*D*F*R*S*T*a*a*d*s;
a3 = a3  + 12.*C*D*F*F*R*T*T*a*a*d*s;
a3 = a3  + 4.*C*F*F*F*T*T*T*a*a*d*r;
a3 = a3  + 4.*C*C*C*F*R*R*R*a*a*d*t;
a3 = a3  + 12.*C*C*D*F*R*R*S*a*a*d*t;
a3 = a3  + 12.*C*C*F*F*R*R*T*a*a*d*t;
a3 = a3  + 12.*C*D*D*F*R*S*S*a*a*d*t;
a3 = a3  + 24.*C*D*F*F*R*S*T*a*a*d*t;
a3 = a3  + 12.*C*F*F*F*R*T*T*a*a*d*t;
a3 = a3  + 4.*C*C*C*H*R*R*R*a*a*d;
a3 = a3  + 12.*C*C*D*H*R*R*S*a*a*d;
a3 = a3  + 12.*C*C*F*H*R*R*T*a*a*d;
a3 = a3  + 12.*C*D*D*H*R*S*S*a*a*d;
a3 = a3  + 24.*C*D*F*H*R*S*T*a*a*d;
a3 = a3  + 12.*C*F*F*H*R*T*T*a*a*d;
a3 = a3  + 4.*D*D*D*D*S*S*S*a*a*d*s;
a3 = a3  + 12.*D*D*D*F*S*S*T*a*a*d*s;
a3 = a3  + 12.*D*D*F*F*S*T*T*a*a*d*s;
a3 = a3  + 4.*D*F*F*F*T*T*T*a*a*d*s;
a3 = a3  + 4.*D*D*D*F*S*S*S*a*a*d*t;
a3 = a3  + 12.*D*D*F*F*S*S*T*a*a*d*t;
a3 = a3  + 12.*D*F*F*F*S*T*T*a*a*d*t;
a3 = a3  + 4.*D*D*D*H*S*S*S*a*a*d;
a3 = a3  + 12.*D*D*F*H*S*S*T*a*a*d;
a3 = a3  + 12.*D*F*F*H*S*T*T*a*a*d;
a3 = a3  + 4.*F*F*F*F*T*T*T*a*a*d*t;
a3 = a3  + 4.*F*F*F*H*T*T*T*a*a*d;
a3 = a3  + 2.*C*C*C*R*R*R*a*b*d;
a3 = a3  + 6.*C*C*D*R*R*S*a*b*d;
a3 = a3  + 6.*C*C*F*R*R*T*a*b*d;
a3 = a3  + 6.*C*D*D*R*S*S*a*b*d;
a3 = a3  + 12.*C*D*F*R*S*T*a*b*d;
a3 = a3  + 6.*C*F*F*R*T*T*a*b*d;
a3 = a3  + 2.*D*D*D*S*S*S*a*b*d;
a3 = a3  + 6.*D*D*F*S*S*T*a*b*d;
a3 = a3  + 6.*D*F*F*S*T*T*a*b*d;
a3 = a3  + 2.*F*F*F*T*T*T*a*b*d;
a3 = a3 -2.*C*C*C*G*R*R*R*a*a*e;
a3 = a3 -6.*C*C*D*G*R*R*S*a*a*e;
a3 = a3 -6.*C*C*F*G*R*R*T*a*a*e;
a3 = a3 -6.*C*D*D*G*R*S*S*a*a*e;
a3 = a3 -12.*C*D*F*G*R*S*T*a*a*e;
a3 = a3 -6.*C*F*F*G*R*T*T*a*a*e;
a3 = a3 -2.*D*D*D*G*S*S*S*a*a*e;
a3 = a3 -6.*D*D*F*G*S*S*T*a*a*e;
a3 = a3 -6.*D*F*F*G*S*T*T*a*a*e;
a3 = a3 -2.*F*F*F*G*T*T*T*a*a*e;
a3 = a3 -4.*C*C*G*G*R*R*S*a*a;
a3 = a3 -8.*C*D*G*G*R*S*S*a*a;
a3 = a3 -8.*C*F*G*G*R*S*T*a*a;
a3 = a3 -4.*D*D*G*G*S*S*S*a*a;
a3 = a3 -8.*D*F*G*G*S*S*T*a*a;
a3 = a3 -4.*F*F*G*G*S*T*T*a*a;
a3 = a3 ;
float a4 = 1.*C*C*C*C*R*R*R*R*a*a*d;
a4 = a4  + 4.*C*C*C*D*R*R*R*S*a*a*d;
a4 = a4  + 4.*C*C*C*F*R*R*R*T*a*a*d;
a4 = a4  + 6.*C*C*D*D*R*R*S*S*a*a*d;
a4 = a4  + 12.*C*C*D*F*R*R*S*T*a*a*d;
a4 = a4  + 6.*C*C*F*F*R*R*T*T*a*a*d;
a4 = a4  + 4.*C*D*D*D*R*S*S*S*a*a*d;
a4 = a4  + 12.*C*D*D*F*R*S*S*T*a*a*d;
a4 = a4  + 12.*C*D*F*F*R*S*T*T*a*a*d;
a4 = a4  + 4.*C*F*F*F*R*T*T*T*a*a*d;
a4 = a4  + 1.*D*D*D*D*S*S*S*S*a*a*d;
a4 = a4  + 4.*D*D*D*F*S*S*S*T*a*a*d;
a4 = a4  + 6.*D*D*F*F*S*S*T*T*a*a*d;
a4 = a4  + 4.*D*F*F*F*S*T*T*T*a*a*d;
a4 = a4  + 1.*F*F*F*F*T*T*T*T*a*a*d;
a4 = a4 ;
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
        if (!(p.x >= xstart && p.x <= xstart + 1.0 && p.y >= ystart && p.y <= ystart+1.0))    
            continue;
        if (roots[i] < dist)    
        {
            dist = roots[i];
            pos = p;
        }

    }
    if (dist < dist_infin)
    {
        //nor = vec3(0.+2.*d*pos.x-2.*K*d*pos.y-6.*C*C*a*d*pos.x*pos.x-8.*C*D*a*d*pos.x*pos.y-8.*C*F*a*d*pos.x*pos.z-8.*C*H*a*d*pos.x-2.*D*D*a*d*pos.y*pos.y-4.*D*F*a*d*pos.y*pos.z-4.*D*H*a*d*pos.y-2.*F*F*a*d*pos.z*pos.z-4.*F*H*a*d*pos.z-2.*H*H*a*d-4.*C*b*d*pos.x-2.*D*b*d*pos.y-2.*F*b*d*pos.z-2.*H*b*d+2.*K*d*f+4.*C*C*K*a*d*pos.x*pos.y+4.*C*D*K*a*d*pos.y*pos.y+4.*C*F*K*a*d*pos.y*pos.z+4.*C*H*K*a*d*pos.y+2.*C*K*b*d*pos.y+4.*C*C*C*C*a*a*d*pos.x*pos.x*pos.x+12.*C*C*C*D*a*a*d*pos.x*pos.x*pos.y+12.*C*C*C*F*a*a*d*pos.x*pos.x*pos.z+12.*C*C*C*H*a*a*d*pos.x*pos.x+12.*C*C*D*D*a*a*d*pos.x*pos.y*pos.y+24.*C*C*D*F*a*a*d*pos.x*pos.y*pos.z+24.*C*C*D*H*a*a*d*pos.x*pos.y+12.*C*C*F*F*a*a*d*pos.x*pos.z*pos.z+24.*C*C*F*H*a*a*d*pos.x*pos.z+12.*C*C*H*H*a*a*d*pos.x+6.*C*C*C*a*b*d*pos.x*pos.x+12.*C*C*D*a*b*d*pos.x*pos.y+12.*C*C*F*a*b*d*pos.x*pos.z+12.*C*C*H*a*b*d*pos.x-4.*C*C*K*a*d*f*pos.x+4.*C*D*D*D*a*a*d*pos.y*pos.y*pos.y+12.*C*D*D*F*a*a*d*pos.y*pos.y*pos.z+12.*C*D*D*H*a*a*d*pos.y*pos.y+12.*C*D*F*F*a*a*d*pos.y*pos.z*pos.z+24.*C*D*F*H*a*a*d*pos.y*pos.z+12.*C*D*H*H*a*a*d*pos.y+6.*C*D*D*a*b*d*pos.y*pos.y+12.*C*D*F*a*b*d*pos.y*pos.z+12.*C*D*H*a*b*d*pos.y-4.*C*D*K*a*d*f*pos.y+4.*C*F*F*F*a*a*d*pos.z*pos.z*pos.z+12.*C*F*F*H*a*a*d*pos.z*pos.z+12.*C*F*H*H*a*a*d*pos.z+6.*C*F*F*a*b*d*pos.z*pos.z+12.*C*F*H*a*b*d*pos.z-4.*C*F*K*a*d*f*pos.z+4.*C*H*H*H*a*a*d+6.*C*H*H*a*b*d-4.*C*H*K*a*d*f+2.*C*C*b*b*d*pos.x+2.*C*D*b*b*d*pos.y+2.*C*F*b*b*d*pos.z+2.*C*H*b*b*d-2.*C*K*b*d*f+4.*C*G*a*e*pos.x+2.*D*G*a*e*pos.y+2.*F*G*a*e*pos.z+2.*G*H*a*e+1.*G*b*e-1.*K*e*e+2.*C*G*K*a*e*pos.y-6.*C*C*C*G*a*a*e*pos.x*pos.x-12.*C*C*D*G*a*a*e*pos.x*pos.y-12.*C*C*F*G*a*a*e*pos.x*pos.z-12.*C*C*G*H*a*a*e*pos.x-6.*C*C*G*a*b*e*pos.x+2.*C*C*K*a*e*e*pos.x-6.*C*D*D*G*a*a*e*pos.y*pos.y-12.*C*D*F*G*a*a*e*pos.y*pos.z-12.*C*D*G*H*a*a*e*pos.y-6.*C*D*G*a*b*e*pos.y+2.*C*D*K*a*e*e*pos.y-6.*C*F*F*G*a*a*e*pos.z*pos.z-12.*C*F*G*H*a*a*e*pos.z-6.*C*F*G*a*b*e*pos.z+2.*C*F*K*a*e*e*pos.z-6.*C*G*H*H*a*a*e-6.*C*G*H*a*b*e+2.*C*H*K*a*e*e-1.*C*G*b*b*e+1.*C*K*b*e*e-2.*C*G*K*a*e*f+8.*C*C*G*G*a*a*f*pos.x+8.*C*D*G*G*a*a*f*pos.y+8.*C*F*G*G*a*a*f*pos.z+8.*C*G*G*H*a*a*f+4.*C*G*G*a*b*f-8.*C*C*G*G*a*a*pos.x*pos.y-8.*C*D*G*G*a*a*pos.y*pos.y-8.*C*F*G*G*a*a*pos.y*pos.z-8.*C*G*G*H*a*a*pos.y-4.*C*G*G*a*b*pos.y, 0.-2.*K*d*pos.x-4.*C*D*a*d*pos.x*pos.x-4.*D*D*a*d*pos.x*pos.y-4.*D*F*a*d*pos.x*pos.z-4.*D*H*a*d*pos.x-2.*D*b*d*pos.x+2.*K*K*d*pos.y+2.*C*C*K*a*d*pos.x*pos.x+8.*C*D*K*a*d*pos.x*pos.y+4.*C*F*K*a*d*pos.x*pos.z+4.*C*H*K*a*d*pos.x+6.*D*D*K*a*d*pos.y*pos.y+8.*D*F*K*a*d*pos.y*pos.z+8.*D*H*K*a*d*pos.y+2.*F*F*K*a*d*pos.z*pos.z+4.*F*H*K*a*d*pos.z+2.*H*H*K*a*d+2.*C*K*b*d*pos.x+4.*D*K*b*d*pos.y+2.*F*K*b*d*pos.z+2.*H*K*b*d-2.*K*K*d*f+4.*C*C*C*D*a*a*d*pos.x*pos.x*pos.x+12.*C*C*D*D*a*a*d*pos.x*pos.x*pos.y+12.*C*C*D*F*a*a*d*pos.x*pos.x*pos.z+12.*C*C*D*H*a*a*d*pos.x*pos.x+6.*C*C*D*a*b*d*pos.x*pos.x+12.*C*D*D*D*a*a*d*pos.x*pos.y*pos.y+24.*C*D*D*F*a*a*d*pos.x*pos.y*pos.z+24.*C*D*D*H*a*a*d*pos.x*pos.y+12.*C*D*F*F*a*a*d*pos.x*pos.z*pos.z+24.*C*D*F*H*a*a*d*pos.x*pos.z+12.*C*D*H*H*a*a*d*pos.x+12.*C*D*D*a*b*d*pos.x*pos.y+12.*C*D*F*a*b*d*pos.x*pos.z+12.*C*D*H*a*b*d*pos.x-4.*C*D*K*a*d*f*pos.x+4.*D*D*D*D*a*a*d*pos.y*pos.y*pos.y+12.*D*D*D*F*a*a*d*pos.y*pos.y*pos.z+12.*D*D*D*H*a*a*d*pos.y*pos.y+12.*D*D*F*F*a*a*d*pos.y*pos.z*pos.z+24.*D*D*F*H*a*a*d*pos.y*pos.z+12.*D*D*H*H*a*a*d*pos.y+6.*D*D*D*a*b*d*pos.y*pos.y+12.*D*D*F*a*b*d*pos.y*pos.z+12.*D*D*H*a*b*d*pos.y-4.*D*D*K*a*d*f*pos.y+4.*D*F*F*F*a*a*d*pos.z*pos.z*pos.z+12.*D*F*F*H*a*a*d*pos.z*pos.z+12.*D*F*H*H*a*a*d*pos.z+6.*D*F*F*a*b*d*pos.z*pos.z+12.*D*F*H*a*b*d*pos.z-4.*D*F*K*a*d*f*pos.z+4.*D*H*H*H*a*a*d+6.*D*H*H*a*b*d-4.*D*H*K*a*d*f+2.*C*D*b*b*d*pos.x+2.*D*D*b*b*d*pos.y+2.*D*F*b*b*d*pos.z+2.*D*H*b*b*d-2.*D*K*b*d*f+2.*D*G*a*e*pos.x+2.*C*G*K*a*e*pos.x+4.*D*G*K*a*e*pos.y+2.*F*G*K*a*e*pos.z+2.*G*H*K*a*e+1.*G*K*b*e-6.*C*C*D*G*a*a*e*pos.x*pos.x-12.*C*D*D*G*a*a*e*pos.x*pos.y-12.*C*D*F*G*a*a*e*pos.x*pos.z-12.*C*D*G*H*a*a*e*pos.x-6.*C*D*G*a*b*e*pos.x+2.*C*D*K*a*e*e*pos.x-6.*D*D*D*G*a*a*e*pos.y*pos.y-12.*D*D*F*G*a*a*e*pos.y*pos.z-12.*D*D*G*H*a*a*e*pos.y-6.*D*D*G*a*b*e*pos.y+2.*D*D*K*a*e*e*pos.y-6.*D*F*F*G*a*a*e*pos.z*pos.z-12.*D*F*G*H*a*a*e*pos.z-6.*D*F*G*a*b*e*pos.z+2.*D*F*K*a*e*e*pos.z-6.*D*G*H*H*a*a*e-6.*D*G*H*a*b*e+2.*D*H*K*a*e*e-1.*D*G*b*b*e+1.*D*K*b*e*e-2.*D*G*K*a*e*f+8.*C*D*G*G*a*a*f*pos.x+8.*D*D*G*G*a*a*f*pos.y+8.*D*F*G*G*a*a*f*pos.z+8.*D*G*G*H*a*a*f+4.*D*G*G*a*b*f-4.*C*C*G*G*a*a*pos.x*pos.x-16.*C*D*G*G*a*a*pos.x*pos.y-8.*C*F*G*G*a*a*pos.x*pos.z-8.*C*G*G*H*a*a*pos.x-4.*C*G*G*a*b*pos.x-12.*D*D*G*G*a*a*pos.y*pos.y-16.*D*F*G*G*a*a*pos.y*pos.z-16.*D*G*G*H*a*a*pos.y-8.*D*G*G*a*b*pos.y-4.*F*F*G*G*a*a*pos.z*pos.z-8.*F*G*G*H*a*a*pos.z-4.*F*G*G*a*b*pos.z-4.*G*G*H*H*a*a-4.*G*G*H*a*b-1.*G*G*b*b, 0.-4.*C*F*a*d*pos.x*pos.x-4.*D*F*a*d*pos.x*pos.y-4.*F*F*a*d*pos.x*pos.z-4.*F*H*a*d*pos.x-2.*F*b*d*pos.x+4.*C*F*K*a*d*pos.x*pos.y+4.*D*F*K*a*d*pos.y*pos.y+4.*F*F*K*a*d*pos.y*pos.z+4.*F*H*K*a*d*pos.y+2.*F*K*b*d*pos.y+4.*C*C*C*F*a*a*d*pos.x*pos.x*pos.x+12.*C*C*D*F*a*a*d*pos.x*pos.x*pos.y+12.*C*C*F*F*a*a*d*pos.x*pos.x*pos.z+12.*C*C*F*H*a*a*d*pos.x*pos.x+6.*C*C*F*a*b*d*pos.x*pos.x+12.*C*D*D*F*a*a*d*pos.x*pos.y*pos.y+24.*C*D*F*F*a*a*d*pos.x*pos.y*pos.z+24.*C*D*F*H*a*a*d*pos.x*pos.y+12.*C*D*F*a*b*d*pos.x*pos.y+12.*C*F*F*F*a*a*d*pos.x*pos.z*pos.z+24.*C*F*F*H*a*a*d*pos.x*pos.z+12.*C*F*H*H*a*a*d*pos.x+12.*C*F*F*a*b*d*pos.x*pos.z+12.*C*F*H*a*b*d*pos.x-4.*C*F*K*a*d*f*pos.x+4.*D*D*D*F*a*a*d*pos.y*pos.y*pos.y+12.*D*D*F*F*a*a*d*pos.y*pos.y*pos.z+12.*D*D*F*H*a*a*d*pos.y*pos.y+6.*D*D*F*a*b*d*pos.y*pos.y+12.*D*F*F*F*a*a*d*pos.y*pos.z*pos.z+24.*D*F*F*H*a*a*d*pos.y*pos.z+12.*D*F*H*H*a*a*d*pos.y+12.*D*F*F*a*b*d*pos.y*pos.z+12.*D*F*H*a*b*d*pos.y-4.*D*F*K*a*d*f*pos.y+4.*F*F*F*F*a*a*d*pos.z*pos.z*pos.z+12.*F*F*F*H*a*a*d*pos.z*pos.z+12.*F*F*H*H*a*a*d*pos.z+6.*F*F*F*a*b*d*pos.z*pos.z+12.*F*F*H*a*b*d*pos.z-4.*F*F*K*a*d*f*pos.z+4.*F*H*H*H*a*a*d+6.*F*H*H*a*b*d-4.*F*H*K*a*d*f+2.*C*F*b*b*d*pos.x+2.*D*F*b*b*d*pos.y+2.*F*F*b*b*d*pos.z+2.*F*H*b*b*d-2.*F*K*b*d*f+2.*F*G*a*e*pos.x+2.*F*G*K*a*e*pos.y-6.*C*C*F*G*a*a*e*pos.x*pos.x-12.*C*D*F*G*a*a*e*pos.x*pos.y-12.*C*F*F*G*a*a*e*pos.x*pos.z-12.*C*F*G*H*a*a*e*pos.x-6.*C*F*G*a*b*e*pos.x+2.*C*F*K*a*e*e*pos.x-6.*D*D*F*G*a*a*e*pos.y*pos.y-12.*D*F*F*G*a*a*e*pos.y*pos.z-12.*D*F*G*H*a*a*e*pos.y-6.*D*F*G*a*b*e*pos.y+2.*D*F*K*a*e*e*pos.y-6.*F*F*F*G*a*a*e*pos.z*pos.z-12.*F*F*G*H*a*a*e*pos.z-6.*F*F*G*a*b*e*pos.z+2.*F*F*K*a*e*e*pos.z-6.*F*G*H*H*a*a*e-6.*F*G*H*a*b*e+2.*F*H*K*a*e*e-1.*F*G*b*b*e+1.*F*K*b*e*e-2.*F*G*K*a*e*f+8.*C*F*G*G*a*a*f*pos.x+8.*D*F*G*G*a*a*f*pos.y+8.*F*F*G*G*a*a*f*pos.z+8.*F*G*G*H*a*a*f+4.*F*G*G*a*b*f-8.*C*F*G*G*a*a*pos.x*pos.y-8.*D*F*G*G*a*a*pos.y*pos.y-8.*F*F*G*G*a*a*pos.y*pos.z-8.*F*G*G*H*a*a*pos.y-4.*F*G*G*a*b*pos.y);
        nor = vec3(1.0);
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
        t = 0.;
    }
    vec3 ro = vec3(0.0, 0.0, 5.5); // camera
    ro = rotateY(-m.x*TAU)*rotateX(-m.y*PI)*ro; //camera rotation
    
    
    const float fl = 1.5; // focal length
    float dist = dist_infin;
    mat3 rota  = rotateZ(t)*rotateY(-t);
    mat3 rota_1  = rotateY(t)*rotateZ(-t);
    mat3 sky = rotateZ(0.0)*rotateX(PI/2.0);
    
   
    vec3 tot = vec3(0.0);
    
    #define AA 1
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
        float xstart = 2.0;
        float ystart = 2.0;
        BEZIER2 xb = tileCurve(xstart, 213.56);
        BEZIER2 yb = tileCurve(ystart, 752.343);
        HIT giper = giper3D(rota*ro, rota*rd, xb, yb, 0.0, 0.0);
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