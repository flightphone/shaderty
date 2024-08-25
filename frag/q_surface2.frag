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
#define nn 60.0
const float  eps = 50.05;
const float dist_infin = 1000.0;
#define AA 1


struct HIT
{
    float dist;
    vec3 nor;
    vec3 pos;
};

const HIT hit_inf = HIT(dist_infin, vec3(0.0), vec3(0.0));


vec3 ccolor(vec3 col_in, vec3 backcol, vec3 rd, vec3 light1, vec3 light2, vec3 nor)
{
    vec3 col = col_in;
    float d = dot(rd, nor);
    if (d < 0.0)
        col = backcol;
    
    nor *= -sign(d);
    float difu1 = dot(nor, light1);
    float difu2 = dot(nor, light2);
    float difu = max(difu1, difu2);
        col *= clamp(difu, 0.1, 1.0);
    return col;   
}

float aafi(vec2 p) {
    float fi = atan(p.y, p.x);
    fi += step(p.y, 0.0)*TAU;
    return fi;
}

//converts a vector on a sphere to longitude and latitude

vec2 lonlat (vec3 p)
{
    float lon = aafi(p.xy)/2.0/PI;
    float lat = aafi(vec2(p.z, length(p.xy)))/PI;
    return vec2(1.0-lon, lat);
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


float qcubic(float a, float b, float c) {
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



vec3 tileCurve(float p_x, float n)
{
    //return vec3(-3.0, 3.0, 0.0);
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

    vec2 a2 = v0 - 2.0*v1 + v2;
    vec2 a1 = (2.0*v1 - 2.0*v0);
    vec2 a0 = v0;
    return vec3(a2.y, a1.y, a0.y);
}



HIT grid(vec3 ro, vec3 rd)
{
    float dist = dist_infin;
    vec3 nor = vec3(0.0, 0.0, -1.0);
    vec3 pos = plane(ro, rd, vec3(0.), nor);
    // if  (length(pos - ro) < nn)
    //     dist = 2.0*length(pos - ro);
    

    for (float j = 0.; j < nn; j++)
    {
        vec3 po = ro + rd*j*0.3;
        vec3 u = tileCurve(po.x, 287.5668);
        vec3 v = tileCurve(po.y, 117.3456);
        float x0 = floor(po.x);
        float y0 = floor(po.y);

        float a = ro.x-x0;
        float b = rd.x;
        float c = ro.y-y0;
        float d = rd.y;
        float e = ro.z;
        float f = rd.z;

        
        float k = u.x;
        float l = u.y;
        float m = u.z;

        float K = v.x;
        float L = v.y;
        float M = v.z;

        float t = 10.0;
        
        float a0 = 1.*K*a*a*c*c*k*t + 1.*L*a*a*c*k*t + 1.*M*a*a*k*t + 1.*K*a*c*c*l*t + 1.*L*a*c*l*t + 1.*M*a*l*t + 1.*K*c*c*m*t + 1.*L*c*m*t + 1.*M*m*t-1.*e;
        float a1 = 2.*K*a*a*c*d*k*t + 1.*L*a*a*d*k*t + 2.*K*a*c*d*l*t + 1.*L*a*d*l*t + 2.*K*c*d*m*t + 1.*L*d*m*t + 2.*K*a*b*c*c*k*t + 2.*L*a*b*c*k*t + 2.*M*a*b*k*t + 1.*K*b*c*c*l*t + 1.*L*b*c*l*t + 1.*M*b*l*t-1.*f;
        float a2 = 1.*K*a*a*d*d*k*t + 1.*K*a*d*d*l*t + 1.*K*d*d*m*t + 4.*K*a*b*c*d*k*t + 2.*L*a*b*d*k*t + 2.*K*b*c*d*l*t + 1.*L*b*d*l*t + 1.*K*b*b*c*c*k*t + 1.*L*b*b*c*k*t + 1.*M*b*b*k*t;
        float a3 = 2.*K*a*b*d*d*k*t + 1.*K*b*d*d*l*t + 2.*K*b*b*c*d*k*t + 1.*L*b*b*d*k*t;
        float a4 = 1.*K*b*b*d*d*k*t;


         vec4 roots = vec4(dist_infin);
         int nroots = quartic(a4, a3, a2, a1, a0, roots);
         for (int i = 0; i < 4; i++)
            {
                if (i >= nroots)
                    break;
                if (roots[i] < 0.0)
                    continue;
                vec3 p = vec3(a, c, e) + roots[i]*rd;
                if (p.x < -0.1 || p.x > 1.1 || p.y < -0.1 || p.y > 1.1)
                    continue;
                if (roots[i] < dist)    
                {
                    dist = roots[i];
                    pos = p;
                }

            }
            if (dist < dist_infin)
            {
                nor = vec3(0.+2.*K*k*t*pos.x*pos.y*pos.y+2.*L*k*t*pos.x*pos.y+2.*M*k*t*pos.x+1.*K*l*t*pos.y*pos.y+1.*L*l*t*pos.y+1.*M*l*t, 0.+2.*K*k*t*pos.x*pos.x*pos.y+1.*L*k*t*pos.x*pos.x+2.*K*l*t*pos.x*pos.y+1.*L*l*t*pos.x+2.*K*m*t*pos.y+1.*L*m*t, 0.-1.);
                nor = normalize(nor);
                return HIT(dist, nor, pos); 
            }
    }

    return hit_inf;
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
    vec3 light = normalize(vec3(0.0, -5.0, -5.0)); //light
    //vec3 light2 = normalize(vec3(0.0, 0.0, 1.0)); //light
    vec3 light2 = -light;
    //vec2 p = (2.0*fragCoord-iResolution.xy)/iResolution.y;
    float t = iTime/4.0;
    vec2 m = vec2(0.0, 0.0);
    //if  (iMouse.z > 0.0)
    {
       m = (-iResolution.xy + 2.0*(iMouse.xy))/iResolution.y;
       //t = 0.;
    }
    vec3 ro = vec3(0.0, 0., 4.); // camera
    const float fl = 1.5; // focal length
    float dist = dist_infin;
    ro = rotateY(-m.x*TAU)*rotateX(-m.y*PI)*ro; //camera rotation
    mat3 rota  = rotateZ(t);
    mat3 rota_1  = rotateZ(-t);
    mat3 rotsp = rotateZ(-iTime)*rotateX(PI/2.0 + PI/5.0); 
    //mat3 rotsp = rotateZ(-iTime*0.6)*rotateX(PI/5.0);
    // mat3 rota  = rotateZ(t + PI/4.0);
    // mat3 rota_1  = rotateZ(-t - PI/4.0);
    //mat3 rota_1  = rotateY(t)*rotateZ(-t);
    
    vec3 tot = vec3(0.0);
    
    //antiblick
    // for( int m=0; m<AA; m++ )
    // for( int n=0; n<AA; n++ )
    // {
        vec3 col = vec3(0.7, 0.7, 0.9); // background    
            // pixel coordinates
        //vec2 o = vec2(float(m),float(n)) / float(AA) - 0.5;
        //vec2 p = (-iResolution.xy + 2.0*(fragCoord+o))/iResolution.y;
        vec2 p = (2.0*fragCoord-iResolution.xy)/iResolution.y;
        vec3 rd = GetRayDir(p, ro, vec3(0,0.,0), fl); //ray direction
        HIT gr = grid(rota*ro, rota*rd);
        if (gr.dist < dist)
        {
            
            col = vec3(0.5, 0.5, 1.0);
            vec3 backcol = col; //vec3(1.0, 0.2, 0.2);
            vec3 nor = rota_1*gr.nor;
            col = ccolor(col, backcol, -rd, light, light2, nor);
        }
        tot += col;
    // }

    // //antiblick
    // tot /= float(AA*AA);
    fragColor = vec4(tot,1.0);
    
}

/////=====================================================================================
void main()
{
    vec4 fragColor = vec4(0);
    mainImage(fragColor,  gl_FragCoord.xy);
    gl_FragColor = fragColor;
}