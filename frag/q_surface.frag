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
const float  eps = 0.05;
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



//=============================2D pic=========================
float sdBuild(vec2 p, float f[30], float r[30], int n)
{
    float fi = aafi(p);
    for (int i = 1; i < 30; i++)
    {
        if (i == n)
            break;
        if (fi >= f[i-1] && fi <= f[i])
        {
            float f1 = fi - f[i-1];
            float f2 = f[i] - fi;
            float h1 = r[i-1] * sin(f1);
            float h2 = r[i] * sin(f2);
            vec2 v1 = vec2(r[i-1] * cos(f[i-1]), r[i-1] * sin(f[i-1]));
            vec2 v2 = vec2(r[i] * cos(f[i]), r[i] * sin(f[i]));
            vec2 v3 = v2 - v1;
            v3 *= h1 / (h1 + h2);
            vec2 v4 = v1 + v3;
            return length(p) - length(v4);
        }
    }
}

vec3 grid_plane(vec2 uv, float n)
{
        float rombr[30]; 
        rombr[0] = 1.0; 
        rombr[1] = 1.0; 
        rombr[2] = 1.0; 
        rombr[3] = 1.0; 
        rombr[4] = 1.0; 

        float rombf[30]; 
        rombf[0] = 0.0; 
        rombf[1] = PI/2.0;
        rombf[2] = PI;
        rombf[3] = 3.0*PI/2.0;
        rombf[4] = 2.0*PI;

    
        vec2 p =  vec2(fract(uv.x*20.0), fract(uv.y*10.0));
        p = (p-0.5)*2.0;
        
        vec3 col = vec3(1.0);
        float pst = sdBuild(p, rombf, rombr, 5);
        col = mix(col, vec3(1.0, 1.0, 0.0), step(pst, 0.0));
        for (int i = 0; i < 5; i++)
        {
            rombr[i] = 0.8;
        }
        pst = sdBuild(p, rombf, rombr, 5);
        col = mix(col, vec3(0.0, 0.0, 1.0), step(pst, 0.0));
        return col;
}
//=============================2D pic=========================


HIT grid(vec3 ro, vec3 rd)
{
    float dist = dist_infin;
    vec3 pos = vec3(0.0);
    vec3 nor = vec3(0.0);

    for (float j = 0.; j < nn; j++)
    {
        vec3 po = ro + rd*j*0.5;

        
        float x0 = floor(po.x);
        float y0 = floor(po.y);
        float z0 = floor(po.z);
        

        float a = ro.x-x0 - 0.5;
        float b = rd.x;
        float c = ro.y-y0 - 0.5;
        float d = rd.y;
        float e = ro.z-z0 - 0.5;
        float f = rd.z;

        
        float t = 0.25;
        float a0 = 1.*a*a + 1.*c*c + 1.*e*e-1.*t*t;
        float a1 = 2.*a*b + 2.*c*d + 2.*e*f;
        float a2 = 1.*b*b + 1.*d*d + 1.*f*f;

         vec2 roots = vec2(dist_infin);
         int nroots = quadratic(a2, a1, a0, roots);  // quartic(a4, a3, a2, a1, a0, roots);
         for (int i = 0; i < 2; i++)
            {
                if (i >= nroots)
                    break;
                if (roots[i] < 0.0)
                    continue;
                vec3 p = vec3(a, c, e) + roots[i]*rd;
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
       t = 0.;
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
    float w = 4.;
    float h = 4.;
    vec3 tot = vec3(0.0);
    
    //antiblick
    for( int m=0; m<AA; m++ )
    for( int n=0; n<AA; n++ )
    {
        vec3 col = vec3(0.7, 0.7, 0.9); // background    
            // pixel coordinates
        vec2 o = vec2(float(m),float(n)) / float(AA) - 0.5;
        vec2 p = (-iResolution.xy + 2.0*(fragCoord+o))/iResolution.y;
        vec3 rd = GetRayDir(p, ro, vec3(0,0.,0), fl); //ray direction
        HIT gr = grid(rota*ro, rota*rd);
        if (gr.dist < dist)
        {
            
            vec2 pt = lonlat(rotsp*gr.pos); //get longitude and latitude
            col = grid_plane(pt, 10.0);
            //col = texture(iChannel0, pt).rgb;
            vec3 backcol = vec3(1.0, 0.2, 0.2);
            vec3 nor = rota_1*gr.nor;
            col = ccolor(col, backcol, -rd, light, light2, nor);
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