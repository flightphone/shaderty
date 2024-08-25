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
const float  eps = 3.0;
const float dist_infin = 1000.0;
const float ngrid = 60.0;
const float rr = 3.;
const float dt = 3.;
#define AA 2


struct HIT
{
    float dist;
    vec3 nor;
    vec3 pos;
};

const HIT hit_inf = HIT(dist_infin, vec3(0.0), vec3(dist_infin));

//https://iquilezles.org/articles/smin/
// polynomial smooth min 2 (k=0.1)
float smin( float a, float b, float k )
{
    float h = max( k-abs(a-b), 0.0 )/k;
    return min( a, b ) - h*h*k*(1.0/4.0);
}

vec2 smin2( float a, float b, float k )
{
    float h =  max( k-abs(a-b), 0.0 )/k;
    float m = h*h*0.5;
    float s = m*k*(1.0/2.0);
    return (a<b) ? vec2(a-s,m) : vec2(b-s,1.0-m);
}
//https://iquilezles.org/articles/smin/


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
        col *= clamp(difu, 0.4, 1.0);
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
vec3 grid_plane(vec2 uv, float n)
{
        /*
        vec3 col = vec3(1.0);
        float h = floor(uv.x*n) + floor(uv.y*n);
        float pst = step(1., mod(h, 2.0));
        col = mix(col, vec3(0., 0., 1.), pst);
        */
        vec3 col = texture(iChannel0, uv).rgb;
        return col;
}
//=============================2D pic=========================
HIT plane(vec3 ro, vec3 rd, vec3 po, vec3 nor)
{
    float t = dot(nor, (po - ro)) / dot(nor, rd);
    if (t < 0.0)
        return hit_inf;
    vec3 pos = ro + t*rd;
    return HIT(t, nor, pos);    
    
}

// df(x)/dx
//analitic function https://www.shadertoy.com/view/4sBGDy
vec3 nTorus( in vec3 pos, vec2 tor )
{
	return normalize( pos*(dot(pos,pos)- tor.y*tor.y - tor.x*tor.x*vec3(1.0,1.0,-1.0)));
}



HIT tor(vec3 ro, vec3 rd, float x0, float y0, float r, float R, float w, float h, float sh)
{
    vec3 pos = vec3(0.0);
    vec3 nor = vec3(0.0);

    float dist = dist_infin;
    float a = ro.x-x0;
    float b = rd.x;
    float c = ro.y-y0;
    float d = rd.y;
    float e = ro.z-r;
    float f = rd.z;
    float a22 = a*a;
    float b22 = b*b;
    float c22 = c*c;
    float d22 = d*d;
    float e22 = e*e;
    float f22 = f*f;
    float a44 = a22*a22;
    float R22 = R*R;
    float r22 = r*r;

    float a0 = a44 + 2.*a22*c22 + 2.*a22*e22-2.*R22*a22-2.*a22*r22 + c22*c22 + 2.*c22*e22-2.*R22*c22-2.*c22*r22 + e22*e22 + 2.*R22*e22-2.*e22*r22 + R22*R22-2.*R22*r22 + r22*r22;
    float a1 = 4.*a22*a*b + 4.*a22*c*d + 4.*a22*e*f + 4.*a*b*c22 + 4.*c*c22*d + 4.*c22*e*f + 4.*a*b*e22 + 4.*c*d*e22 + 4.*e*e22*f-4.*R22*a*b-4.*R22*c*d + 4.*R22*e*f-4.*a*b*r22-4.*c*d*r22-4.*e*f*r22;
    float a2 = 6.*a22*b22 + 2.*a22*d22 + 2.*a22*f22 + 2.*b22*c22 + 6.*c22*d22 + 2.*c22*f22 + 2.*b22*e22 + 2.*d22*e22 + 6.*e22*f22-2.*R22*b22-2.*R22*d22 + 2.*R22*f22-2.*b22*r22-2.*d22*r22-2.*f22*r22 + 8.*a*b*c*d + 8.*a*b*e*f + 8.*c*d*e*f;
    float a3 = 4.*a*b*b22 + 4.*a*b*d22 + 4.*a*b*f22 + 4.*b22*c*d + 4.*c*d*d22 + 4.*c*d*f22 + 4.*b22*e*f + 4.*d22*e*f + 4.*e*f*f22;
    float a4 = 1.*b22*b22 + 2.*b22*d22 + 2.*b22*f22 + d22*d22 + 2.*d22*f22 + f22*f22;
    vec4 roots = vec4(dist_infin);
    int nroots = quartic(a4, a3, a2, a1, a0, roots);  // quartic(a4, a3, a2, a1, a0, roots);
    
    for (int i = 0; i < 4; i++)
    {
        if (i >= nroots)
            break;
        if (roots[i] < 0.0)
            continue;
            
        vec3 p = vec3(a, c, e) + roots[i]*rd;
        if (abs(p.x + x0) > w || abs(p.y + y0) > h)
            continue;
        if (length(p.xy) > R)
          continue;        
        if (p.z > sh)  
          continue;
        if (roots[i] < dist)    
        {
            dist = roots[i];
            pos = p;
        }

    }
    if (dist < dist_infin)
    {
        nor = nTorus(pos, vec2(R, r));
        nor = normalize(nor);
        pos.x += x0;
        pos.y += y0;
        pos.z += r;
        return HIT(dist , nor, pos);
    }
    return hit_inf;
}

HIT sphere(vec3 ro, vec3 rd, float x0, float y0, float t, float dlt, float w, float h)
{
    vec3 pos = vec3(0.0);
    vec3 nor = vec3(0.0);

    float dist = dist_infin;
    float a = ro.x-x0;
    float b = rd.x;
    float c = ro.y-y0;
    float d = rd.y;
    float e = ro.z + dlt;
    float f = rd.z;

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
        if (p.z - dlt < 0.) 
            continue;
        if (abs(p.x + x0) > w || abs(p.y + y0) > h)
            continue;
        
            

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
        pos.x += x0;
        pos.y += y0;
        pos.z -= dlt;
        return HIT(dist , nor, pos);
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
    
    //vec2 p = (2.0*fragCoord-iResolution.xy)/iResolution.y;
    float t = iTime;
    vec2 m = vec2(0.0, 0.0);
    //if  (iMouse.z > 0.0)
    {
       m = (-iResolution.xy + 2.0*(iMouse.xy))/iResolution.y;
       //t = 0.;
    }
    
    
    
    const float fl = 1.5; // focal length
    float dist = dist_infin;
    //ro = rotateY(-m.x*TAU)*rotateX(-m.y*PI)*ro; //camera rotation
    mat3 rota  = rotateZ(PI/6.0);
    mat3 rota_1  = rotateZ(-PI/6.0);
    

    vec3 ro = vec3(0.0, -10., 4.); // camera
    vec3 light = normalize(vec3(0.0, 0., 1.0)); //light
    //vec3 light2 = normalize(vec3(0.0, 0.0, 1.0)); //light
    vec3 light2 = -light;

    vec3 rdm = GetRayDir(m, ro, vec3(0,0.,0), fl); //ray direction mouse
    HIT gr = plane(rota*ro, rota*rdm, vec3(0.), vec3(0., 0.0, 1.0));
    float x0 = gr.pos.x;
    float y0 = gr.pos.y;

    float w = 20.;
    float h = 12.;
    
    vec3 tot = vec3(0.0);
    
    float dlt = dt*cos(t);
    float rseg = sqrt(rr*rr - dlt*dlt);
    float Rtor = rseg * eps;
    float rtor = (rr*rr - Rtor*Rtor - dlt*dlt)/(dlt - rr)/2.0;
    float sh = -(dlt + rtor)/(rr + rtor)*rtor;

    //r = (h^2 - R^2 - d^2)/(d-h)
    
    //antiblick
    for( int m=0; m<AA; m++ )
    for( int n=0; n<AA; n++ )
    {
        dist = dist_infin;
        vec3 col = vec3(0.7, 0.7, 0.9); // background    
        vec3 pos = vec3(0.);
        vec3 nor = vec3(0.0);
            // pixel coordinates
        vec2 o = vec2(float(m),float(n)) / float(AA) - 0.5;
        vec2 p = (-iResolution.xy + 2.0*(fragCoord+o))/iResolution.y;
        vec3 rd = GetRayDir(p, ro, vec3(0,0.,0), fl); //ray direction
       
        
        HIT gr = plane(rota*ro, rota*rd, vec3(0.), vec3(0., 0.0, 1.0));
        if (gr.dist < dist)
        {
            if (abs(gr.pos.x) <= w && abs(gr.pos.y) <= h)
            if (length(gr.pos.xy - vec2(x0, y0)) > rseg)    
            {
                dist = gr.dist;
                nor = gr.nor;
            }
        }
        
            
        gr =  sphere(rota*ro, rota*rd, x0, y0, rr, dlt, w, h);
        if (gr.dist < dist) 
        {
            nor = gr.nor;
            dist = gr.dist;
        }
        
        gr =  tor(rota*ro, rota*rd, x0, y0, rtor, Rtor, w, h, sh);
        if (gr.dist < dist) 
        {
            nor = gr.nor;
            dist = gr.dist;
        }
        
        if (dist < dist_infin)
        {
            pos = rota*ro + rota*rd*dist;
            col = grid_plane(0.5 + pos.xy/vec2(2.*w, 2.*h), ngrid);
            
            //col = mix(col, vec3(1.0, 0.0, 0.0), pst);
            vec3 backcol = vec3(1.0, 0.0, 0.);

            vec3 nor = rota_1*nor;
            col = ccolor(col, col, -rd, light, light2, nor);
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