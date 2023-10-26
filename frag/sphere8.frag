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
#define iChannel0 u_tex0
#define iChannel1 u_tex1

#define texture texture2D
#define iMouse u_mouse


/////=====================================================================================

#define PI 3.14159265359
#define TAU 6.28318530718
#define nn 20.
const float dist_infin = 1000.0;
#define AA 2


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
        col *= clamp(difu, 0.5, 1.0);
    return col;   
}

float aafi(vec2 p) {
    float fi = atan(p.y, p.x);
    fi += step(p.y, 0.0)*TAU;
    return fi;
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



float plane(vec3 ro, vec3 rd, vec3 po, vec3 nor)
{
    float t = dot(nor, (po - ro)) / dot(nor, rd);
    if (t < 0.)
        t = dist_infin;
    return t;
}


HIT grid(vec3 ro, vec3 rd)
{
    float dist = dist_infin;
    vec3 pos = vec3(0.0);
    vec3 nor = vec3(0.0);

    float a = ro.x;
    float b = rd.x;
    float c = ro.y;
    float d = rd.y;
    float e = ro.z;
    float f = rd.z;

    float t = 2.5;
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
        float fi = aafi(pos.xy);
        float f = fract(fi/TAU * nn);
        float fn = floor(fi/TAU * nn);
        
        if (f < 0.03)
        {
            nor = vec3(0.+2.*pos.x, 0.+2.*pos.y, 0.+2.*pos.z);
            nor = normalize(nor);
            pos = vec3(0., 0., 1.);
            // if (mod(fn, 2.0) == 0.)
            //     pos = vec3(0., 1., 0.);
            // else    
            //     pos = vec3(0., 0., 1.);

            return HIT(dist, nor, pos); 
        }
        
        else
        {
            dist = dist_infin;
            float tdist = dist;
            float p1 = fn*TAU/nn;
            float p2 = fn*TAU/nn + TAU/nn;
            vec3 po = vec3(0.);
            vec3 nr = vec3(-sin(p1), cos(p1), 0.);
            tdist = plane(ro, rd, po, nr);
            pos = ro + rd*tdist;
            if (length(pos) > t)
                tdist = dist_infin;
            if (tdist < dist && tdist > 0.0)
            {
                nor = nr;
                dist = tdist;
            }
           
            nr = vec3(-sin(p2), cos(p2), 0.);
            tdist = plane(ro, rd, po, nr);
            pos = ro + rd*tdist;
            if (length(pos) > t)
                tdist = dist_infin;
            if (tdist < dist && tdist > 0.0)
            {
                nor = nr;
                dist = tdist;
            }

            if (dist < dist_infin)
            {
                pos = vec3(1., 1., 0.);
                return HIT(dist, nor, pos);
            }
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
    vec3 light = normalize(vec3(5.0, 5.0, 0.0)); //light
    vec3 light2 = -light;
    float t = iTime/2.0;
    vec2 m = vec2(0.0, 0.0);
    //if  (iMouse.z > 0.0)
    {
       m = (-iResolution.xy + 2.0*(iMouse.xy))/iResolution.y;
       //t = 0.;
    }
    vec3 ro = vec3(0.0, -3., 5.); // camera
    const float fl = 1.5; // focal length
    
    ro = rotateY(-m.x*TAU)*rotateX(-m.y*PI)*ro; //camera rotation
    mat3 rota  = rotateZ(t);
    mat3 rota_1  = rotateZ(-t);
    vec3 tot = vec3(0.0);
    
    //antiblick
    for( int m=0; m<AA; m++ )
    for( int n=0; n<AA; n++ )
    {
        float dist = dist_infin;
        vec3 col = vec3(0.7, 0.7, 0.9); // background    
            // pixel coordinates
        vec2 o = vec2(float(m),float(n)) / float(AA) - 0.5;
        vec2 p = (-iResolution.xy + 2.0*(fragCoord+o))/iResolution.y;
        vec3 rd = GetRayDir(p, ro, vec3(0,0.,0), fl); //ray direction
        HIT gr = grid(rota*ro, rota*rd);
        if (gr.dist < dist)
        {
            
            
            col = gr.pos;//vec3(0.5, 0.5, 1.0);
            vec3 backcol = col;//vec3(1.0, 0.2, 0.2);
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