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
//Collection of implicit surfaces. implicit surfaces, raytracing, binary search
/*
Rendering implicit surfaces. Using raytracing and binary searchy. 
Here, these same surfaces are obtained by creating grids using an algorithm 
3D Marching Cubes: https://flightphone.github.io/paramgeometry.html
*/
#define PI  3.14159265359
#define TAU 6.28318530718
#define rot(f) mat2(cos(f), -sin(f), sin(f), cos(f))
#define nn 128.
#define newton 10

float dist_infin =2.2;


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


vec3 calcNormal(in vec3 p) {
    const float eps = 0.0001;
    vec2 q = vec2(0.0, eps);
	vec3 res =  vec3(map(p+q.yxx) - map(p-q.yxx), 
			    map(p+q.xyx) - map(p-q.xyx),
			    map(p+q.xxy) - map(p-q.xxy));
    return normalize(res);
}

vec3 getPoint(vec3 a, vec3 b, float v0, float v1) {
            vec3 m;
            //binary search with  n iterations, n = newton
            for (int i = 0; i < newton; i++) {
                m = (a+b)*0.5;
                float v = map(m);
                if (v == 0.)
                    break;

                if (sign(v) * sign(v0) <= 0.) {
                    v1 = v;
                    b = m;
                }
                else {
                    v0 = v;
                    a = m;
                }
            }
            return m;
        }

vec3 GetRayDir(vec2 uv, vec3 p, vec3 l, float z) {
    vec3 f = normalize(l - p), r = normalize(vec3(f.z, 0, -f.x)), u = cross(f, r), c = f * z, i = c + uv.x * r + uv.y * u;
    return normalize(i);
}

vec3 calccolor(vec3 col_in, vec3 backcol, vec3 rd, vec3 light1, vec3 light2, vec3 nor) {
    vec3 col = col_in;
    float d = dot(rd, nor);
    if(d < 0.0)
        col = backcol;

    float difu1 = dot(nor, light1);
    float difu2 = dot(nor, light2);
    float difu = max(difu1, difu2);
   

    vec3 R1 = reflect (light1, nor);
    vec3 R2 = reflect (light2, nor);
    float shininess=20.0;
    float specular1    =  pow(max(dot(R1, rd), 0.), shininess);
    float specular2    =  pow(max(dot(R2, rd), 0.), shininess);
    float specular = max(specular1, specular2);
    col = col*(col*clamp(difu, 0., 1.0) + 0.3) + vec3(.5)*specular*specular;
    // gamma
    col = pow( col, vec3(0.4545) );  
    return col;
}

/*
#if HW_PERFORMANCE==0
#define AA 1
#else
#define AA 2
#endif
*/
#define AA 1
void mainImage(out vec4 fragColor, in vec2 fragCoord) {
    
    //csurf = glz();
    dist_infin = 1.2;
    float    hh =3.5;
        
    vec3 light = normalize(vec3(0.0, 1.0, -2.5)); //light
    vec3 light2 = normalize(vec3(0.0, -1.0, 2.5)); //light
    vec2 mo = 1.5*cos(1.5*iTime + vec2(0,11));
    //if  (iMouse.z > 0.0)
    {
        mo = (-iResolution.xy + 2.0 * (iMouse.xy)) / iResolution.y;
    }
    vec3 ro = vec3(0.0, 0.0, hh ); // camera
    //camera rotation
    ro.yz *= rot(mo.y);
    ro.xz *= rot(-mo.x - 1.57);

    const float fl = 1.5; // focal length
    float dist = dist_infin;

    vec3 b1 = vec3(0.23529411764705882, 0.4235294117647059, 0.7725490196078432), b2 = vec3(0.3686274509803922, 0.5725490196078431, 0.8941176470588236);
    vec3 bg = 2.0*mix(b2, b1*b1, 0.);
    vec3 col1 = vec3(0.7304607400847158,0.5906188409113381,0.3005437944049895);
    vec3 col2 = vec3(0.7230551289161951, 0.0060488330203860696, 0.0060488330203860696);
 
    
    //antialiasing
    vec3 tot = vec3(0.0);
    for(int m = 0; m < AA; m++) 
    for(int n = 0; n < AA; n++) {
            vec2 o = vec2(float(m), float(n)) / float(AA) - 0.5;
            vec2 p = (-iResolution.xy + 2.0 * (fragCoord + o)) / iResolution.y;
            vec3 rd = GetRayDir(p, ro, vec3(0, 0., 0), fl); //ray direction
            vec3 col = bg * bg; // background  
            
            
            //STEP 1. Calculating bounding sphere
            float d = length(cross(ro, rd));
            if (d >= dist)
            {
                 tot += col;
                 continue;
            }
            /*
            STEP 2.
            ray tracing inside the bounding sphere, 
            searching for a segment with different signs of the function value 
            at the ends of the segment
            */
            float td = abs(dot(ro, rd));
            d = sqrt(dist*dist - d*d);
            vec3 pos0 = ro + rd * (td - d);
            vec3 pos1 = ro + rd * (td + d);
            vec3 rd0 = pos1 - pos0;
            vec3 pos = pos0;
            float val0 = map(pos0);
            for(float i = 1.; i < nn; i++) {
                pos1 = pos0 + rd0 * i / (nn - 1.);
                float val1 = map(pos1);
                if (sign(val0)*sign(val1) <= 0.)
                {
                    //different signs of the function value  at the ends of the segment
                    //STEP 3. binary search to clarify the intersection of a ray with a surface.
                    col = col1;
                    pos = getPoint(pos, pos1, val0, val1);
                    vec3 nor = calcNormal(pos);
                    col = calccolor(col, col2, -rd, light, light2, nor);
                    break;
                }
                val0 = val1;
                pos = pos1;
            }
            tot += col;
        }
    tot = tot / float(AA)/ float(AA);
    //antialiasing
    fragColor = vec4(tot, 1.0);
}
/////=====================================================================================
void main() {
    vec4 fragColor = vec4(0);
    mainImage(fragColor, gl_FragCoord.xy);
    gl_FragColor = fragColor;
}