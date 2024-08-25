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
#define nn 50.
#define newton 5
float csurf = 0.5;
float dist_infin =2.2;


float isf(vec3 p) {
   float x = p.x, y = -p.z, z = p.y; 
   return (2. * y * (y * y - 3. * x * x) * (1. - z * z) + (x * x + y * y) * (x * x + y * y) - (9. * z * z - 1.) * (1. - z * z));// IMPLICIT SURFACE Function
}


float cassinian(vec3 p) {

        p.xy*=rot(PI/4.);
        p.yz*=rot(PI/4.);
        float res = 1., b = 1.5, r = 1.8;
        
        for (float i = 0.; i < 2.; i++)
            for (float j = 0.; j < 2.; j++)
                for (float k = 0.; k < 2.; k++)
                    res *= length(p - vec3(r * (i - 0.5), r * (j - 0.5), r * (k - 0.5)));
        return res - b * b * b * b * r * r * r * r;
    }    


float glz() {
    float t = iTime / 5.;
    float st = mod(floor(t), 4.);
    float res;
    if(st == 0.)
        res = 1.;
    if(st == 1.)
        res = cos(fract(t) * PI / 2.);//(1.- fract(t))*(1.- fract(t));
    if(st == 2.)
        res = 0.;
    if(st == 3.)
        res = sin(fract(t) * PI / 2.); //fract(t)*fract(t);   
    return res;
}


float map(vec3 p) {
    if (csurf == 0.0)
        return isf(p);
    else
    if (csurf == 1.0)
        return cassinian(p);
    else       
        return cassinian(p)*csurf + (1.-csurf)*isf(p);
    
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
    
    csurf = glz();
    dist_infin = 2.2;
    float    hh =3.5;
        
    vec3 light = normalize(vec3(0.0, 1.0, -2.5)); //light
    vec3 light2 = normalize(vec3(0.0, -1.0, 2.5)); //light
    vec2 mo = 1.5*cos(1.5*iTime + vec2(0,11));
    //if  (iMouse.z > 0.0)
    {
        //mo = (-iResolution.xy + 2.0 * (iMouse.xy)) / iResolution.y;
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