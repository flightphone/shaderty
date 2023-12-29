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

#define PI  3.14159265359
#define TAU 6.28318530718
#define rot(f) mat2(cos(f), -sin(f), sin(f), cos(f))


const float dist_infin = 10.0;
#define nn 128
const float eps = 0.001;

vec3 sdfColor;
vec3 resColor;
vec3 col1 = vec3(.2);
vec3 col2 = vec3(0.4745098039215686, 0.26666666666666666, 0.23137254901960785);
vec3 col3 = vec3(1., 0.8431, 0.);
float resReflect = 0.5;
float sdfReflect;

float mtime = 0.;
//https://mathcurve.com/courbes2d/larme/larme.shtml
float larme2(vec2 p, float a)
{
    float k = 0.3;
    p.x += a*(k - 0.95)/2.0;
    float zoom = 1.2;
    float x = clamp(p.x, -a*0.95, a*k);
    float f = acos(x/a);
    float y = a*sin(f)*pow(sin(f/2.), 2.)*zoom;
    float d = length(p - vec2(x, y));
    if (abs(p.y) < abs(y))
    {
        sdfReflect = 0.;
        sdfColor = col2;
    }
    else
    {
        sdfReflect = 0.05;
        sdfColor = col1;
    }
    
    f = acos(k);
    y = a*sin(f)*pow(sin(f/2.), 2.)*zoom;
    float y2 = clamp(p.y, -y, y);
    float d2 = length(p - vec2(k*a, y2));
    if (d2 < d)
    {
        d = d2;
        if (p.x > x)
        {
            sdfReflect = 0.05;
            sdfColor = col1;
        }
        else
        {
            sdfReflect = 0.;
            sdfColor = col2;
        }
    }
    return d;
}
float larme(vec3 p, float a)
{
    float l = length(p.xy);
    float d = larme2(vec2(p.z, l), a);
    return d*0.3 - 0.02;
}

float ep = 0.3;
float men2(vec2 p)
{
    p*=vec2(8., 11.);
    vec2 a = floor(p);
    float res = 0.;
    /*
    if (a.y == 0. || a.y == 10.)
        res = 1.0;
    */
    res += smoothstep(0., ep, p.y)*smoothstep(1., 1.-ep, p.y);
    res += smoothstep(10., 10.+ep, p.y)*smoothstep(11., 11.-ep, p.y);
    
    a.y = a.y - 2.;
    float mx, mi;
    
    if (mod(a.y, 2.) == 0.)
    {
        float h = abs(3.-a.y);
        float d = (a.y <=2.)? 0.:1.;
        mi = 5.-h-d;
        mx = 5.+h-2.*d;
        //res = smoothstep(mi, mi+ep, p.x)*smoothstep(mx+1., mx+1.-ep, p.x);
        if (a.x <= mx && a.x >= mi) 
            res = 1.;
        
    }    
    
    
    if (mod(a.x, 2.) == 0.)
    {
        mx = 0.;
        mi = 1.;
        if (a.x <= 4.)
        {
            float h = 3.-a.x;
            mi = clamp(2.-h, 0., 6.);
            mx = 3.+ h;
            

        }

        if (a.x == 6.0)
        {
            mi = 3.;
            mx = 6.;
            
        }
        
        if (a.y <= mx && a.y >= mi) 
            res = 1.;
        
        
        
    }
    
    return res;
}



//https://www.shadertoy.com/view/4ljBWd
float greekwave(vec2 U)
{
    U = fract(U)-.5;                           // local tile coords in [-.5,.5]²
    
    float v, d = length(U),                    // distance to tile center
          t = 4.5*PI,                          // note the time delay with position
          a = t * max(0.,.5-d );               // angle ~time, and decrease with d
    U *= rot(a);   // rotate frame by angle(t,d)
    v = U.y;   
    float O = smoothstep(-1.,1.,v*500. );  
    return O;     
}


float men3(vec2 p)
{
    p*=vec2(6.,7.);
    vec2 a = floor(p);
    float res = 1.0;
    if (a.y == 0. || a.y == 6.)
        res = 0.;
    if (a.x == 0. && a.y < 6. && a.y > 2.)    
        res = 0.;
    if (a.y == 2. && (a.x < 2. || a.x == 5.))    
        res = 0.;
    if (a.y == 4. && a.x >= 2. && a.x <=4.)    
        res = 0.;
    if (a.x == 3. && a.y < 4.)    
        res = 0.;
    return res;    
}


float map(vec3 p) {
    p.yz *= rot(-PI/2.);
    p.yz *= rot(mtime/2.);
    p.xz *= rot(mtime);
    
    
    
    
    
    float a = 3.;
    float d = larme(p, a);

    if (d <= 10.0*eps && sdfColor == col1)
    {
        float n = 20., dlon = TAU/n, lon = mod(atan(p.y, p.x), TAU), dx = fract(lon/dlon);
        float k = 0.3;
        p.z += a*(k - 0.95)/2.0;
        float h = clamp(p.z, -a*0.95, a*k);
        if (h > -a*0.9 && h < -a*0.8)
        {
            float dy = (h + a*0.9)/0.1/a, r = men2(vec2(1.-dx, dy));
            sdfColor = mix(col1, col3, r);
            sdfReflect = mix(0.05, 0.5, r);
        }

        if (h > -a*0.8 && h < -a*0.7)
        {
            float dy = (h + a*0.8)/0.1/a, r = greekwave(vec2(1.-dx, dy));
            if (r == 1.)
            {
                sdfReflect = 0.5;
                sdfColor = col3;
            }
        }
        if (h > -a*0.7 && h < -a*0.5)
        {
            float dy = (h + a*0.7)/0.2/a, r =1.-men2(vec2(1.-dx, dy));
            sdfColor = mix(col1, col3, r);
            sdfReflect = mix(0.05, 0.5, r);
            /*
            if (r == 1.)
            {
                sdfReflect = 0.5;
                sdfColor = col3;
            }
            */
        }

        if (h > a*0.2 && h < a*0.3)
        {
            float dy = (h - a*0.2)/0.1/a, r = men3(vec2(dx, dy));
            if (r == 0.)
            {
                sdfReflect = 0.5;
                sdfColor = col3;
            }
        }

        
    }


    resColor = sdfColor;
    resReflect = sdfReflect;
    return d;
}


vec3 csky(vec3 p) {
    float n = 6., m = 6., dlat = PI / n, dlon = TAU / m;
    float lon = mod(atan(p.y, p.x), TAU), lat = atan(length(p.xy), p.z);
    float fo = fract(lon / dlon), fa = fract(lat / dlat);

    float pst = fo * fa * (1. - fo) * (1. - fa);
    pst=pow(pst, 2.);

    pst = smoothstep(0.0, pow(0.0625, 2.), pst);
    pst = clamp(pst, 0.1, 1.0);
    return vec3(pst);
}


// https://iquilezles.org/articles/normalsSDF
vec3 calcNormal(in vec3 pos) {
    const float h = 0.0001; // replace by an appropriate value
    const vec2 k = vec2(1, -1);
    return normalize(k.xyy * map(pos + k.xyy * h) +
        k.yyx * map(pos + k.yyx * h) +
        k.yxy * map(pos + k.yxy * h) +
        k.xxx * map(pos + k.xxx * h));
}

vec3 GetRayDir(vec2 uv, vec3 p, vec3 l, float z) {
    vec3 f = normalize(l - p), r = normalize(vec3(f.z, 0, -f.x)), u = cross(f, r), c = f * z, i = c + uv.x * r + uv.y * u;
    return normalize(i);
}
/*
#if HW_PERFORMANCE==0
#define AA 1
#else
#define AA 2
#endif
*/
#define AA 2
/*
vec3 ccolor(vec3 col, vec3 rd, vec3 light, vec3 nor) {
    float d = dot(rd, nor);
    nor *= -sign(d);
    float difu = dot(nor, light);
    col *= clamp(difu, 0.3, 1.0);
    return col;
}
*/
//new color 29.12.2023
vec3 ccolor(vec3 col, vec3 rd, vec3 light, vec3 nor) {
    //float d = dot(rd, nor);
    //nor *= -sign(d);

    vec3 R = reflect (light, nor);
    float shininess=10.0;
    float specular    =  pow(max(dot(R, rd), 0.), shininess);
    
    float difu = dot(nor, light);
    //col = col*(col*clamp(difu, 0., 1.0) + 0.3);// + vec3(.5)*specular*specular;
    col = col*clamp(difu, 0.3, 1.0);// + vec3(.5)*specular*specular;
    return col;
}

vec3 calccolor(vec3 col, vec3 backcol, vec3 rd, vec3 light1, vec3 light2, vec3 nor) {
    float difu1 = dot(nor, light1);
    float difu2 = dot(nor, light2);
    float difu = max(difu1, difu2);
    //col *= clamp(difu, 0.3, 1.0);

    vec3 R1 = reflect (light1, nor);
    vec3 R2 = reflect (light2, nor);
    float shininess=10.0;
    float specular1    =  pow(max(dot(R1, rd), 0.), shininess);
    float specular2    =  0.0;//pow(max(dot(R2, rd), 0.), shininess);
    float specular = max(specular1, specular2);

    //col+= vec3(.5)*specular*specular;
    col = col*(col*clamp(difu, 0., 1.0) + 0.3) + vec3(.5)*specular*specular;
    //col = col*clamp(difu, 0.3, 1.0) + vec3(.5)*specular*specular;
    

    return col;
}

// https://knarkowicz.wordpress.com/2016/01/06/aces-filmic-tone-mapping-curve/
vec3 ACESFilm(vec3 x){
    return clamp((x * (2.51 * x + 0.03)) / (x * (2.43 * x + 0.59) + 0.14), 0.0, 1.0);
}




void mainImage(out vec4 fragColor, in vec2 fragCoord) {
    vec3 light = normalize(vec3(1.0, -1.0, -1.)); //light
    vec3 light2 = normalize(vec3(0., 1.0, 0.)); //light
    vec2 mo = vec2(0.0, 0.0);
    mtime = iTime;
    //if  (iMouse.z > 0.0)
    {
        mo = (-iResolution.xy + 2.0 * (iMouse.xy)) / iResolution.y;
        mtime = 0.;
    }
    vec3 ro = vec3(0.0, 0.0, 6.0); // camera
    //camera rotation
    ro.yz *= rot(mo.y * PI);
    ro.xz *= rot(-mo.x * TAU);

    const float fl = 1.5; // focal length
    float dist = dist_infin;

    vec3 b1 = vec3(0.23529411764705882, 0.4235294117647059, 0.7725490196078432), b2 = vec3(0.3686274509803922, 0.5725490196078431, 0.8941176470588236);
    vec3 bg = mix(b2, b1*b1, fragCoord.y / iResolution.y);   
    //vec3 bg = vec3(0.);
    //antialiasing
    vec3 tot = vec3(0.0);
    for(int m = 0; m < AA; m++) for(int n = 0; n < AA; n++) {
            vec2 o = vec2(float(m), float(n)) / float(AA) - 0.5;
            vec2 p = (-iResolution.xy + 2.0 * (fragCoord + o)) / iResolution.y;
            vec3 rd = GetRayDir(p, ro, vec3(0, 0., 0), fl); //ray direction
            vec3 col = bg; // background  
            //==========================raymatch=============================
            float td = 0.;
            vec3 pos = vec3(0.);
            for(int i = 0; i < nn; i++) {
                pos = ro + rd * td;
                float h = map(pos);
                if(h < eps || td >= dist_infin)
                    break;
                td += h;
            }
            if(td < dist_infin) {
                col = resColor;
                vec3 nor = calcNormal(pos);

                //reflection

                vec3 psk = reflect(rd, nor);
                vec3 c2 = csky(psk);
                
                //col = ccolor(col, rd, light, nor);
                col = calccolor(col, col, rd, light, light2, nor);
                
                col = mix(col, c2, resReflect);

                //col += c2*0.1;

            }
            //==========================raymatch=============================
            tot += col;
        }
    tot = tot / float(AA);    
    //tot = pow(tot, vec3(0.4545)) / float(AA);
    //tot = ACESFilm(tot);
    //tot = pow(tot, vec3(0.7)) / float(AA);
    
    //tot = pow(tot, vec3(0.4545));
    //antialiasing
    fragColor = vec4(tot, 1.0);
}
/////=====================================================================================
void main() {
    vec4 fragColor = vec4(0);
    mainImage(fragColor, gl_FragCoord.xy);
    gl_FragColor = fragColor;
}