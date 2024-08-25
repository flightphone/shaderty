//https://www.johndcook.com/blog/2017/02/16/simulating-seashells/

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

mat3 rotateX(float f) {
    return mat3(vec3(1.0, 0.0, 0.0), vec3(0.0, cos(f), -sin(f)), vec3(.0, sin(f), cos(f)));
}

mat3 rotateZ(float f) {
    return mat3(vec3(cos(f), -sin(f), 0.0), vec3(sin(f), cos(f), 0.0), vec3(0.0, 0.0, 1.0));

}

mat3 rotateY(float f) {
    return mat3(vec3(cos(f), 0.0, sin(f)), vec3(0.0, 1.0, 0.0), vec3(-sin(f), 0.0, cos(f)));
}

struct HIT {
    float dist;
    vec3 nor;
    vec3 pos;
};

#define TAU 6.28318530718
float aafi(vec2 p) {
    float fi = atan(p.y, p.x);
    fi += step(p.y, 0.0)*TAU;
    return fi;
}

//converts a vector on a sphere to longitude and latitude
vec2 lonlat(vec3 p) {
    float lon = aafi(p.xy) / TAU;
    float lat = aafi(vec2(p.z, length(p.xy))) / PI;
    return vec2(1.0 - lon, lat);
}

const float dist_infin = 20.0;
const HIT hit_inf = HIT(dist_infin, vec3(0.0), vec3(0.0));
#define nn 128
const float eps = 0.01;
const float kk = 0.9;
#define nu 8.0

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



vec3 culccolor(vec3 col_in, vec3 backcol, vec3 rd, vec3 light1, vec3 light2, vec3 nor) {
    vec3 col = col_in;
    float d = dot(rd, nor);
    float pst = 1.0 - smoothstep(-0.05, -0.01, d);
    col = mix(col, backcol, pst);
    // if(d < -0.05)
    //     col = backcol;

    nor *= -sign(d);
    float difu1 = dot(nor, light1);
    float difu2 = dot(nor, light2);
    float difu = max(difu1, difu2);
    col *= clamp(difu, 0.5, 1.0);
    return col;
}


vec4 near_dist(vec3 point, vec3 shell, out vec4 coo) {
    float d = dist_infin;
    vec3 nor = vec3(0.0);
    float fi = aafi(point.xy);
    float r = shell.x;
    float r2 = shell.y;
    float shift = (shell.z * TAU * nu + 2.0 * r2)/2.0;
    
    for(float i = 0.; i < nu; i++) {
        float t = fi + TAU * i;
        float z = shell.z * t - shift;
        vec2 v1 = normalize(point.xy) * r;
        vec3 nr = point - vec3(v1.x, v1.y, z);
        float d2 = length(nr) - r2;
        d2 = abs(d2);
        if(d2 < d) {
            d = d2;
            nor = nr;
            coo = vec4(r, r2, z, fi);
        }
    }

    //if (point.x > 0.0 && point.y > 0.0)
    {
        float t = TAU * nu;
        float z = shell.z * t - shift;
        vec3 op = vec3(r, 0.0, z);
        vec3 v1 = point - op;
        vec2 v2 = normalize(v1.xz) * r2;
        vec3 v3 = vec3(v2.x, 0.0, v2.y) + op;
        vec3 nr = point - v3;
        float d2 = length(nr);
        if(d2 < d) {
            d = d2;
            nor = nr;
        }
    }

    {
        float z = -shift;
        vec3 op = vec3(r, 0.0, z);
        vec3 v1 = point - op;
        vec2 v2 = normalize(v1.xz) * r2;
        vec3 v3 = vec3(v2.x, 0.0, v2.y) + op;
        vec3 nr = point - v3;
        float d2 = length(nr);
        if(d2 < d) {
            d = d2;
            nor = nr;
        }
    }

    return vec4(nor, d*kk);

}


HIT giper3D(vec3 ro, vec3 rd, vec3 shell, out vec4 coo) {
    float t = 0.;
    float f = 0.;
    float sn = 0.;
    
    for(int i = 0; i < nn; i++) {
        vec3 pos = ro + rd * t;
        vec4 d = near_dist(pos, shell, coo);
        if(d[3] < eps) {
            vec3 nor = normalize(d.xyz);
            return HIT(d[3], nor, pos);
        }
        

        t+=d[3];
        if(t >= dist_infin)
            return hit_inf;
    }

    return hit_inf;
}

vec3 GetRayDir(vec2 uv, vec3 p, vec3 l, float z) {
    vec3 f = normalize(l - p), r = normalize(cross(vec3(0, 1, 0), f)), u = cross(f, r), c = f * z, i = c + uv.x * r + uv.y * u;
    return normalize(i);
}

void mainImage(out vec4 fragColor, in vec2 fragCoord) {
    //surface (x+y+z-a)(xy+yz+zx) - kxyz = 0
    vec3 light = normalize(vec3(0.0, 0.0, -1.0)); //light
    vec3 light2 = normalize(vec3(0.0, 0.0, 1.0)); //light

    float t =  iTime/2.0;
    vec2 m = vec2(0.0, 0.0);
    //if  (iMouse.z > 0.0)
    // {
        //m = (-iResolution.xy + 2.0 * (iMouse.xy)) / iResolution.y;
        //t = 0.;
    // }
    vec3 ro = vec3(0.0, 0.0, 8.5); // camera
    ro = rotateY(-m.x * TAU) * rotateX(-m.y * PI) * ro; //camera rotation

    const float fl = 1.5; // focal length
    float dist = dist_infin;
    

    // mat3 rota = rotateZ(t)*rotateX(PI/4.0)*rotateY(PI); 
    // mat3 rota_1 = rotateY(-PI)*rotateX(-PI/4.0)*rotateZ(-t);
    mat3 rota = rotateX(t)*rotateY(t); 
    mat3 rota_1 = rotateY(-t)*rotateX(-t);
    
    //mat3 sky = rotateZ(0.0)*rotateX(PI/2.0);

    vec3 shell = vec3(2.0, 0.3, 0.67/TAU);
    //vec3 shift = getShift(shell);
    //vec3 shift = vec3(0.0);
    vec3 col = vec3(0.7, 0.7, 0.9); // background        

    vec3 tot = vec3(0.0);

    #define AA 2
    //antiblick
    for(int m = 0; m < AA; m++) for(int n = 0; n < AA; n++) {
            vec3 col = vec3(0.7, 0.7, 0.9); // background        
            // pixel coordinates
            vec2 o = vec2(float(m), float(n)) / float(AA) - 0.5;
            vec2 p = (-iResolution.xy + 2.0 * (fragCoord + o)) / iResolution.y;
            //vec2 p = (2.0*fragCoord-iResolution.xy)/iResolution.y;
            //vec3 rd = normalize( vec3(p,fl) ); // ray direction
            vec3 rd = GetRayDir(p, ro, vec3(0, 0., 0), fl); //ray direction
            vec3 backcol = vec3(1.0, 0.5, 0.5);
            vec4 coo = vec4(0.0);
            HIT giper = giper3D(rota * (ro), rota * rd, shell, coo);
            if(giper.dist < dist) {
                float x = coo[3] /TAU;
                x = fract(x);
                float y = aafi(vec2(coo[0] - length(giper.pos.xy), giper.pos.z-coo[2])) /TAU;
                col = texture(iChannel0, vec2(x,y)).rgb;
                backcol = texture(iChannel1, vec2(x,y)).rgb;
                //col = vec3(0.5, 0.5, 1.0);
                vec3 nor = rota_1 * giper.nor;
                col = culccolor(col, backcol, -rd, light, light2, nor);
            }
            tot += col;
        }
    //antiblick
    tot /= float(AA * AA);
    fragColor = vec4(tot, 1.0);

}

/////=====================================================================================
void main() {
    vec4 fragColor = vec4(0);
    mainImage(fragColor, gl_FragCoord.xy);
    gl_FragColor = fragColor;
}