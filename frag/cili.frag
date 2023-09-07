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


/////=====================================================================================

#define PI 3.14159265359
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

vec3 culccolor(vec3 col_in, vec3 backcol, vec3 rd, vec3 light, vec3 nor)
{
    vec3 col = col_in;
    float d = dot(rd, nor);
    if (d < 0.0)
        col = backcol;
    
    nor *= -sign(d);
    float difu = dot(nor, light);
    col *= clamp(difu, 0.05, 1.0);
    return col;   
}
HIT giperI(vec3 ro, vec3 rd, float ax, float by, float h, float t)
{
    if (t > 0.0)
        return hit_inf;
    vec3 pos = ro + rd*t;
    if (abs(pos.z) > h)
        return hit_inf;
    float dist = length(ro - pos);
    vec3 nor = vec3(pos.x/ax/ax, pos.y/by/by, 0.0);
    nor = normalize(nor);
    return HIT(dist, nor, pos);
}

HIT giper3D(vec3 ro1, vec3 rd1, float ax, float by,  float h, float ra)
{
    
    
    
    vec3 ro = vec3(ro1.x/ax, ro1.y/by, ro1.z);
    vec3 rd = vec3(rd1.x/ax, rd1.y/by, rd1.z);
    float a = rd.x*rd.x + rd.y*rd.y;
    float b = 2.0 * (ro.x*rd.x + ro.y*rd.y);
    float c = ro.x*ro.x + ro.y*ro.y - ra*ra;
    float d = b*b - 4.0*a*c;
    if (d < 0.0)
        return hit_inf;

    
    d = pow(d, 0.5);
    float t1 = (-b + d)/2.0/a;
    float t2 = (-b - d)/2.0/a;

    HIT r = hit_inf;
    HIT r1 = giperI(ro1, rd1,  ax,  by,  h,  t1);
    if (r1.dist < r.dist)
        r = r1;
    
    HIT r2 = giperI(ro1, rd1,  ax,  by,  h,  t2);
    if (r2.dist < r.dist)
        r = r2;

    return r;
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    
    
    vec3 light = normalize(vec3(2.0, 1.0, -1.0)); //light
	vec3 ro = vec3(0.0, 0.0, 25.0); // camera	
    const float fl = 2.5; // focal length
    float t = iTime;

    
    
    float dist = dist_infin;
    mat3 rota  = rotateZ(PI/3.0+t)*rotateX(PI/2.0*3.0/4.0+t)*rotateY(-t/3.0);
    mat3 rota_1  = rotateY(t/3.0)*rotateX(-PI/2.0*3.0/4.0-t)*rotateZ(-PI/3.0-t);
    mat3 sky = rotateZ(0.0)*rotateX(PI/2.0);

    vec2 p = (2.0*fragCoord-iResolution.xy)/iResolution.y;
    vec3 rd = normalize( vec3(p,fl) ); // ray direction
    vec3 col = vec3(0.7, 0.7, 0.9); // background


    HIT giper = giper3D(rota*(ro), rota*rd, 1.0, 1.0, 6.0, 3.0);
    if (giper.dist < dist)
    {
        
        //col = backcol;
        vec3 nor = rota_1*giper.nor;
        vec3 pos = giper.pos;
        float v = 2.0;
        float f = (pos.z + 3.0)*v;
        float d = dot(normalize(vec2(pos.x, pos.y)), vec2(cos(f), sin(f)));
        float pst = smoothstep(0.7, 1.0, d);
        
        if (pst > 0.0) 
        //col = mix(col, vec3(1.0, 0.0, 0.0), pst);
        {
            col = vec3(1.0, 1.0, 0.0);
            vec3 backcol = vec3(0.0, 1.0, 1.0);
            col = culccolor(col, backcol, rd, light, nor);
            // gamma
            col = pow( col, vec3(0.4545) ); 
        }
        //reflect
        //col = calcSkyReflect(rd, nor, sky);
    }
    
   
    fragColor = vec4(col,1.0);
    
}

/////=====================================================================================
void main()
{
    vec4 fragColor = vec4(0);
    mainImage(fragColor,  gl_FragCoord.xy);
    gl_FragColor = fragColor;
}