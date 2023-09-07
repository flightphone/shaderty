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
    col *= clamp(difu, 0.1, 1.0);
    // gamma
    //col = pow( col, vec3(0.4545) );
    return col;   
}
HIT hyperI(vec3 ro, vec3 rd, float ax, float by, float cz, float h, float t)
{
    if (t > 0.0)
        return hit_inf;
    vec3 pos = ro + rd*t;
    if (abs(pos.z) > h)
        return hit_inf;
    float dist = length(ro - pos);
    vec3 nor = vec3(pos.x/ax/ax, pos.y/by/by, -pos.z/cz/cz);
    nor = normalize(nor);
    return HIT(dist, nor, pos);
}

HIT hyper3D(vec3 ro1, vec3 rd1, float ax, float by, float cz, float h)
{
    
    
    
    vec3 ro = vec3(ro1.x/ax, ro1.y/by, ro1.z/cz);
    vec3 rd = vec3(rd1.x/ax, rd1.y/by, rd1.z/cz);
    float a = rd.x*rd.x + rd.y*rd.y - rd.z*rd.z;
    float b = 2.0 * (ro.x*rd.x + ro.y*rd.y - ro.z*rd.z);
    float c = ro.x*ro.x + ro.y*ro.y - ro.z*ro.z - 1.0;
    float d = b*b - 4.0*a*c;
    if (d < 0.0)
        return hit_inf;

    
    d = pow(d, 0.5);
    float t1 = (-b + d)/2.0/a;
    float t2 = (-b - d)/2.0/a;

    HIT r = hit_inf;
    HIT r1 = hyperI(ro1, rd1,  ax,  by,  cz,  h,  t1);
    if (r1.dist < r.dist)
        r = r1;
    
    HIT r2 = hyperI(ro1, rd1,  ax,  by,  cz,  h,  t2);
    if (r2.dist < r.dist)
        r = r2;

    return r;
}

HIT sphereI(vec3 ro, vec3 rd, float t)
{
    if (t > 0.0)
        return hit_inf;
    vec3 pos = ro + rd*t;
    float dist = length(ro - pos);
    vec3 nor = vec3(pos.x, pos.y, pos.z);
    nor = normalize(nor);
    return HIT(dist, nor, pos);
}

HIT sphere3D(vec3 ro, vec3 rd,  float ra)
{
    float a = rd.x*rd.x + rd.y*rd.y + rd.z*rd.z;
    float b = 2.0 * (ro.x*rd.x + ro.y*rd.y + ro.z*rd.z);
    float c = ro.x*ro.x + ro.y*ro.y + ro.z*ro.z - ra*ra;
    float d = b*b - 4.0*a*c;
    if (d < 0.0)
        return hit_inf;

    
    d = pow(d, 0.5);
    float t1 = (-b + d)/2.0/a;
    float t2 = (-b - d)/2.0/a;

    HIT r = hit_inf;
    HIT r1 = sphereI(ro, rd, t1);
    if (r1.dist < r.dist)
        r = r1;
    
    HIT r2 = sphereI(ro, rd, t2);
    if (r2.dist < r.dist)
        r = r2;

    return r;
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    #define AA 1
    //antiblick
    vec3 tot = vec3(0.0);
    for( int m=0; m<AA; m++ )
    for( int n=0; n<AA; n++ )
    {
        // pixel coordinates
        vec2 o = vec2(float(m),float(n)) / float(AA) - 0.5;
        vec2 p = (-iResolution.xy + 2.0*(fragCoord+o))/iResolution.y;    
    //vec2 p = (2.0*fragCoord-iResolution.xy)/iResolution.y;
    vec3 light = normalize(vec3(2.0, 1.0, -1.0)); //light

	vec3 ro = vec3(0.0, 0.0, 15.0); // camera	
    const float fl = 2.5; // focal length

    vec3 rd = normalize( vec3(p,fl) ); // ray direction
    float t = iTime;
	
    vec3 col = vec3(0.7, 0.7, 0.9); // background
    //t = 0.0;
    
    float dist = dist_infin;
    vec3 nor = vec3(0.0);
    float ra = 4.0;
    mat3 sky = rotateZ(0.0)*rotateX(PI/2.0);
    float h = 0.5;

    #define ni 6.0
    #define nj 6.0
    for (float i = 0.0; i < ni; i++)
    for (float j = 0.0; j < nj; j++)
    {
        mat3 rota  = rotateX(PI/2.0 + 2.0*PI*i/ni)*rotateY(2.0*PI*j/nj)*rotateY(t)*rotateX(t);
        mat3 rota_1  = rotateX(-t)*rotateY(-t)*rotateY(-2.0*PI*j/nj)*rotateX(-PI/2.0 - 2.0*PI*i/ni);
        vec3 shift = rota_1*vec3(0.0, 0.0, ra + 0.7*h);

        HIT hyper = hyper3D(rota*(ro + shift), rota*rd, 0.3, 0.3, 0.2, h);
        if (hyper.dist < dist)
        {
            nor = rota_1*hyper.nor;
            dist = hyper.dist;
        }
    }
    HIT sphere = sphere3D(ro, rd, ra);
    if (sphere.dist < dist)
    {
        nor = sphere.nor;
        dist = sphere.dist;
    }


    vec3 backcolor = vec3(1.0, 0.0, 0.0);
    vec3 surfcolor = vec3(0.4, 0.4, 1.0);
    if (dist < dist_infin)
        //col = calcSkyReflect(rd, nor, sky);
        col = culccolor(surfcolor, backcolor, rd, light, nor);

    
    //antiblick
    col = sqrt( col );
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