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

const float v = 0.35;

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
    return col;   
}
HIT giperI(vec3 ro, vec3 rd, float h, float t)
{
    
    if (t > 0.0)
        return hit_inf;
    float dist = dist_infin;
    vec3 nor = vec3(0.0);
    vec3 pos = ro + rd*t;
    
    if (abs(pos.z) > h)
        return hit_inf;
    
    
    float f = (pos.z)*v;
    float d = dot(normalize(vec2(pos.x, pos.y)), vec2(cos(f), sin(f)));    
    if (d > 0.98)
        dist = length(ro - pos);

    f = (pos.z)*v + PI;
    d = dot(normalize(vec2(pos.x, pos.y)), vec2(cos(f), sin(f)));    
    if (d > 0.98)
        dist = length(ro - pos);    
    
    if (dist < dist_infin)
    {
        nor = vec3(pos.x, pos.y, 0.0);
        nor = normalize(nor);
    }
    return HIT(dist, nor, pos);
}

HIT giper3D(vec3 ro, vec3 rd, float h, float ra)
{
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
    HIT r1 = giperI(ro, rd,   h,  t1);
    if (r1.dist < r.dist)
        r = r1;
    
    HIT r2 = giperI(ro, rd,    h,  t2);
    if (r2.dist < r.dist)
        r = r2;

    return r;
}

HIT planeY(vec3 ro, vec3 rd)
{
    vec3 col = vec3(0.0);
    vec3 pos = vec3(0.0);

    float t = -(ro.y)/rd.y;
    if (t >= 0.0)
		return hit_inf;

    pos = ro + rd*t;
    float dist = length(pos - ro);
    return HIT(dist, vec3(0.0, 1.0, 0.0), pos);      
}


HIT planeW(vec3 ro, vec3 rd, float ra, float h)
{
    float f = h*v;
    vec3 ro1 = rotateZ(f)*ro;
    vec3 rd1 = rotateZ(f)*rd;
    HIT plane = planeY(ro1, rd1);
    if (plane.pos.z > h + ra*0.05 || plane.pos.z < h - ra*0.05 || abs(plane.pos.x) > ra || abs(plane.pos.x) < ra*0.05 )
        return hit_inf;
    plane.nor = rotateZ(-f)*plane.nor;    
    return plane;
    
}



void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    
    
    vec3 light = normalize(vec3(2.0, 1.0, -1.0)); //light

	vec3 ro = vec3(0.0, 0.0, 25.0); // camera	
    const float fl = 2.5; // focal length

    
    float t = iTime;
	//t = 3.0;

    
    
    float dist = dist_infin;
    mat3 rota  =   rotateZ(t)*rotateX(PI/2.0);
    mat3 rota_1  = rotateX(-PI/2.0)*rotateZ(-t);
    //mat3 sky = rotateZ(0.0)*rotateX(PI/2.0);

    #define AA 1
    //antiblick
    vec3 tot = vec3(0.0);
    for( int m=0; m<AA; m++ )
    for( int n=0; n<AA; n++ )
    {
    vec3 col = vec3(0.7, 0.7, 0.9); // background
    vec3 backcol = col;
    vec3 nor = vec3(0.0);    
    //vec2 p = (2.0*fragCoord-iResolution.xy)/iResolution.y;
    vec2 o = vec2(float(m),float(n)) / float(AA) - 0.5;
    vec2 p = (-iResolution.xy + 2.0*(fragCoord+o))/iResolution.y;
    
    vec3 rd = normalize( vec3(p,fl) ); // ray direction

    vec3 pos = vec3(0.0);
    vec3 col1 = vec3(1.0, 0.5, 0.5);
    vec3 col2 =  vec3(0.5, 0.5, 1.0);

    float ra = 4.0;
    float h = 12.0;
    #define nn 17.0
    HIT giper = giper3D(rota*(ro), rota*rd, h, ra);
    if (giper.dist < dist)
    {
        
        dist = giper.dist;
        nor = rota_1*giper.nor;
        
        col = vec3(1.0, 1.0, 0.0);
        backcol = vec3(0.0, 1.0, 1.0);
        col = culccolor(col, backcol, rd, light, nor);
        // gamma
        //col = pow( col, vec3(0.4545) ); 
        //col1 = calcSkyReflect(rd, nor, sky);
        //col = mix(col, col1, 0.5);
        //reflect
        //col = calcSkyReflect(rd, nor, sky);
    }
    float dist2 = dist;
    
    for (float i = 0.0; i <= nn; i++)
    {
        float hh = (i - nn/2.0)/nn*2.0*h;
        HIT plane = planeW(rota*(ro), rota*rd, ra, hh);
        if (plane.dist < dist2)
        {
            dist2 = plane.dist;
            pos = plane.pos;
            nor = rota_1*plane.nor;

            col1 = vec3(1.0, 0.2, 0.2);
            col2 =  vec3(0.2, 0.2, 1.0);
            
            if (mod(i, 3.0) < 1.0)    
            {
                col1 = vec3(0.5, 0.5, 0.5);
                col2 = vec3(0.2, 1.0, 0.2);
            }

            if (mod(i, 4.0) < 1.0)    
            {
                col1 = vec3(1.0, 0.2, 0.2);
                col2 = vec3(0.2, 1.0, 0.2);
            }
            
            
        }
    }

    if (dist2 < dist)
    {
        dist = dist2;
        col = col1;
        if (pos.x < 0.0)
            col = col2;
        backcol = col;
        col = culccolor(col, backcol, rd, light, nor);
        //col = pow( col, vec3(0.4545) ); 
    }

    // if (dist < dist_infin)
    // {
    //     col1 = calcSkyReflect(rd, nor, sky);
    //     col = mix(col, col1, 0.5);
    // }
        col = sqrt( col );
        tot += col;
    }
    //antiblick
    tot /= float(AA*AA);
   
    fragColor = vec4(tot,1.0);
    //fragColor = vec4(col,1.0);
    
}

/////=====================================================================================
void main()
{
    vec4 fragColor = vec4(0);
    mainImage(fragColor,  gl_FragCoord.xy);
    gl_FragColor = fragColor;
}