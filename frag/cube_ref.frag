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
    vec3 col;
    vec3 pos;
};

const float dist_infin = 1000.0;

vec3 culccolor(vec3 col_in, vec3 backcol, vec3 rd, vec3 light, vec3 nor)
{
    vec3 col = col_in;
    float d = dot(rd, nor);
    if (d < 0.0)
        col = backcol;
    
    //nor *= -sign(d); sign(d)*
    nor *= -sign(d);
    float difu = dot(nor, light);
    col *= clamp(difu, 0.05, 1.0);
    return col;   
}

HIT planeX(vec3 ro, vec3 rd, float h, float w)
{
    vec3 col = vec3(0.0);
    vec3 p = vec3(0.0);

    float dist = dist_infin;
    float t = -(ro.x)/rd.x;
    if (t >= 0.0)
		return HIT(dist, col, p);

    p = ro + rd*t;
    if (p.z < -w || p.y < -h || p.z > w  || p.y > h)
        return HIT(dist, col, p);  

    dist = length(p - ro);
    return HIT(dist, col, p);      
}


HIT planeY(vec3 ro, vec3 rd, float h, float w)
{
    vec3 col = vec3(0.0);
    vec3 p = vec3(0.0);

    float dist = dist_infin;
    float t = -(ro.y)/rd.y;
    if (t >= 0.0)
		return HIT(dist, col, p);

    p = ro + rd*t;
    if (p.x < -w || p.z < -h || p.x > w  || p.z > h)
        return HIT(dist, col, p);  

    dist = length(p - ro);
    return HIT(dist, col, p);      
}





HIT planeW(vec3 ro, vec3 rd, float h, float w)
{
    vec3 col = vec3(0.0);
    vec3 p = vec3(0.0);

    float dist = dist_infin;
    float t = -(ro.z)/rd.z;
    if (t >= 0.0)
		return HIT(dist, col, p);

    p = ro + rd*t;
    if (p.x < -w || p.y < -h || p.x > w  || p.y > h)
        return HIT(dist, col, p);  

    dist = length(p - ro);
    return HIT(dist, col, p);      
}

HIT cube3D(vec3 ro, vec3 rd, float h)
{
    
    vec3 rnor = vec3(0.0);
    vec3 nor = vec3(0.0);
    vec3 p = vec3(0.0);
    float dist = dist_infin;
    
    

    nor = vec3(0.0, 0.0, -1.0);
    HIT plan = planeW(ro - nor*h, rd,  h, h);
    if (plan.dist < dist)
    {
        dist = plan.dist;
        p = plan.pos;
        rnor = nor;
    }
    
    nor = vec3(0.0, 0.0, 1.0);
    plan = planeW(ro - nor*h, rd, h, h);
    if (plan.dist < dist)
    {
        dist = plan.dist;
        p = plan.pos;
        rnor = nor;
    }

    nor = vec3(-1.0, 0.0, 0.0);
    plan = planeX(ro - nor*h, rd, h, h);
    if (plan.dist < dist)
    {
        dist = plan.dist;
        p = plan.pos;
        rnor = nor;
    }

    nor = vec3(1.0, 0.0, 0.0);
    plan = planeX(ro - nor*h, rd, h, h);
    if (plan.dist < dist)
    {
        dist = plan.dist;
        p = plan.pos;
        rnor = nor;
    }

    nor = vec3(0.0, -1.0, 0.0);
    plan = planeY(ro - nor*h, rd,  h, h);
    if (plan.dist < dist)
    {
        dist = plan.dist;
        p = plan.pos;
        rnor = nor;
    }
    
    
    nor = vec3(0.0, 1.0, 0.0);
    plan = planeY(ro - nor*h, rd,  h, h);
    if (plan.dist < dist)
    {
        dist = plan.dist;
        p = plan.pos;
        rnor = nor;
    }
    
    return HIT(dist, rnor, p);      

}

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

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    
    float t = iTime*0.5;
    //t = 0.0;
    vec2 p = (2.0*fragCoord-iResolution.xy)/iResolution.y;

    vec3 light = normalize(vec3(0, 1.0, -1.0)); //light

	vec3 ro = vec3(0.0, 0.0, 25.0); // camera	
    const float fl = 2.5; // focal length

    vec3 rd = normalize( vec3(p,fl) ); // ray direction
    mat3 sky = rotateZ(t)*rotateX(PI/2.0);
    
    
	
    vec2 fon = lonlat(sky*rd); //get longitude and latitude
    vec3 col = texture(iChannel0, fon).rgb;  //get color from texture
    col = vec3(0.7, 0.7, 0.9); 
    
    float dist = dist_infin;

    
    mat3 rotaCub  = rotateZ(PI/3.0+t)*rotateX(PI/2.0*3.0/4.0+t)*rotateY(t);
    mat3 rotaCub_1  = rotateY(-t)*rotateX(-PI/2.0*3.0/4.0-t)*rotateZ(-PI/3.0-t);
    //mat3 rotaCub = rotateZ(0.0);
    //mat3 rotaCub_1 = rotateZ(0.0);
    float h = 4.0;
    HIT cube = cube3D(rotaCub*(ro), rotaCub*rd, h);
    if (cube.dist < dist)
    {
        dist = cube.dist;
        vec3 nor = cube.col;
        
        //texture cube=========================================
        // vec2 pos = vec2(0.0);
        // if (abs(nor.z) > 0.0)
        // {
        //     pos = cube.pos.xy;
        // }
        // if (abs(nor.y) > 0.0)
        // {
        //     pos = cube.pos.xz;
        // }
        // if (abs(nor.x) > 0.0)
        // {
        //     pos = cube.pos.zy;
        // }
        // col = texture(iChannel1, vec2(0.5 + pos.x/2.0/h, 0.5 - pos.y/2.0/h)).rgb;
        //texture cube=========================================
        nor = rotaCub_1*nor;
        vec3 rd_r = reflect(rd, nor);
        vec2 fon = lonlat(sky*rd_r); //get longitude and latitude
        vec3 col_r = texture(iChannel0, fon).rgb;  //get color from texture
        //col = mix(col, col_r, 0.5);
        col = col_r;
        
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