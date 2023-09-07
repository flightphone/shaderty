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

struct RAY
{
    vec3 ro;
    vec3 rd;
};
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

HIT cube3D(vec3 ro, vec3 rd, vec3 light)
{
    
    vec3 col = vec3(1.0, 1.0, 0.0);
    vec3 backcol = vec3(0.0, 1.0, 1.0);
    vec3 rnor = vec3(0.0);
    vec3 nor = vec3(0.0);
    vec3 p = vec3(0.0);
    float dist = dist_infin;
    float h = 2.0;
    

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
    
    /*
    nor = vec3(0.0, 1.0, 0.0);
    plan = planeY(ro - nor*h, rd,  h, h);
    if (plan.dist < dist)
    {
        dist = plan.dist;
        p = plan.pos;
        rnor = nor;
    }
    */
    //col = texture(iChannel1, vec2(0.5 - plan.pos.x/2.0/w, 0.5 - plan.pos.y/2.0/h)).rgb;
    float d = dot(rd, rnor);
    if (d > 0.0)
    {
        vec2 pos = vec2(0.0);
        if (abs(rnor.z) > 0.0)
        {
            pos = p.xy;
        }
        if (abs(rnor.y) > 0.0)
        {
            pos = p.xz;
        }
        if (abs(rnor.x) > 0.0)
        {
            pos = p.zy;
        }
        col = texture(iChannel1, vec2(0.5 + pos.x/2.0/h, 0.5 - pos.y/2.0/h)).rgb;
    }
    col = culccolor(col, backcol, rd, light, rnor);
    return HIT(dist, col, p);      

}


void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    
    vec2 p = (2.0*fragCoord-iResolution.xy)/iResolution.y;

    vec3 light = normalize(vec3(0, 1.0, -1.0)); //light

	vec3 ro = vec3(0.0, 0.0, 25.0); // camera	
    const float fl = 2.5; // focal length

    vec3 rd = normalize( vec3(p,fl) ); // ray direction
    float t = iTime;
	
    vec3 col = vec3(0.7, 0.7, 0.9); // background
    float dist = dist_infin;

    //render plan
    //t = 0.0;
    mat3 rotaPlan = rotateZ(PI/3.0)*rotateX(PI/2.0*3.0/4.0);
    mat3 rotaCub  = rotateZ(PI/3.0+t)*rotateX(PI/2.0*3.0/4.0+t)*rotateY(t);
    vec3 shift = vec3(2.0, 4.0, -4.0);
    vec3 shiftPlan = vec3(.0, .0, .0);


    HIT cube = cube3D(rotaCub*(ro + shift), rotaCub*rd, rotaCub*light);
    if (cube.dist < dist)
    {
        dist = cube.dist;
        col = cube.col;
        // gamma
        col = pow( col, vec3(0.4545) ); 
    }
    float w = 6.0;
    float h = 8.0;
    HIT plan = planeW(rotaPlan*(ro + shiftPlan), rotaPlan*rd, h, w); 
    if (plan.dist < dist)
    {
        dist = plan.dist;
        col = vec3(1.0, 0.0, 0.0);
        vec3 backcol = vec3(0.0, 1.0, 0.0);
        col = texture(iChannel0, vec2(0.5 + plan.pos.x/2.0/w, 0.5 - plan.pos.y/2.0/h)).rgb;
        col = culccolor(col, backcol, rotaPlan*rd, rotaPlan*light, vec3(0.0, 0.0, 1.0));

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