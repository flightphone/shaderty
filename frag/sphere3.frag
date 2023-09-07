precision highp float;
precision highp int;
precision lowp sampler2D;

// #ifdef GL_ES
// precision mediump float;
// #endif

uniform vec2 u_resolution;
uniform vec2 u_mouse;
uniform float u_time;
uniform sampler2D u_tex0;

#define iResolution u_resolution
#define iTime u_time
#define iChannel0 u_tex0
#define texture texture2D

/////=====================================================================================
//https://iquilezles.org/articles/distfunctions2d/
#define PI 3.14159265359
#define TOOP 6.28318530718
#define INV_TAU = 0.1591549430918953
#define INV_PI = 0.31830988618
//#define PI 3.142

mat3 rotateX(float f)
{
    return mat3(
    vec3(1.0,    0.0,      0.0),
    vec3(0.0,	 cos(f),  -sin(f)), 	
	vec3(.0, sin(f), cos(f))
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

mat3 rotateZ(float f)
{
    return mat3(
    vec3(cos(f),    -sin(f),  0.0),
    vec3(sin(f),	 cos(f),  0.0), 	
	vec3(0.0, 0.0, 1.0)
    );
    
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
    
    float lon = aafi(p.xy)/TOOP;
    float lat = aafi(vec2(p.z, length(p.xy)))/PI;
    lon = 1.0 - lon;
    return vec2(lon, lat);

    /*
    vec3 HPN = normalize(p);
    return vec2(
        (0.5 + atan(HPN.x, HPN.y) * INV_TAU),
        1.0 - (0.5 - asin(-HPN.z) * INV_PI));
    */

    
}
mat2 rotate2d(float _angle){
    return mat2(cos(_angle),-sin(_angle),
                sin(_angle),cos(_angle));
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    float r = 0.7;
    vec3 col = vec3(0.7, 0.7, 0.9);
    //vec3 col = vec3(1.0, 0.0, 0.0);
    vec3 light = vec3(-10.0, 0.0, 10.0);
    
    vec2 p = (2.0*fragCoord-iResolution.xy)/iResolution.y;
    float l = length(p);
    if (l > r ) 
    {
         fragColor = vec4(col, 1.0);
         return;
    }
    //get vector on sphere in Orthographic Camera
    float z = r*sin(acos(l/r));
    vec3 sp = vec3(p, z);


    vec3 sp_rot = rotateZ(-iTime*0.6)*rotateX(PI/2.0 + PI/5.0)*sp; //rotate sphere
    
    vec2 pt = lonlat(sp_rot); //get longitude and latitude
    
    col = texture(iChannel0, pt).rgb;  //get color from texture

    //fix problems
    vec3 col2 = texture(iChannel0, vec2(0.001, pt.y)).rgb;
    float pst = smoothstep(pow(pt.y, 0.5)/150.0, 0.0, pt.x); 
    col = mix(col, col2, pst);      
    pst = smoothstep(1.0-pow(pt.y, 0.5)/150.0, 1.0, pt.x);
    col = mix(col, col2, pst);      
    

    //Lambertian
    float diffuse = dot(normalize(sp), normalize(light));
    diffuse = clamp(diffuse, 0.3, 1.0);
    col *= diffuse;

    //gamma
    col = pow( col, vec3(0.4545) );

    fragColor = vec4(col, 1.0);
    
}

/////=====================================================================================
void main()
{
    vec4 fragColor = vec4(0);
    mainImage(fragColor,  gl_FragCoord.xy);
    gl_FragColor = fragColor;
}