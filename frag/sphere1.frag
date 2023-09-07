#ifdef GL_ES
precision mediump float;
#endif

uniform vec2 u_resolution;
uniform vec2 u_mouse;
uniform float u_time;
uniform sampler2D u_tex0;

#define iResolution u_resolution
#define iTime u_time
#define iChannel0 u_tex0
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

mat3 nmCamera()
{
	return mat3( vec3(1.0, 0.0, 0.0),
                 vec3(0.0, 1.0, 0.0),   
                 vec3(.0, 0.0, 1.0)
    
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
    float lon = aafi(p.xy)/2.0/PI;
    float lat = aafi(vec2(p.z, length(p.xy)))/PI;
    return vec2(1.0-lon, lat);
}

//  Function from Iñigo Quiles
//  https://www.shadertoy.com/view/MsS3Wc
vec3 hsb2rgb( in vec3 c ){
    vec3 rgb = clamp(abs(mod(c.x*6.0+vec3(0.0,4.0,2.0),
                             6.0)-3.0)-1.0,
                     0.0,
                     1.0 );
    rgb = rgb*rgb*(3.0-2.0*rgb);
    return c.z * mix(vec3(1.0), rgb, c.y);
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    float t = iTime*0.3;
    float r = 0.4;
    vec3 col = vec3(0.7, 0.7, 0.9);
    vec3 light = vec3(-10.0, 0.0, 10.0);
    
    vec2 p = (2.0*fragCoord-iResolution.xy)/iResolution.y;
    //sky box
    // camera-to-world transformation
    mat3 ca = nmCamera();
    ca = rotateZ(-t - PI)*rotateX(PI/2.0)*ca;

    // focal length
    const float fl = 2.5;
    // ray direction
    vec3 rd = ca * normalize( vec3(p,fl) );
    vec2 fon = lonlat(rd); //get longitude and latitude
    col = texture(iChannel0, fon).rgb;  //get color from texture
    //fix problems
    //vec3 col2 = texture(iChannel0, vec2(0.0, fon.y)).rgb;
    //float pst = smoothstep( 0.01, 0.0, fon.x); 
    // col = mix(col, col2, pst);      
    // pst = smoothstep( 0.99, 1.0, fon.x);
    // col = mix(col, col2, pst);      
    //fix problems


    
    
    float l = length(p);
    if (l > r ) 
    {
         fragColor = vec4(col, 1.0);
         return;
    }
    //get vector on sphere in Orthographic Camera
    float z = r*sin(acos(l/r));
    vec3 sp = vec3(p, z);

    sp = normalize(sp);
    sp = reflect(vec3(0.0, 0.0, -1.0), sp);


    vec3 sp_rot = rotateZ(-t)*rotateX(PI/2.0)*sp; //rotate sphere
    vec2 pt = lonlat(sp_rot); //get longitude and latitude
    //col = mix(col, texture(iChannel0, pt).rgb, 0.4);  //get color from texture
    col = texture(iChannel0, pt).rgb;
    //fix problems
    // col2 = texture(iChannel0, vec2(0.0, pt.y)).rgb;
    // pst = smoothstep( 0.1, 0.0, pt.x); 
    // col = mix(col, col2, pst);      
    // pst = smoothstep( 0.9, 1.0, pt.x);
    // col = mix(col, col2, pst);      
    //fix problems

    
    //Lambertian
    float diffuse = dot(normalize(sp), normalize(light));
    diffuse = clamp(diffuse, 0.3, 1.0);
    //col *= diffuse;

    // gamma
    //col = pow( col, vec3(0.4545) );

    // float f = asin(pow(l/r, 2.0));
    // f =  smoothstep(0., PI/2.0, f);
    // vec3 rad = hsb2rgb(vec3(f,1.0,1.0));
    // col = mix(col, rad, 0.1);

    fragColor = vec4(col, 1.0);
    
}

/////=====================================================================================
void main()
{
    vec4 fragColor = vec4(0);
    mainImage(fragColor,  gl_FragCoord.xy);
    gl_FragColor = fragColor;
}