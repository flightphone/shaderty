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
//https://iquilezles.org/articles/distfunctions2d/
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


float rand(float t, float n, float shift)
{
    float res =  fract(sin(t + shift)*n);
    res = -1.0 + 2.0*res;
    return res;
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

vec2 lonlat (vec3 p)
{
    float lon = aafi(p.xy)/2.0/PI;
    float lat = aafi(vec2(p.z, length(p.xy)))/PI;
    return vec2(1.0-lon, lat);
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    
    vec2 p = (2.0*fragCoord-iResolution.xy)/iResolution.y;
    float aspect = iResolution.x/iResolution.y;

    const float fl = 2.5;
    // ray direction
    vec3 rd = normalize( vec3(p,fl) );
    float t = iTime/7.0;
    mat3 rota = rotateZ(t)*rotateX(PI/2.0);
    mat3 rota1 = rotateZ(t+PI)*rotateX(PI/2.0);
    vec2 fon = lonlat(rota*rd); //get longitude and latitude
    vec3 col = texture(iChannel0, fon).rgb;  //get color from texture


    
    float nn = 10.0;
    float nshow = 0.0;
    for (float i = 10.0; i < 30.0; i++)
    {
        
        if (i > 10.0 + nn)
            break;
        float n1 = fract(sin(i)*213313.343)*100000.0;
        float n2 = fract(sin(i)*907629.243)*100000.0;
        
        float scale1 = fract(sin(i)*34545.56243)*5.0+4.0;
        float scale2 = fract(sin(i)*23554.243)*5.0 + 4.0;

        float i1 = floor(iTime/scale1);  // integer
        float f1 = fract(iTime/scale1);  // fraction

        float i2 = floor(iTime/scale2);  // integer
        float f2 = fract(iTime/scale2);  // fraction


        

        

        

        float r = 0.2*fract(sin(i)*25465.55) + 0.1;

        float vis = floor(rand(i1, n1, 8.0)*20.0);
        
        if (vis < 2.0 && !(i == 9.0 + nn && nshow == 0.0))
            continue;
        nshow += 1.0;

        float x = mix(rand(i1, n1, 5.0), rand(i1+1.0, n1, 5.0), smoothstep(0.,1.,f1));
        float y = mix(rand(i2, n2, 15.0), rand(i2+1.0, n2, 15.0), smoothstep(0.,1.,f2));
        x*=aspect;
        y*=aspect;
        float l = length(p - vec2(x,y));
        if (l > r)
            continue;

        float z = r*sin(acos(l/r));
        vec3 sp = normalize(vec3(p - vec2(x,y), z));
        sp = reflect(vec3(0.0, 0.0, -1.0), sp);
        vec2 fon1 = lonlat(rota1*sp); //get longitude and latitude
        col = mix(col, texture(iChannel0, fon1).rgb, 0.5);  //get color from texture
        //col = texture(iChannel0, fon1).rgb;

      
        float f = asin(pow(l/r, 2.0));
        f =  smoothstep(0., PI/2.0, f);
        vec3 rad = hsb2rgb(vec3(f,1.0,1.0));
        col = mix(col, rad, 0.1);
        
    }

    
    fragColor = vec4(col, 1.0);
    
}

/////=====================================================================================
void main()
{
    vec4 fragColor = vec4(0);
    mainImage(fragColor,  gl_FragCoord.xy);
    gl_FragColor = fragColor;
}