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
#define TAU 6.283185
//used https://www.shadertoy.com/view/mtdyD2 as a short template

#define rot2(a)      mat2(cos(a+vec4(0,11,33,0)))                 // rotation




//https://www.shadertoy.com/view/tt23RR
float sdHexagram( in vec2 p, in float r )
{
    const vec4 k = vec4(-0.5,0.86602540378,0.57735026919,1.73205080757);
    
    p = abs(p);
    p -= 2.0*min(dot(k.xy,p),0.0)*k.xy;
    p -= 2.0*min(dot(k.yx,p),0.0)*k.yx;
    p -= vec2(clamp(p.x,r*k.z,r*k.w),r);
    return length(p)*sign(p.y);
}

//Extrussion
//https://www.shadertoy.com/view/4lyfzw
float sdHexagram3( in vec3 p, in float h, in float r )
{
    float d = sdHexagram(p.xy, r);
    vec2 w = vec2( d, abs(p.z) - h );
    return min(max(w.x,w.y),0.0) + length(max(w,0.0)) - 0.05;
}

float sdf(vec3 pos) {
    return sdHexagram3(pos, 0.15, 0.35);

}

float aafi(vec2 p) {
    float fi = atan(p.y, p.x);
    fi += step(p.y, 0.0)*TAU;
    return fi;
}


//converts a vector on a sphere to longitude and latitude
vec2 lonlat (vec3 p)
{
    float lon = aafi(p.xy)/2.0/PI;
    float lat = aafi(vec2(p.z, length(p.xy)))/PI;
    return vec2(1.0-lon, lat);
}

vec4 csky(vec3 p)
{
    vec2 fon = lonlat(p); //get longitude and latitude
    return texture(iChannel0, fon);
}

void mainImage(out vec4 O, vec2 U)
{
    float t=9.,r;
    vec3  R = vec3(iResolution, 0); 
    vec3  e = vec3(1e-3,0,0), 
          N,
          D = normalize(vec3(U+U, -18.*R.y) - R),          // ray direction
          p = vec3(0,0,15), q,                             // marching point along ray 
          //C = iMouse.z > 0. ? 8.*iMouse.xyz/R -4.             // camera control
          //                : 3.* cos(.3*iTime + vec3(0,11,0)); // demo mode
          C = 3.* cos(.3*iTime + vec3(0,11,0));
       // S = ceil( vec3( 5.*iTime+vec2(1,0),0)   );
    
    p.yz *= rot2(-C.y),                                    // rotations
    p.xz *= rot2(-C.x-1.57),
    D.yz *= rot2(-C.y),                              
    D.xz *= rot2(-C.x-1.57);
    //for ( O=vec4(1); O.x > 0. && t > .01; O-=.01 )       // march scene
    O=vec4(1);
    for (int i = 0; i < 100; i++)
    {
        O -= .01;
        if (O.x <= 0.)
            break;

        q = p, 
        t = min(t, sdf(q) ),                               // solenoïd
        p += .5*t*D;                                       // step forward = dist to obj    
        if (t <= 0.01)
          break;
    }

    N = vec3( sdf(q+e), sdf(q+e.yxy), sdf(q+e.yyx) ) - t ; // normal
    O.x < 0. ? O = .5*csky(D) :              // uncomment to display environment 
    O *= O*O*O*3.* csky(reflect(D,N/length(N))); // reflect of environment map
}



/////=====================================================================================
void main()
{
    vec4 fragColor = vec4(0);
    mainImage(fragColor,  gl_FragCoord.xy);
    gl_FragColor = fragColor;
}