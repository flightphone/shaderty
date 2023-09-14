#version 300 es

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
//#define texture texture2D




/////=====================================================================================

#define PI 3.14159265359
#define TAU 6.283185


vec3 uhash3( uvec3 x ) {         // iq version https://shadertoy.com/view/XlXcW4
    const uint k = 1103515245U;  // GLIB C version
    x = ((x>>8U)^x.yzx)*k;
    x = ((x>>8U)^x.yzx)*k;
    x = ((x>>8U)^x.yzx)*k;    
    return vec3(x)/float(0xffffffffU);
}
#define iH3(f)  uhash3( floatBitsToUint(f) ) // FabriceN version https://www.shadertoy.com/view/NtjyWw
#define iH32(f) iH3(f).xy

/*
void mainImage( out vec4 O, vec2 u )
{
  vec2 pt = u/iResolution.xy;
  O = texture(iChannel0, pt);
}
*/

void mainImage( out vec4 O, vec2 u ) // ----------------------------------------------------------------------
{
    float d = 9., r,t, n;                                   // n x n cells
    vec2  R = iResolution.xy;
    r = R.x/R.y;                                            // aspect ratio
    t = floor( 16.*16. * ( 
                  // target point density
                             10.* ( .5+.5*sin(iTime) ) )  // mouse vs demo mode.
                          // : float(iFrame)
             );
    n = max(1., floor(sqrt(t/5.)) );                        // 5 = target point per cell. Choose n accordingly.
    t = round(t/(n*n))*n*n;                                 // if we want same number of points everywhere.
    
    vec2  U = n * u / R.y, D,
          I = floor(U),                                     // cell id
          F = U-I;                                          // coords in cells  
    O *= 0.;   
   
    if ( u.x < R.x/2. )                                     // --- left: white noise in n x n cells
      for (float i=0.; i < t/(n*n); i++ )
          D = F - iH32(vec3(I,i)),                          // distance vector to point
          d = min(d, dot(D,D) );                            // keep smallest d²
       // O += smoothstep(1.5*n/R.y, 0., length(D) );       // draw point
          
    else                                                    // --- right: white noise in whole screen.
      for (float i=0.; i < t*r/2.; i++ )                    // unstructured: we need to test all points
          D = U - n*( vec2(r/2.,0) + vec2(r/2.,1)* iH32(i*vec3(1,17,71)) ),
          d = min(d, dot(D,D) );
    
    O += smoothstep(1.5*n/R.y, 0., sqrt(d) );               // draw closest point
    
    if ( int(u)==int(R)/2 ) O = vec4(1,0,0,0);              // red separator
}

/////=====================================================================================



out vec4 FragColor;
void main()
{
    mainImage(FragColor,  gl_FragCoord.xy);
}