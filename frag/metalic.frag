//used https://www.shadertoy.com/view/mtdyD2 as a short template

#define rot2(a)      mat2(cos(a+vec4(0,11,33,0)))                 // rotation
#define TAU 6.283185

float sdCosN(vec3 p, float a, float n) {
    
    float fi = atan(p.y, p.x);//aafi(p.xy);
    float L = length(p.xy);
    float d = 10.;
    fi += step(p.y, 0.0)*TAU;
   
    for (float i = 0.; i < 5.; i++)
    {   
        //if (fract(n*i) == 0. && i > 0.)
        //    break;  
        float r = a * cos(n*fi + n*i*TAU);
        d = min(abs(L - r), d);
    }
    float f = acos(L/a)/n;
    for (float i = 0.; i < 6.; i++)
    {
        d = min(2.0 * abs(sin((fi - (f + i*TAU/n)) / 2.0)) * L, d);    
        d = min(2.0 * abs(sin((fi + (f + i*TAU/n)) / 2.0)) * L, d);    
        d = min(2.0 * abs(sin((fi - (f - i*TAU/n)) / 2.0)) * L, d);    
        d = min(2.0 * abs(sin((fi + (f - i*TAU/n)) / 2.0)) * L, d);    
    }

    d = sqrt(p.z * p.z  + d * d);
    d *= .7;
    d -= 0.04;
    return d;
}

float sdf(vec3 pos) {
    return sdCosN(pos, 1.0, 2.2);

}

void mainImage(out vec4 O, vec2 U)
{
    float t=9.,r;
    vec3  R = iResolution, e = vec3(1e-3,0,0), N,
          D = normalize(vec3(U+U, -18.*R.y) - R),          // ray direction
          p = vec3(0,0,20), q,                             // marching point along ray 
          C = iMouse.z > 0. ? 8.*iMouse.xyz/R -4.             // camera control
                          : 3.* cos(.3*iTime + vec3(0,11,0)); // demo mode
    p.yz *= rot2(-C.y),                                    // rotations
    p.xz *= rot2(-C.x-1.57),
    D.yz *= rot2(-C.y),                              
    D.xz *= rot2(-C.x-1.57);
    for ( O=vec4(1); O.x > 0. && t > .01; O-=.01 )         // march scene
        q = p, 
        t = min(t, sdf(q) ),                               // solenoïd
        p += t*D;                                       // step forward = dist to obj    
    
    N = vec3( sdf(q+e), sdf(q+e.yxy), sdf(q+e.yyx) ) - t ; // normal
    O.x < 0. ? O =  vec4(0.5, 0.5, 0.7, 1.0) :              // uncomment to display environment 
    O = texture(iChannel0, reflect(D,N/length(N) ) ); // reflect of environment map
}

