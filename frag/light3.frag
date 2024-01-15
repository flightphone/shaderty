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

#define PI  3.14159265359
#define TAU 6.28318530718


// IQ's vec2 to float hash.
float hash21(vec2 p){  
    return fract(sin(mod(dot(p, vec2(27.609, 57.583)), 6.2831853))*43758.5453); 
}

float hash1(float t)
{
    return  fract(sin(t)*5347.3987);
}

float star(vec2 p, vec2 c, float r)
{
    if (abs(p.x - c.x) > r || abs(p.y - c.y) > r)
        return 0.;
    float pst = smoothstep(-0.001, 0.0, length(p - (c + vec2(r, r))) - r);
    pst *=  smoothstep(-0.001, 0.0, length(p - (c + vec2(r, -r))) - r);   
    pst *=  smoothstep(-0.001, 0.0, length(p - (c + vec2(-r, r))) - r);   
    pst *=  smoothstep(-0.001, 0.0, length(p - (c + vec2(-r, -r))) - r);   
    return pst;
}



vec3 hsb2rgb( in vec3 c )
{
    vec3 rgb = clamp(abs(mod(c.x*6.0+vec3(0.0,4.0,2.0),
                             6.0)-3.0)-1.0,
                     0.0,
                     1.0 );
    rgb = rgb*rgb*(3.0-2.0*rgb);
    return (c.z * mix( vec3(1.0), rgb, c.y));
}

void mainImage(out vec4 fragColor, in vec2 fragCoord) 
{
    vec3 col1 = vec3(0., 0., 0.3);
    vec3 tot = col1;
    
    vec2 p = fragCoord/iResolution.xy;
    float n = 2000.;
    float k = 30.;
    vec3 c = texture(iChannel0, vec2(fract(iTime/100.), p.y)).rgb;
    float pst = smoothstep(0., 1., c.g);
    //float pst = c.g;
    //pst = floor(pst*k)/k;
    //if (hash1(y) > 0.8)
    pst *=exp(-3.*p.x/(pst)); //smoothstep(0.5/n, 0., abs(y-p.y))*
    //pst *=exp(-3.*p.x);
    tot = mix(tot, vec3(1.0), pst);
    
    /*
    vec2 p = (-iResolution.xy + 2.0 * fragCoord) / iResolution.y;
    float n = 150., dt = TAU/n, ep = 0.01;
    float t = float(iTime/300000.);
    float fi = mod(atan(p.y, p.x), TAU), f0 = floor(fi/dt)*dt, f1 = mod(f0 + dt, TAU);
    float pst = 0., level = 0.8;
    if (hash21(vec2(f0, t)) > level)
    {
        vec2 v0 = vec2(cos(f0), sin(f0));
        float d = length(p - v0 * abs(dot(p, v0)));
        pst = smoothstep(ep, 0., d)*exp(-6.*length(p));
        tot = mix(tot, vec3(1.0), pst);
    }
    
    if (hash21(vec2(f1, t)) > level)
    {
        vec2 v1 = vec2(cos(f1), sin(f1));
        float d = length(p - v1 * abs(dot(p, v1)));
        pst = smoothstep(ep, 0., d)*exp(-6.*length(p));;
        tot = mix(tot, vec3(1.0), pst);
    }
    */

    fragColor = vec4(tot, 1.0);
}
    
/////=====================================================================================
void main() {
    vec4 fragColor = vec4(0);
    mainImage(fragColor, gl_FragCoord.xy);
    gl_FragColor = fragColor;
}