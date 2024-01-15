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


float sdCosN(vec2 p, float a, float n) {

    float fi = atan(p.y, p.x);
    float L = length(p.xy);
    float d1 = 100.;
    float d2 = 100.;

    for(float i = 0.; i < 10.; i++) {
        if(fract(n * i) == 0. && i > 0.)
            break;
        float r = a * cos(n * fi + n * i * TAU);
        d1 = min(abs(L - r), d1);
    }
    float f = acos(L / a) / n;
    for(float i = 0.; i < 10.; i++) {

        d2 = min(2.0 * abs(sin((fi - (f + i * TAU / n)) / 2.0)) * L, d2);
        d2 = min(2.0 * abs(sin((fi + (f + i * TAU / n)) / 2.0)) * L, d2);
        d2 = min(2.0 * abs(sin((fi - (f - i * TAU / n)) / 2.0)) * L, d2);
        d2 = min(2.0 * abs(sin((fi + (f - i * TAU / n)) / 2.0)) * L, d2);

    }
    
    float d = min(d1, d2);
    return d;
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
   /*
    vec3 col1 = vec3(0., 0., 0.15);
    float np[3]; float lv[3]; float ep[3];
    np[0] = 200.; np[1] =  100.; np[2] = 50.; 
    lv[0] = 0.7; lv[1] = .6; lv[2] = .5; 
    ep[0] = 0.1; ep[1] = 0.05; ep[2] = 0.01; 
   */
    vec3 col1 = vec3(0., 0., 0.3);
    float np[3]; float lv[3]; float ep[3];
    np[0] = 200.; np[1] =  100.; np[2] = 50.; 
    lv[0] = 0.8; lv[1] = .9; lv[2] = .9; 
    ep[0] = 0.02; ep[1] = 0.05; ep[2] = 0.1; 
 

    vec2 p = fragCoord / iResolution.y, c = vec2(iResolution.x/iResolution.y*0.5, 0.5), pc = p-c;
    vec3 tot = col1;
    float r = 0.8;
    float w = 1., v = 0.8, fi = mod(atan(pc.y, pc.x), TAU), head = w*iTime;
    float shift = mod((mod(head, TAU) - fi), TAU)/w;
    if (head < TAU && head < fi)
    {
        fragColor = vec4(col1, 1.0);
        return;
    }
    float d = sdCosN(pc, r, 2.);
    for (int i = 0; i < 3; i++)
    {
        float npp = np[i], level = lv[i], rd = 1./npp;  
        vec2 pp = floor(p*npp)/npp;
        //d > -ep[i]*shift*v && 
        if (d < ep[i]*shift*v)
        {
            float fil = hash21(pp);
            if (fil > level)
            {
                float pst = star(p, pp+rd*0.5, rd*.5);
                pst *= exp(-shift*0.7)*(1.+cos(shift*TAU*2. + fract(fil*500.)*3.*PI))/2.;
                vec3 col = hsb2rgb(vec3(fract(fil*1000.)*3., 1., 2.));
                tot = mix(tot, col, pst);
                
            }
        }
    }
  
    fragColor = vec4(tot, 1.0);
}
    
/////=====================================================================================
void main() {
    vec4 fragColor = vec4(0);
    mainImage(fragColor, gl_FragCoord.xy);
    gl_FragColor = fragColor;
}