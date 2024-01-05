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
#define rot(f) mat2(cos(f), -sin(f), sin(f), cos(f))

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
    return d*0.3;
}


void mainImage(out vec4 fragColor, in vec2 fragCoord) 
{
    vec3 col3 = vec3(1., 1., 0.75);
    vec3 col1 = vec3(0., 0., 0.3);
    //col3 = pow(col3, vec3(.001));
    float np[3]; float lv[3]; float ep[3];
    np[0] = 200.; np[1] =  100.; np[2] = 50.; 
    lv[0] = 0.8; lv[1] = .9; lv[2] = .9; 
    ep[0] = 0.02; ep[1] = 0.05; ep[2] = 0.1; 

    

    vec2 p = fragCoord / iResolution.y, c = vec2(iResolution.x/iResolution.y*0.5, 0.5), pc = p-c;

    vec3 tot = col1;

    
    //float r = iResolution.x/iResolution.y*0.45;
    float r = 0.5;
    float w = 1., v = 0.3, fi = mod(atan(pc.y, pc.x), TAU), head = w*iTime;
    float shift = mod((mod(head, TAU) - fi), TAU)/w;
    if (head < TAU && head < fi)
    {
        fragColor = vec4(col1, 1.0);
        return;
    }
    //float d = length(pc) - r;
    //float pst = smoothstep(-0.01, 0., d) * smoothstep(0.01, 0., d);
    float d = sdCosN(pc, r, 3.5);
    float pst = smoothstep(0.002, 0., d);
    pst *= exp(-shift);
    tot = mix(tot, col3, pst);
    
    
    for (int i = 0; i < 3; i++)
    {
        float npp = np[i], level = lv[i], rd = 1./npp;  
        vec2 pp = floor(p*npp)/npp;
        if (d > -ep[i]*shift*v && d < ep[i]*shift*v)
        if (hash21(pp) > level)
        {
            float pst = star(p, pp+rd*0.5, rd*.5);
            pst *= exp(-shift);
            tot = mix(tot, col3, pst);
            
        }
    }
    //tot*=exp(-shift);
        
    
   
    
    
    //float pst = smoothstep(-0.01, 0., d) * smoothstep(0.01, 0., d);
    //float pst = star(p, vec2(0.0), r);
    //float pst = star(p, vec2(0.0), 0.3);
    //tot = mix(tot, vec3(1.), pst);
    
    fragColor = vec4(tot, 1.0);

    
}
    
/////=====================================================================================
void main() {
    vec4 fragColor = vec4(0);
    mainImage(fragColor, gl_FragCoord.xy);
    gl_FragColor = fragColor;
}