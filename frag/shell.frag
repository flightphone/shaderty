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


void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    
    vec2 p = (2.0*fragCoord-iResolution.xy)/iResolution.y;
    float l = length(p);
    float r = 0.95;
    float pst = step(l, r);
    vec3 fon = vec3(0.9,0.6,0.3);
    vec3 cir = vec3(0.65,0.85,1.0);
    vec3 col = mix(fon, cir, pst);
	col = mix( col, fon, 1.0-smoothstep(0.0,0.01,abs(l-r)));

    float n = 20.0;
    float fi = aafi(p);
    float f = acos(l/(r*1.0));
    
    float fip = mod((fi - f + iTime), (2.0*PI));
    fip = floor(fip*n/2.0/PI);
    float pst3 = step(1.0, mod(fip, 2.0));
    col = mix(col, vec3(0.5, 1.0, 0.5), 0.5*pst*pst3);
    
    fip = mod((fi + f - iTime), (2.0*PI));
    fip = floor(fip*n/2.0/PI);
    pst3 = step(1.0, mod(fip, 2.0));
    col = mix(col, vec3(1.0, 0.5, 0.5), 0.5*pst*pst3);
    
    fragColor = vec4(col,1.0);
}

/////=====================================================================================
void main()
{
    vec4 fragColor = vec4(0);
    mainImage(fragColor,  gl_FragCoord.xy);
    gl_FragColor = fragColor;
}