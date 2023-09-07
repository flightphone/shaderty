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
const float sp = 0.5;
const float num_samples = 256.0;
const float hh = 1.0;
const float kstep = 0.2;
const float ksum = 10.0;
const float eatt = 10.0;

vec3 attenuation(vec3 c, float d)
{
    return c/exp(eatt*d);
}


vec3 radial_blur_filter(in vec2 origin, in vec2 point)
{
    vec4 fon = texture(iChannel0, point);
    vec3 col = fon.rgb * fon.a;
    vec3 light = vec3(0.0);
    vec2 v = point - origin;
    float l = length(v);
    v = normalize(v); 
    float dz =  hh / num_samples * ksum;
    float l1 = l/(1.0 + hh);
    float al = (l-l1);
    for(float s = 0.0; s < num_samples; s++)
    {
        float l2 = l1 + al*s/num_samples;
        vec2 p = origin + v*l2;
        
        vec3 c = texture(iChannel0, p).rgb;
        float pst = smoothstep(kstep, 1.0, c.r);
        c *= pst;

        float fi = atan(l2, 1.0);
        float d = (l - l2)/sin(fi);
        c = attenuation(c, d);
        light += c;
    }
    light *= dz;
    float pst = smoothstep(0., 1.0, col.r);  
	light*=(1.0 - pst);
    
    col += light;
    return col;
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    vec2 pt = fragCoord.xy/iResolution.xy;
    vec4 col = texture(iChannel0, pt);
    vec2 mouse = vec2(0.5 - cos(iTime*sp),0.3);
    col.rgb = radial_blur_filter(mouse, pt);
    fragColor = col;
    
}

/////=====================================================================================
void main()
{
    vec4 fragColor = vec4(0);
    mainImage(fragColor,  gl_FragCoord.xy);
    gl_FragColor = fragColor;
}