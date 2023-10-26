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
const float nn = 6.0;
const vec3 fon = vec3(0., 0., 0.3);
const vec3 bound = vec3(1., 1., 0.);
const vec3 ball = vec3(1., 0.8, 0.8);

float rand(float t, float n, float shift)
{
    float res =  fract(sin(t + shift)*n);
    res = -1.0 + 2.0*res;
    return res;
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    vec2 uv = (2.0*fragCoord-iResolution.xy)/iResolution.y;
    float aspect = iResolution.x/iResolution.y;
    float tot = 0.;
    vec3 col = fon;
    for (float i = 10.0; i < 30.0; i++)
    {
        
        if (i > 10.0 + nn)
            break;
        float n1 = fract(sin(i)*213313.343)*100000.0;
        float n2 = fract(sin(i)*907629.243)*100000.0;
        
        float scale1 = fract(sin(i)*34545.56243)*5.0+4.0;
        float scale2 = fract(sin(i)*23554.243)*5.0 + 4.0;

        float i1 = floor(iTime/scale1/4.0);  // integer
        float f1 = fract(iTime/scale1/4.0);  // fraction

        float i2 = floor(iTime/scale2);  // integer
        float f2 = fract(iTime/scale2);  // fraction
        

        float r = 0.1*fract(sin(i)*25465.55) + 0.05;
        float x = mix(rand(i1, n1, 5.0), rand(i1+1.0, n1, 5.0), smoothstep(0.,1.,f1));
        float y = mix(rand(i2, n2, 15.0), rand(i2+1.0, n2, 15.0), smoothstep(0.,1.,f2));
        x*=aspect;
        y*=aspect;

        float t = iTime/r/5.;
        float e = 1.0 + 0.3*cos(t);
        float a = r * e;
        float b = r / e;
        float c = sqrt(abs(a*a - b*b));
        float s = 2.0*a;

        vec2 fc1 = vec2(x-c, y);
        vec2 fc2 = vec2(x+c, y);
        if (a < b)
        {
            fc1 = vec2(x, y-c);
            fc2 = vec2(x, y+c);
            s = 2.0*b;
        }
        float l = distance(uv, fc1) + distance(uv, fc2);
        tot += (s/l);
    }
    
    float pst1 = smoothstep(0.98, 1., tot);
    float pst2 = smoothstep(0.98, 1.05, tot);
    col = mix(col, bound, pst1);
    col = mix(col, ball, pst2);
    fragColor = vec4(col, 1.0);
    
}

/////=====================================================================================
void main()
{
    vec4 fragColor = vec4(0);
    mainImage(fragColor,  gl_FragCoord.xy);
    gl_FragColor = fragColor;
}