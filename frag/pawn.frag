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
const vec3 fon = vec3(0., 0., 0.4);
const vec3 bound = vec3(1., 1., 0.);
const vec3 ball = vec3(1., 0.8, 0.8);

float rect(float x, float y, float a, float b, vec2 p)
{
    float h = max(2.0 * abs(x + a/2.0 - p.x) / a,  2.0 * abs(y + b/2. - p.y) / b);
    return 1./h/h/h;
}

float ellips (float x, float y, float a, float b, vec2 p)
{
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
        float l = distance(p, fc1) + distance(p, fc2);
        return (s*s*s*s/l/l/l/l);
}
void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    vec2 uv = (2.0*fragCoord-iResolution.xy)/iResolution.y;
    float aspect = iResolution.x/iResolution.y;
    float tot = 0.;
    vec3 col = fon;
    tot += rect(-0.4, -0.8, 0.8, 0.03, uv);
    tot += ellips(0., -0.6, 0.6, 0.15, uv);
    tot += rect(-0.1, -.4, 0.2, 0.4, uv);
    tot += ellips(0., 0.3, 0.4, 0.03, uv);
    tot += ellips(0., 0.68, 0.15, 0.15, uv);
    
    float pst1 = smoothstep(0.99, 1., tot);
    float pst2 = smoothstep(0.99, 1.1, tot);
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