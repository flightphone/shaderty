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
/*
national ornament
tile, pixel, patern
national ornament
*/
#define PI  3.14159265359
#define TAU 6.28318530718

vec3 draw_vert(vec2 p, vec2 a, float l, vec3 col, vec3 col1)
{
    float s = 0.;
    if (p.x == a.x && p.y - a.y + 1. <= l && p.y >= a.y)
        s = 1.;
    col = mix(col, col1, s); 
    return col;    
}
vec3 draw_hr(vec2 p, vec2 a, float l, vec3 col, vec3 col1)
{
    float s = 0.;
    if (p.y == a.y && p.x - a.x + 1. <= l && p.x >= a.x)
        s = 1.;
    col = mix(col, col1, s); 
    return col;    
}

vec3 draw_box(vec2 p, vec2 a, vec2 l, vec3 col, vec3 col1)
{
    float s = 0.;
    if (
        p.x - a.x + 1. <= l.x && p.x >= a.x &&
        p.y - a.y + 1. <= l.y && p.y >= a.y
    )
        s = 1.;
    col = mix(col, col1, s); 
    return col;    
}



vec3 ornament(vec2 p)
{
    vec2 xo = vec2(1., -1.), yo = vec2(1, 1);
    vec3 col0 = vec3(1.), col1 = vec3(1., 0., 0.), col = col0;
    float dlt = 1./55., x = floor(p.x/dlt), y = floor(p.y/dlt);
    vec2 pp = x*xo + y*yo;
    float n = 25., n2 = n*2.;
    x = mod(pp.x, n2);
    y = mod(pp.y, n2);
    
    
    col = draw_hr(vec2(x, y), vec2(0., 0.), n2, col, col1);
    col = draw_hr(vec2(x, y), vec2(0., n2-2.), n2, col, col1);
    col = draw_vert(vec2(x, y), vec2(0., 0.), n2, col, col1);
    col = draw_vert(vec2(x, y), vec2(n2 - 2., 0.), n2, col, col1);

    col = draw_box(vec2(x, y), vec2(4.), vec2(2.*(n-4.) - 1.), col, col1);
    
    
    col = draw_box(vec2(x, y), vec2(3., 13.), vec2(7.), col, col0);
    col = draw_box(vec2(x, y), vec2(3., 29.), vec2(7.), col, col0);

    col = draw_box(vec2(x, y), vec2(13., 3.), vec2(7.), col, col0);
    col = draw_box(vec2(x, y), vec2(29., 3.), vec2(7.), col, col0);

    col = draw_box(vec2(x, y), vec2(39., 13.), vec2(7.), col, col0);
    col = draw_box(vec2(x, y), vec2(39., 29.), vec2(7.), col, col0);

    col = draw_box(vec2(x, y), vec2(13., 39.), vec2(7.), col, col0);
    col = draw_box(vec2(x, y), vec2(29., 39.), vec2(7.), col, col0);
    

    col = draw_hr(vec2(x, y), vec2(0., 16.), 7., col, col1);
    col = draw_hr(vec2(x, y), vec2(41., 16.), 7., col, col1);
    col = draw_vert(vec2(x, y), vec2( 16., 0.), 7., col, col1);
    col = draw_vert(vec2(x, y), vec2( 16., 41.), 7., col, col1);
    

    col = draw_hr(vec2(x, y), vec2(0., 32.), 7., col, col1);
    col = draw_hr(vec2(x, y), vec2(41., 32.), 7., col, col1);
    col = draw_vert(vec2(x, y), vec2( 32., 0.), 7., col, col1);
    col = draw_vert(vec2(x, y), vec2( 32., 41.), 7., col, col1);


    col = draw_box(vec2(x, y), vec2(13.), vec2(2.*(n-13.) - 1.), col, col0);
    col = draw_box(vec2(x, y), vec2(20.), vec2(2.*(n-20.) - 1.), col, col1);
    col = draw_box(vec2(x, y), vec2(23.), vec2(3.), col, col0);
    col = draw_box(vec2(x, y), vec2(7., 23.), vec2(3.), col, col0);
    col = draw_box(vec2(x, y), vec2(23., 7.), vec2(3.), col, col0);
    col = draw_box(vec2(x, y), vec2(39., 23.), vec2(3.), col, col0);
    col = draw_box(vec2(x, y), vec2(23., 39.), vec2(3.), col, col0);
    col = draw_box(vec2(x, y), vec2(39., 7.), vec2(3.), col, col0);
    col = draw_box(vec2(x, y), vec2(7., 39.), vec2(3.), col, col0);
    col = draw_box(vec2(x, y), vec2(7., 7.), vec2(3.), col, col0);
    col = draw_box(vec2(x, y), vec2(39., 39.), vec2(3.), col, col0);


    col = draw_hr(vec2(x, y), vec2(16., 16.), 2.*(n-16.), col, col1);
    col = draw_hr(vec2(x, y), vec2(16., 32.), 2.*(n-16.), col, col1);
    col = draw_vert(vec2(x, y), vec2(16., 16.), 2.*(n-16.), col, col1);
    col = draw_vert(vec2(x, y), vec2(32., 16.), 2.*(n-16.), col, col1);
    //col = draw_vert(vec2(x, y), vec2(0., 0.), n2, col, col1);
    //col = draw_vert(vec2(x, y), vec2(n2 - 2., 0.), n2, col, col1);
    
    return col;

}

void mainImage(out vec4 fragColor, in vec2 fragCoord) {
	
    vec2 p = fragCoord / iResolution.y;
    p.y -= 0.05;
    p.x -= 0.2;
    //if  (iMouse.z > 0.0)
    {
        vec2 mo = (-iResolution.xy + 2.0 * (iMouse.xy)) / iResolution.y;
        p -= mo;
    }
    vec3 col = ornament(p);
    fragColor = vec4(col, 1.0);
}
/////=====================================================================================
void main() {
    vec4 fragColor = vec4(0);
    mainImage(fragColor, gl_FragCoord.xy);
    gl_FragColor = fragColor;
}