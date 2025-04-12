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

/*
braid rings
2D-SDF, tile, braid
simple braid pattern
*/
#define PI  3.14159265359
#define TAU 6.28318530718
#define rot(f) mat2(cos(f), -sin(f), sin(f), cos(f))

vec3 curColorN(float d1, vec3 col, vec3 col1, float dd)
{
    float s1 = smoothstep(0., -0.01, d1);
    float cs = dd*dd - (dd-abs(d1))*(dd-abs(d1));
    if (cs < 0.)
    {
        cs = 1.;
        d1 = 0.;
    }
    else
    {
        cs = sqrt(cs)/dd;  
        d1 = d1/dd;
    }
    vec3 norm = vec3(d1, 0., cs);
    vec3 light = vec3(0., 0., 1.);
    vec3 rd = vec3(0., 0., -1.);
    float difu = dot(norm , light);   
    vec3 R1 = reflect (light, norm);
    float shininess=25.0;
    float specular    =  pow(max(dot(R1, rd), 0.), shininess);
    vec3 colccs = col1*difu + 0.8*specular;      
    //return  (s1>0.) ? s1*col1*cs: col;
    return mix(col, colccs, vec3(s1));
}

vec3 vi5(vec2 p, vec3 col) {
    p.x *= TAU;
    vec3 col0 = vec3(0.5, 0.5, 1.);
    vec3 col1 = vec3(1., 0.5, 0.5);
    vec3 col2 = vec3(0.5, 1., 0.5);
    float dd = 0.15, k = 0.25;
    float fa0 = 0., fa1 = -2. / 3. * PI, fa2 = -4. / 3. * PI;
    vec3 facol0 = col0, facol1 = col1, facol2 = col2;
    if(p.x >= PI / 3.) {
        fa0 = -4. / 3. * PI, fa1 = 0., fa2 = -2. / 3. * PI;
        facol0 = col2, facol1 = col0, facol2 = col1;
    }
    if(p.x >= 2. * PI / 3.) {
        fa1 = -4. / 3. * PI, fa2 = 0., fa0 = -2. / 3. * PI;
        facol1 = col2, facol2 = col0, facol0 = col1;
    }
    if(p.x >= PI) {
        fa2 = -4. / 3. * PI, fa0 = 0., fa1 = -2. / 3. * PI;
        facol2 = col2, facol0 = col0, facol1 = col1;
    }
    if(p.x >= 4. * PI / 3.) {
        fa0 = -4. / 3. * PI, fa1 = 0., fa2 = -2. / 3. * PI;
        facol0 = col2, facol1 = col0, facol2 = col1;
    }
    if(p.x >= 5. * PI / 3.) {
        fa1 = -4. / 3. * PI, fa2 = 0., fa0 = -2. / 3. * PI;
        facol1 = col2, facol2 = col0, facol0 = col1;
    }
    float y = k * cos(p.x + fa0), alf = PI / 2. - (atan(-sin(p.x + fa0), 1.));
    float d1 = abs(p.y - y) - dd / sin(alf);
    col = curColorN(d1, col, facol0, dd / sin(alf));
    y = k * cos(p.x + fa1), alf = PI / 2. - (atan(-sin(p.x + fa1), 1.));
    d1 = abs(p.y - y) - dd / sin(alf);
    col = curColorN(d1, col, facol1, dd / sin(alf));
    y = k * cos(p.x + fa2), alf = PI / 2. - (atan(-sin(p.x + fa2), 1.));
    d1 = abs(p.y - y) - dd / sin(alf);
    col = curColorN(d1, col, facol2, dd / sin(alf));
    return col;
}

vec3 vi6(vec2 p)
{
    float h = 0.2, n = 10., d = TAU/n, r = floor(length(p)/h) * h + h/2.;
    r = max(r, 0.5);
    r = min(r, 0.9);
    if (r == 0.7)
        p.xy *= rot(iTime*0.1);

    float x = fract(mod(atan(p.y, p.x), TAU)/d);
    float y = (length(p) - r)/h;
    vec3 col = vec3(0.78, 0.78, 0.53);
    col = vi5(vec2(x, y), col );
    return col;

}
vec3 draw_braid(vec2 p, float r, vec3 col)
{
    float h = 0.2, n = 10., d = TAU/n;
    float x = fract(mod(atan(p.y, p.x), TAU)/d);
    float y = (length(p) - r)/h;
    //vec3 col = vec3(0.78, 0.78, 0.53);
    col = vi5(vec2(x, y), col );
    return col;
}

vec3 draw(vec2 p)
{
    float n = 3., d = TAU/n, alf = fract(mod(atan(p.y, p.x), TAU)/d)*TAU;
    float r = 0.7;// + 0.2*cos(alf);
    vec3 col = vec3(0.78, 0.78, 0.53);
    col = draw_braid(p, r, col);
    return col;
}


void mainImage(out vec4 fragColor, in vec2 fragCoord) {
    /*
    vec2 p = fragCoord / iResolution.xy;
    float n = 3.;
    float sign = (mod(floor(p.y*n), 2.0) - 0.5) * 2.0;
    p.x += sign * iTime*0.05;
    p = fract(p * vec2(3., n));
    p.y -= 0.5;
    vec3 col = vec3(0.78, 0.78, 0.53);
    col = vi5(p, col);
    */
    
    vec2 p = (2.0*fragCoord-iResolution.xy)/iResolution.y;
    vec3 col = draw(p); //vi6(p);
    
    
    fragColor = vec4(col, 1.0);
}
/////=====================================================================================
void main() {
    vec4 fragColor = vec4(0);
    mainImage(fragColor, gl_FragCoord.xy);
    gl_FragColor = fragColor;
}