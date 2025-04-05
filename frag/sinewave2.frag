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
stained glass windows
2D-SDF, tile, noise
stained glass windows
*/
#define PI  3.14159265359
#define TAU 6.28318530718

vec3 curColor(float d1, vec3 col, vec3 col1, float dd) {
    float s1 = smoothstep(0., -0.003, d1);
    float cs = dd * dd - (dd - abs(d1)) * (dd - abs(d1));
    if(cs < 0.)
        cs = 1.;
    else
        cs = sqrt(cs) / dd;
    vec3 colccs = col1 * cs;      
    return mix(col, colccs, vec3(s1));
}

vec3 vi5(vec2 p) {
    p.x *= TAU;
    p.y -= 0.5;
    vec3 col = vec3(0.78, 0.78, 0.53);
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
    col = curColor(d1, col, facol0, dd / sin(alf));
    y = k * cos(p.x + fa1), alf = PI / 2. - (atan(-sin(p.x + fa1), 1.));
    d1 = abs(p.y - y) - dd / sin(alf);
    col = curColor(d1, col, facol1, dd / sin(alf));
    y = k * cos(p.x + fa2), alf = PI / 2. - (atan(-sin(p.x + fa2), 1.));
    d1 = abs(p.y - y) - dd / sin(alf);
    col = curColor(d1, col, facol2, dd / sin(alf));
    return col;
}

void mainImage(out vec4 fragColor, in vec2 fragCoord) {
	//vec2 p = (-iResolution.xy + 2.0 * fragCoord) / iResolution.xy;
    vec2 p = fragCoord / iResolution.xy;
    float n = 3.;
    float sign = (mod(floor(p.y*n), 2.0) - 0.5) * 2.0;
    p.x += sign * iTime*0.05;
    p = fract(p * vec2(3., n));
    vec3 col = vi5(p);
    fragColor = vec4(col, 1.0);
}
/////=====================================================================================
void main() {
    vec4 fragColor = vec4(0);
    mainImage(fragColor, gl_FragCoord.xy);
    gl_FragColor = fragColor;
}