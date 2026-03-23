#define PI  3.14159265359
#define TAU 6.28318530718



float lens(vec2 p, vec2 a, vec2 b, float d)
{
    vec2 c = (a+b)/2.0;
    vec2 norm = (b-a);
    float l = length(norm)/2.;
    norm = normalize(vec2(-norm.y, norm.x));
    float x = d*l, y = l/d, r = (x+y)/2., h = r - x;
    vec2 c1 = c + norm*h, c2 = c - norm*h;
    float di = min((r - length(p-c1)), (r - length(p-c2)));
    float hh = 0.001;
    di = smoothstep(0., hh, di)*di;
    di = pow(di, 0.25);
    x = pow(x, 0.25);
    float res = di/x;
    return res;
}


void mainImage(out vec4 fragColor, in vec2 fragCoord) {
	vec2 p = (-iResolution.xy + 2.0 * fragCoord) / iResolution.y;
	vec2 a = vec2(0., 0.9), b = vec2(0.0, -0.9);
    float r0 = 0.95, res = 0.;
    for (float i = 0.0; i < 13.0; i++)
    {
        vec2 a = r0*vec2(cos(TAU/13.0 * i), sin(TAU/13.0 * i)), 
        b = r0*vec2(cos(TAU/13.0 * (i + 5.0)), sin(TAU/13.0 * (i + 5.0)));
        res = max(res, lens(p, a, b, 0.05));

    }
    
    vec3 col =  vec3(res);
	fragColor = vec4(col, 1.0);
}

