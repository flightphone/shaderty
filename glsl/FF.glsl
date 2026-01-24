#define PI  3.14159265359
#define TAU 6.28318530718


float FF(vec2 p)
{
	vec2 a = vec2(0.0, 0.5);  
    vec2 b = vec2(0.0, -0.5);
    float l1 = length(p-a);
    float l2 = length(p-b);
    float eps = 0.01, res = 1.0/eps;
    if (l1 > eps && l2 > eps)
    {
        res = 1.0/l1 - 1.0/l2;
    }
    res =  abs(res);
    if (res > .7)
        res = log(res)/3.;
    
    return res;    
}

float EE(vec2 p)
{
	vec2 a = vec2(0.0, 0.5);  
    vec2 b = vec2(0.0, -0.5);
    vec2 l1 = normalize(p-a);
    vec2 l2 = normalize(p-b);
    float cos1 = dot(vec2(0., 1.), l1);
    float cos2 = dot(vec2(0., 1.), l2);
    return cos1 - cos2;    
}


vec3 lines(vec2 p)
{
	
    vec3 col = vec3(1., 1., 1.),colline = vec3(1., 0., 0.),colline2 = vec3(0., 0., 1.);

	float y = fract(((FF(p)))*5.), h = 0.07, eps = 15./iResolution.y, s1 = smoothstep(1. - h - eps, 1.-h, y),	
	s2 = smoothstep(h, h-eps, y);
	col = mix(col, colline, s1);
	col = mix(col, colline, s2);

    y = fract(((EE(p)))*5.); 
    h = 0.05; 
    eps = 15./iResolution.y; 
    s1 = smoothstep(1. - h - eps, 1.-h, y);	
	s2 = smoothstep(h, h-eps, y);
	col = mix(col, colline2, s1);
	col = mix(col, colline2, s2);
	return col;	
}


void mainImage(out vec4 fragColor, in vec2 fragCoord) {
    //vec2 p = vec2(fragCoord.x/iResolution.x, fragCoord.y/iResolution.y); //(-iResolution.xy + 2.0 * fragCoord) / iResolution.y;
	vec2 p = (-iResolution.xy + 2.0 * fragCoord) / iResolution.y;
    //if  (iMouse.z > 0.0)
    {
        //vec2 mo = (-iResolution.xy + 2.0 * (iMouse.xy)) / iResolution.y;
        //p -= mo;
    }
	vec3 col =  lines(p);
	fragColor = vec4(col, 1.0);
}