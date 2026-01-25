#define PI  3.14159265359
#define TAU 6.28318530718

//@FabriceNeyret2  comment
#define S(y) smoothstep(fwidth(y), 0., abs(fract(y)-0.5) )


float FF(vec2 p)
{
    /*
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
    if (res >= 1.)
        res = log(res)/3.;
    return res;    
    */
    //@FabriceNeyret2  comment
    vec2 a = vec2(0, .5);  
    float l1 = length(p-a),
          l2 = length(p+a),
         eps = .001, res = 1./eps;
    if (l1 > eps && l2 > eps)
        res = abs( 1./l1 - 1./l2 ); 
    res = log((1. +  res));  //replace log(res)/3. to log((1. +  res));
    return res;    
}


float EE(vec2 p)
{
    /*
    float eps = 0.001;
	vec2 a = vec2(0.0, 0.5);  
    vec2 b = vec2(0.0, -0.5);
    vec2 l1 = normalize(p-a);
    vec2 l2 = normalize(p-b);
    vec2 z = normalize(a-b);
    float cos1 = dot(z, l1);
    float cos2 = dot(z, l2);
    float res = (cos1 - cos2);    
    return res;
    */
    //@FabriceNeyret2  comment
    vec2 a = vec2(0, .5);  
    vec2 l1 = normalize(p-a),
         l2 = normalize(p+a);
    return l1.y - l2.y;
}


/*
float dFF (vec2 p)
{
    float eps = 0.001;
    float y = (FF(p) - FF(vec2(p.x, p.y + eps)))/eps;
    float x = (FF(p) - FF(vec2(p.x+eps, p.y)))/eps;
    float res = sqrt(x*x + y*y);
    return res;
}


float dEE (vec2 p)
{
    float eps = 0.001;
    float y = (EE(p) - EE(vec2(p.x, p.y + eps)))/eps;
    float x = (EE(p) - EE(vec2(p.x+eps, p.y)))/eps;
    float res = sqrt(x*x + y*y);
    return res;
}
*/

vec3 lines(vec2 p)
{
	/*
    vec3 col0 = vec3(1., 1., 1.),colline = vec3(1., 0., 0.),colline2 = vec3(0., 0., 1.);
    float df = dFF(p);
    float de = dEE(p);

	float y = fract(FF(p)*8.); 
    float h = 0.03*df; 
    float eps = 15./iResolution.y*df; 
    float s1 = smoothstep(1. - h - eps, 1.-h, y),	
	s2 = smoothstep(h, h-eps, y);
	vec3 col1 = mix(col0, colline, s1);
	col1 = mix(col1, colline, s2);

    y = fract(EE(p)*10. - iTime ); 
    h = 0.04*de; 
    eps = 15./iResolution.y*de; 
    s1 = smoothstep(1. - h - eps, 1.-h, y);	
	s2 = smoothstep(h, h - eps, y);
	vec3 col2 = mix(col0, colline2, s1);
	col2 = mix(col2, colline2, s2);
	return mix(col1, col2, 0.5);	
    */
    //@FabriceNeyret2  comment
    float y = FF(p)*8. + iTime;	
	vec3 col1 = mix(vec3(1), vec3(1,0,0), S(y) );

    y = EE(p)*10. - iTime; 
	vec3 col2 = mix(vec3(1), vec3(0,0,1), S(y) );
	
	return mix(col1, col2, .5);	
    
}


void mainImage(out vec4 fragColor, in vec2 fragCoord) {
	vec2 p = (-iResolution.xy + 2.0 * fragCoord) / iResolution.y;
    if  (iMouse.z > 0.0)
    {
        vec2 mo = (-iResolution.xy + 2.0 * (iMouse.xy)) / iResolution.y;
        p -= mo;
    }
	vec3 col =  lines(p);
	fragColor = vec4(col, 1.0);
}

