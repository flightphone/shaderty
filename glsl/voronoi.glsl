vec2 hash( vec2 p )
{
    p = vec2(dot(p,vec2(127.1,311.7)),
             dot(p,vec2(269.5,183.3)));
    return fract(sin(p)*18.5453);
}


void mainImage(out vec4 fragColor, in vec2 fragCoord) {
	vec2 p = (fragCoord) / iResolution.y;
	float nn = 7.;
    p*= nn;
    

    vec2 n = floor( p );

    //ищем ближайшую точку
    float res = 100.;
    
    for (int i = -1; i <= 1; i++)
    for (int j = -1; j <= 1; j++)
    {
        vec2  g = vec2( float(i), float(j) );
        vec2  o = hash( n + g );
        vec2  r = n + g  + o - p; 
	    
		float d = dot( r, r );
        res = min(res, d);
        
    }
    vec3 col = vec3(res/2.);	
    fragColor = vec4(col, 1.0);
}

