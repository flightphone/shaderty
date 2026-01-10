#iChannel1 "file://../img/street.jpg"


float mlen(vec2 d, float aspect)
{
    float r = (d.x * d.x) + d.y * d.y /(aspect * aspect);
    return r;
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    
    float r = 0.2; //magnifying radius
    float a = .5; //magnifying zoom, may be negative
    
    float aspect = iResolution.x / iResolution.y;
    
    
    vec4 mc = vec4(0.0, 0.0, 1.0, 1.0);
    
    vec2 uv = fragCoord/iResolution.xy;
    vec2 pos = vec2(0.5, 0.5); 
    vec2 mouse = vec2(iMouse)/vec2(iResolution);
    if (!(mouse.x == 0.0 || mouse.y == 0.0))
        {pos = mouse;}
    else
        {pos = pos + vec2(0.2*cos(iTime), 0.2*sin(iTime));}
    
    
    vec2 d = uv - pos;
    float pst = step(mlen(d, aspect), r*r) * 0.08;
    vec4 c = texture(iChannel1, uv);
    if (pst > 0.0)
    {
        
        d = d*a;
        d = d + pos;
        if (d.x < 0.0 || d.x >1.0 || d.y < 0.0 || d.y >1.0)
        {
            c = vec4(1.0, 1.0, 1.0, 1.0);
        }
        else
        {
            c = texture(iChannel1, d);
        }
        c = mix(c, mc, pst);
    }
    
    fragColor = c;
}