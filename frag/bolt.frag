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

float aafi(vec2 p)
{
    float l = length(p);
    float fi = asin(abs(p.y)/l);
    float pst = step(0.0, p.y)*step(p.x, 0.0);
    fi = fi + pst*(PI - 2.0*fi);
    pst = step(p.y, 0.0)*step(p.x, 0.0);
    fi = fi + pst*PI;
    pst = step(p.y, 0.0)*step(0.0, p.x);
    fi = fi + pst*(2.0*PI - 2.0*fi);
    return fi;    
}


void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    
    vec2 p = (2.0*fragCoord-iResolution.xy)/iResolution.y;
    float r = 0.2;
    float r2 = r*0.8;
    float v = 10.0;
    vec3 col = vec3(0.5,0.5, 1.0);
    if (abs(p.y) < r)
    {
        
        float pst = 0.0;
        vec2 p1 = vec2(0.0);
        vec2 p2 = vec2(0.0);
        float l = 0.0;
        float x = 0.0;
        float f = 0.0;
        float f2 = 0.0;
        float dx1 = pow((r*r - p.y*p.y), 0.5); 
        float dx2 = pow((r2*r2 - p.y*p.y), 0.5);
        float dx = (dx1 + dx2)/2.0; 
        float d = 0.0;

        
        
        
        p1 = vec2(dx1, p.y);
        f = mod(v*(p.x-dx1), 2.0*PI);
        
        

        p2 = vec2(dx2, p.y);
        f2 = mod(v*(p.x-dx2), 2.0*PI);
        
        if (aafi(p1) >= f && aafi(p2) <= f2)
            pst = 1.0;
        
        p1 = vec2(-dx1, p.y);
        f = mod(v*(p.x+dx1), 2.0*PI);
        p2 = vec2(-dx2, p.y);
        f2 = mod(v*(p.x+dx2), 2.0*PI);
        
        if (aafi(p1) <= f && aafi(p2) >= f2)
            pst = 1.0;

        /*
        p1 = vec2(dx2, p.y);
        f = v*(p.x-dx2);
        d = dot(normalize(p1), vec2(cos(f), sin(f)));
        pst += step(0.9, d) ;

        p1 = vec2(-dx2, p.y);
        f = v*(p.x+dx2);
        d = dot(normalize(p1), vec2(cos(f), sin(f)));
        pst += step(0.9, d) ;
        */

        /*
        p1 = p + vec2(dx1, 0);
        f = -v*p1.x;
        d = dot(normalize(p1), vec2(cos(f), sin(f)));
        pst += step(0.99, d) ;
        */
        /*
        x = p.x - dx1;
        f = -v*x;
        p2 = vec2(x + r*cos(f), r*sin(f));
        l = length(p-p2);
        pst += step(l, r*0.2);
        
        x = p.x - dx2;
        f = -v*x;
        p2 = vec2(x + r*cos(f), r*sin(f));
        l = length(p-p2);
        pst += step(l, r*0.2);
        
        x = p.x + dx1;
        f = -v*x;
        p2 = vec2(x + r*cos(f), r*sin(f));
        l = length(p-p2);
        pst += step(l, r*0.2);
        
        x = p.x + dx2;
        f = -v*x;
        p2 = vec2(x + r*cos(f), r*sin(f));
        l = length(p-p2);
        pst += step(l, r*0.2);
        */

        col = mix(col, vec3(1.0), pst);
        

    }
    
    fragColor = vec4(col,1.0);
}

/////=====================================================================================
void main()
{
    vec4 fragColor = vec4(0);
    mainImage(fragColor,  gl_FragCoord.xy);
    gl_FragColor = fragColor;
}