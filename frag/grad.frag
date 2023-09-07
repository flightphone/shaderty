#ifdef GL_ES
precision mediump float;
#endif

uniform vec2 u_resolution;
uniform vec2 u_mouse;
uniform float u_time;
uniform sampler2D u_tex0;

#define iResolution u_resolution
#define iTime u_time
#define iChannel0 u_tex0
#define texture texture2D

/////=====================================================================================
//https://iquilezles.org/articles/distfunctions2d/
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

//=============================2D pic=========================
float sdBuild(vec2 p, float[30] f, float[30] r, int n)
{
    float fi = aafi(p);
    
    for (int i = 1; i < 30; i++)
    {
        if (i == n)
            break;
        if (fi >= f[i-1] && fi <= f[i])
        {
            float f1 = fi - f[i-1];
            float f2 = f[i] - fi;
            float h1 = r[i-1] * sin(f1);
            float h2 = r[i] * sin(f2);
            vec2 v1 = vec2(r[i-1] * cos(f[i-1]), r[i-1] * sin(f[i-1]));
            vec2 v2 = vec2(r[i] * cos(f[i]), r[i] * sin(f[i]));
            vec2 v3 = v2 - v1;
            v3 *= h1 / (h1 + h2);
            vec2 v4 = v1 + v3;
            return length(p) - length(v4);
        }
    }
}

vec3 grid(vec2 uv, float n)
{
        float rombr[30]; 
        rombr[0] = 0.5; 
        rombr[1] = 0.5; 
        rombr[2] = 0.5; 
        rombr[3] = 0.5; 
        rombr[4] = 0.5; 

        float rombf[30]; 
        rombf[0] = 0.0; 
        rombf[1] = PI/2.0;
        rombf[2] = PI;
        rombf[3] = 3.0*PI/2.0;
        rombf[4] = 2.0*PI;

    
        //vec2 p =  vec2(fract(uv.x*20.0), fract(uv.y*10.0));
        vec2 p = uv;
        //p = (p-0.5)*2.0;
        
        vec3 col = vec3(0.0);
        float pst = sdBuild(p, rombf, rombr, 5);
        //col = mix(col, vec3(1.0, 1.0, 1.0), step(pst, 0.0));
        if (pst > 0.)
        {
            pst  = clamp(exp(-100.0*abs(cos(iTime/10.0))*pst), 0.0, 1.0);
            col = mix(col, vec3(1.0, 1.0, 1.0), pst);
        }
        else
            col = vec3(1.0, 1.0, 1.0);
        /*
        for (int i = 0; i < 5; i++)
        {
            rombr[i] = 0.8;
        }
        pst = sdBuild(p, rombf, rombr, 5);
        col = mix(col, vec3(0.0, 0.0, 1.0), step(pst, 0.0));

        
        for (int i = 0; i < 5; i++)
        {
            rombr[i] = 0.6;
        }
        pst = sdBuild(p, rombf, rombr, 5);
        col = mix(col, vec3(1.0, 0.0, 0.0), step(pst, 0.0));

        for (int i = 0; i < 5; i++)
        {
            rombr[i] = 0.4;
        }
        pst = sdBuild(p, rombf, rombr, 5);
        col = mix(col, vec3(0.0, 1.0, 0.0), step(pst, 0.0));
        */
        // float grad = smoothstep(0.2, 1.2, length(p));
        // col = mix(col, vec3(0.0,0.0,0.0), grad);
        return col;
}
//=============================2D pic=========================

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    float r = 0.7;
    vec3 col = vec3(0.7, 0.7, 0.9);
    vec3 light = vec3(-10.0, 0.0, 10.0);
    
    vec2 p = (2.0*fragCoord-iResolution.xy)/iResolution.y;
    col = grid(p, 2.0);
    
    
    fragColor = vec4(col, 1.0);
    
}

/////=====================================================================================
void main()
{
    vec4 fragColor = vec4(0);
    mainImage(fragColor,  gl_FragCoord.xy);
    gl_FragColor = fragColor;
}