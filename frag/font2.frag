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

#define PI 3.14159265359
#define TAU 6.283185

float aafi(vec2 p) {
    float fi = atan(p.y, p.x);
    fi += step(p.y, 0.0)*TAU;
    return fi;
}


//converts a vector on a sphere to longitude and latitude
vec2 lonlat (vec3 p)
{
    float lon = aafi(p.xy)/2.0/PI;
    float lat = aafi(vec2(p.z, length(p.xy)))/PI;
    return vec2(1.0-lon, lat);
}
vec3 csky(vec3 p)
{
    vec2 fon = lonlat(p); //get longitude and latitude
    return texture(iChannel1, fon).rgb;
    //return texture(iChannel1, p).rgb;
}

const float dist_infin = 4.0;
#define nn 128
const float eps = 0.01;


float dot2( in vec3 v ) { return dot(v,v); }

vec4 texChar(float char, vec2 uv) {
    vec2 pt = uv/16.0;
    pt.x += char/16.0;
    pt.y += 12.0/16.0;
    return texture(iChannel0, pt);
}


float sdBox( in vec2 p, in vec2 b )
{
    vec2 d = abs(p)-b;
    return length(max(d,0.0)) + min(max(d.x,d.y),0.0);
    
}

float sdTextBox ( vec2 p, vec2 b, float char ) {
    float l = sdBox(p,b);
    vec2 pn = (p.xy / b.xy) * .5 + .5;
    float lt = (texChar(char, pn).w - 0.5);
    return max(lt,l); 
}


float sdBox3( in vec3 p, in vec2 b, float h, float char)
{
    float d = sdTextBox(p.xy, b, char);
    vec2 w = vec2( d, abs(p.z) - h );
    return min(max(w.x,w.y),0.0) + length(max(w,0.0)) - 0.02;
}

float map( in vec3 pos, float char )
{
    return sdBox3(pos, vec2(1., 1.), 0.2, char);
    
}

// https://iquilezles.org/articles/normalsSDF
vec3 calcNormal( in vec3 pos, float char )
{
    const float h = 0.0001; // replace by an appropriate value
    const vec2 k = vec2(1,-1);
    return normalize( k.xyy*map(pos + k.xyy*h, char ) + 
                      k.yyx*map(pos + k.yyx*h, char ) + 
                      k.yxy*map(pos + k.yxy*h, char ) + 
                      k.xxx*map(pos + k.xxx*h, char ) );
}

float steps(vec3 ro, vec3 rd, float char)
{
    float t  = 0.;
    for (int i = 0; i < nn; i++)
    {
        vec3 pos = ro + rd*t;
        float h = map(pos, char);
        if (h < eps || t >= dist_infin)
            break;
        t += h;  
    }    

    if (t >= dist_infin)
        return dist_infin;
    return t;    
}

vec3 GetRayDir(vec2 uv, vec3 p, vec3 l, float z) {
    vec3 
        f = normalize(l-p),
        r = normalize(cross(vec3(0,1,0), f)),
        u = cross(f,r),
        c = f*z,
        i = c + uv.x*r + uv.y*u;
    return normalize(i);
}

mat3 rotateX(float f)
{
    return mat3(
    vec3(1.0,    0.0,      0.0),
    vec3(0.0,	 cos(f),  -sin(f)), 	
	vec3(.0, sin(f), cos(f))
    );
}


mat3 rotateZ(float f)
{
    return mat3(
    vec3(cos(f),    -sin(f),  0.0),
    vec3(sin(f),	 cos(f),  0.0), 	
	vec3(0.0, 0.0, 1.0)
    );
    
}


mat3 rotateY(float f)
{
    return mat3(
    vec3(cos(f), 0.0,  sin(f)),
    vec3(0.0,	 1.0,  0.0), 	
	vec3(-sin(f), 0.0, cos(f))
    );
}


void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    float t = iTime;
    vec2 m = 3.* cos(.1*t + vec2(.5,.5)) + vec2(0.01, 0.01);
    //if  (iMouse.z > 0.0)
    //    m = (-iResolution.xy + 2.0*(iMouse.xy))/iResolution.y;
    vec3 ro = vec3(0.0, 0.0, 2.5); // camera
    ro = rotateY(-m.x*TAU)*rotateX(-m.y*PI)*ro; //camera rotation
    
    
    const float fl = 1.5; // focal length
    vec3 tot = vec3(0.0);
    
    #define AA 2
    //antiblick
    for( int m=0; m<AA; m++ )
    for( int n=0; n<AA; n++ )
    {
        float dist = dist_infin;
        vec2 o = vec2(float(m),float(n)) / float(AA) - 0.5;
        vec2 p = (-iResolution.xy + 2.0*(fragCoord+o))/iResolution.y;
        vec3 rd = GetRayDir(p, ro, vec3(0,0.,0), fl); //ray direction
        vec3 col = vec3(0.5, 0.5, 1.);//csky(rd);
        vec3 shift = vec3(1.5, 0., 0.);
        vec3 sh = vec3(0.0);
        float shx = 1.0;
        vec3 nor = vec3(0.);
        float char = 6.;
        float giper = steps((ro + shift), rd, 2.0);
        if (giper < dist)
        {
            dist = giper;
            char = 2.0;
            sh = shift;
        }
        
        shift.x -= shx;
        giper = steps((ro + shift), rd, .0);
        if (giper < dist)
        {
            dist = giper;
            char = 0.;
            sh = shift;
        }

        
        shift.x -= shx;
        giper = steps((ro + shift), rd, 2.);
        if (giper < dist)
        {
            dist = giper;
            char = 2.;
            sh = shift;
        }
        
        shift.x -= shx;
        giper = steps((ro + shift), rd, 4.);
        if (giper < dist)
        {
            dist = giper;
            char = 4.;
            sh = shift;
        }

        if (dist < dist_infin)
        {
            vec3 pos = (ro + sh) + dist*rd;
            vec3 nor = calcNormal(pos, char);
            col = csky(reflect(rd,normalize(nor)));
        }

        // gamma        
        //col = sqrt( col );
	    tot += col;
    }
    //antiblick
    tot /= float(AA*AA);
    fragColor = vec4(tot,1.0);
    
}

/////=====================================================================================
void main()
{
    vec4 fragColor = vec4(0);
    mainImage(fragColor,  gl_FragCoord.xy);
    gl_FragColor = fragColor;
}