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


struct HIT
{
    float dist;
    vec3 nor;
    vec3 pos;
};

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

//converts a vector on a sphere to longitude and latitude
vec2 lonlat (vec3 p)
{
    float lon = aafi(p.xy)/2.0/PI;
    float lat = aafi(vec2(p.z, length(p.xy)))/PI;
    return vec2(1.0-lon, lat);
}

const float dist_infin = 5.0;
const HIT hit_inf = HIT(100000.0, vec3(0.0), vec3(0.0));
#define nn 128
const float eps = 0.001;

/*
vec3 calcSkyReflect(vec3 rd, vec3 nor, mat3 sky)
{
    vec3 n = nor;
    float d = dot(rd, nor);
    n = nor*sign(d);
    vec3 r = reflect(rd, n);
    vec2 fon = lonlat(sky*r); //get longitude and latitude
    vec3 col = texture(iChannel0, fon).rgb;
    return col;

}
*/
vec3 culccolor(vec3 col_in, vec3 backcol, vec3 rd, vec3 light1, vec3 light2, vec3 nor)
{
    vec3 col = col_in;
    float d = dot(rd, nor);
    if (d < 0.0)
        col = backcol;
    
    nor *= -sign(d);
    float difu1 = dot(nor, light1);
    float difu2 = dot(nor, light2);
    float difu = max(difu1, difu2);
        col *= clamp(difu, 0.3, 1.0);
    return col;   
}

float near_dist(vec3 poin, vec2 tor)
{
  vec2 v1 = normalize(poin.xy)*tor.x;
  float d = length(poin - vec3(v1, 0.0)) - tor.y;
  return d;
}

// df(x)/dx
//analitic function https://www.shadertoy.com/view/4sBGDy
vec3 nTorus( in vec3 pos, vec2 tor )
{
	return normalize( pos*(dot(pos,pos)- tor.y*tor.y - tor.x*tor.x*vec3(1.0,1.0,-1.0)));
}



HIT giper3D(vec3 ro, vec3 rd, vec2 tor)
{
    float t  = 0.;
    
    for (int i = 0; i < nn; i++)
    {
      vec3 pos = ro + rd*t;
      float d = near_dist(pos, tor);
      if (d < eps)
      {
        vec3 nor = nTorus( pos, tor );
        return HIT(d, nor, pos);
      }
      t += d;
      if (t > dist_infin)
        return hit_inf;

    }
    return hit_inf;
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

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    //surface (x+y+z-a)(xy+yz+zx) - kxyz = 0
    vec3 light = normalize(vec3(0.0, 0.0, -1.0)); //light
    vec3 light2 = normalize(vec3(0.0, 0.0, 1.0)); //light

    

    float t = iTime;
    vec2 m = vec2(0.0, 0.0);
    //if  (iMouse.z > 0.0)
    {
        m = (-iResolution.xy + 2.0*(iMouse.xy))/iResolution.y;
        //t = 0.;
    }
    vec3 ro = vec3(0.0, 0.0, 2.5); // camera
    ro = rotateY(-m.x*TAU)*rotateX(-m.y*PI)*ro; //camera rotation
    
    
    const float fl = 1.5; // focal length
    float dist = dist_infin;
    float fi = PI/4.5;
    
    mat3 rota  = rotateZ(t)*rotateY(-t);
    mat3 rota_1  = rotateY(t)*rotateZ(-t);
    //mat3 sky = rotateZ(0.0)*rotateX(PI/2.0);
    
    vec2 tor = vec2(1.0,0.3);
    vec3 col = vec3(0.7, 0.7, 0.9); // background        

    vec3 tot = vec3(0.0);
    
    #define AA 3
    //antiblick
    for( int m=0; m<AA; m++ )
    for( int n=0; n<AA; n++ )
    {
        vec3 col = vec3(0.7, 0.7, 0.9); // background        
            // pixel coordinates
        vec2 o = vec2(float(m),float(n)) / float(AA) - 0.5;
        vec2 p = (-iResolution.xy + 2.0*(fragCoord+o))/iResolution.y;
        //vec2 p = (2.0*fragCoord-iResolution.xy)/iResolution.y;

        //vec3 rd = normalize( vec3(p,fl) ); // ray direction
        vec3 rd = GetRayDir(p, ro, vec3(0,0.,0), fl); //ray direction
        HIT giper = giper3D(rota*ro, rota*rd, tor);
        if (giper.dist < dist)
        {
           col = vec3(0.5, 0.5, 1.0);
            vec3 backcol = vec3(1.0, 0.2, 0.2);
            vec3 nor = rota_1*giper.nor;
            col = culccolor(col, backcol, -rd, light, light2, nor);
            // gamma
            //col = pow( col, vec3(0.4545) ); 
            //reflect
            //col = calcSkyReflect(-rd, nor, sky);
        }
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