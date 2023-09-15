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

const float dist_infin = 100000.0;
const HIT hit_inf = HIT(100000.0, vec3(0.0), vec3(0.0));

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

vec3 culccolor(vec3 col_in, vec3 backcol, vec3 rd, vec3 light, vec3 nor)
{
    vec3 col = col_in;
    float d = dot(rd, nor);
    if (d < 0.0)
        col = backcol;
    
    nor *= -sign(d);
    float difu = dot(nor, light);
    col *= clamp(difu, 0.7, 1.0);
    return col;   
}



// f(x) = (|x|² + R² - r²)² - 4·R²·|xy|² = 0
//analitic function https://www.shadertoy.com/view/4sBGDy
float iTorus( in vec3 ro, in vec3 rd, in vec2 tor )
{
    float po = 1.0;
    
    float Ra2 = tor.x*tor.x;
    float ra2 = tor.y*tor.y;
	
    float m = dot(ro,ro);
    float n = dot(ro,rd);

    // bounding sphere
    {
	float h = n*n - m + (tor.x+tor.y)*(tor.x+tor.y);
	if( h<0.0 ) return -1.0;
	//float t = -n-sqrt(h); // could use this to compute intersections from ro+t*rd
    }
    
	// find quartic equation
    float k = (m - ra2 - Ra2)/2.0;
    float k3 = n;
    float k2 = n*n + Ra2*rd.z*rd.z + k;
    float k1 = k*n + Ra2*ro.z*rd.z;
    float k0 = k*k + Ra2*ro.z*ro.z - Ra2*ra2;
	
    #if 1
    // prevent |c1| from being too close to zero
    if( abs(k3*(k3*k3 - k2) + k1) < 0.01 )
    {
        po = -1.0;
        float tmp=k1; k1=k3; k3=tmp;
        k0 = 1.0/k0;
        k1 = k1*k0;
        k2 = k2*k0;
        k3 = k3*k0;
    }
	#endif

    float c2 = 2.0*k2 - 3.0*k3*k3;
    float c1 = k3*(k3*k3 - k2) + k1;
    float c0 = k3*(k3*(-3.0*k3*k3 + 4.0*k2) - 8.0*k1) + 4.0*k0;

    
    c2 /= 3.0;
    c1 *= 2.0;
    c0 /= 3.0;
    
    float Q = c2*c2 + c0;
    float R = 3.0*c0*c2 - c2*c2*c2 - c1*c1;
    
	
    float h = R*R - Q*Q*Q;
    float z = 0.0;
    if( h < 0.0 )
    {
    	// 4 intersections
        float sQ = sqrt(Q);
        z = 2.0*sQ*cos( acos(R/(sQ*Q)) / 3.0 );
    }
    else
    {
        // 2 intersections
        float sQ = pow( sqrt(h) + abs(R), 1.0/3.0 );
        z = sign(R)*abs( sQ + Q/sQ );
    }		
    z = c2 - z;
	
    float d1 = z   - 3.0*c2;
    float d2 = z*z - 3.0*c0;
    if( abs(d1) < 1.0e-4 )
    {
        if( d2 < 0.0 ) return -1.0;
        d2 = sqrt(d2);
    }
    else
    {
        if( d1 < 0.0 ) return -1.0;
        d1 = sqrt( d1/2.0 );
        d2 = c1/d1;
    }

    //----------------------------------
	
    float result = 1e20;

    h = d1*d1 - z + d2;
    if( h > 0.0 )
    {
        h = sqrt(h);
        float t1 = -d1 - h - k3; t1 = (po<0.0)?2.0/t1:t1;
        float t2 = -d1 + h - k3; t2 = (po<0.0)?2.0/t2:t2;
        if( t1 > 0.0 ) result=t1; 
        if( t2 > 0.0 ) result=min(result,t2);
    }

    h = d1*d1 - z - d2;
    if( h > 0.0 )
    {
        h = sqrt(h);
        float t1 = d1 - h - k3;  t1 = (po<0.0)?2.0/t1:t1;
        float t2 = d1 + h - k3;  t2 = (po<0.0)?2.0/t2:t2;
        if( t1 > 0.0 ) result=min(result,t1);
        if( t2 > 0.0 ) result=min(result,t2);
    }

    return result;
}

// df(x)/dx
//analitic function https://www.shadertoy.com/view/4sBGDy
vec3 nTorus( in vec3 pos, vec2 tor )
{
	return normalize( pos*(dot(pos,pos)- tor.y*tor.y - tor.x*tor.x*vec3(1.0,1.0,-1.0)));
}


HIT giper3D(vec3 ro, vec3 rd, vec2 torus)
{
    
    float dist = iTorus( ro, rd, torus );
    if( dist>0.0 )
    {
        vec3 pos = ro + dist*rd;
        vec3 nor = nTorus( pos, torus );
        return HIT(dist, nor, pos);
    }
    else
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
    
    
    vec3 light = normalize(vec3(0.0, 0.0, -1.0)); //light

		
    //
    
    float t = iTime/2.;
    vec2 m = vec2(0.0, 0.0);
    //if  (iMouse.z > 0.0)
    //m = (-iResolution.xy + 2.0*(iMouse.xy))/iResolution.y;
    vec3 ro = vec3(0.0, 0.0, 1.5); // camera
    ro = rotateX(PI/2.0)*ro; //camera rotation
    light = rotateX(PI/2.0)*light;
    ro = rotateY(-m.x*TAU)*rotateX(-m.y*PI)*ro; //camera rotation

    
    
    
    const float fl = 1.5; // focal length
    float dist = dist_infin;
    mat3 rota  = rotateZ(t/2.0)*rotateY(0.);
    mat3 rota_1  = rotateY(0.)*rotateZ(-t/2.0);
    //mat3 sky = rotateZ(0.0)*rotateX(PI/2.0);
    
    vec2 torus = vec2(1.0,0.8);
    vec3 tot = vec3(0.0);
    
    #define AA 2
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

        HIT giper = giper3D(rota*ro, rota*rd, torus);
        if (giper.dist < dist)
        {
            //col = vec3(1.0, 1.0, 0.0);
            vec3 backcol = vec3(0.0, 1.0, 1.0);
            vec3 nor = rota_1*giper.nor;
            float x = aafi(vec2(giper.pos.x, giper.pos.y)) /TAU;
            //x = fract(x);
            float y = aafi(vec2(torus.x - length(giper.pos.xy), giper.pos.z)) /TAU;
            y = 1.0 - y;
            col = texture(iChannel0, vec2(x,y)).rgb;
            //col = culccolor(col,col, -rd, light, nor);
            // gamma
            //col = pow( col, vec3(0.4545) ); 
            //reflect
            //col = calcSkyReflect(-rd, nor, sky);
        }
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