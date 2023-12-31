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
    vec3 R = reflect (light1, nor);
    
    float difu1 = dot(nor, light1);
    float difu2 = dot(nor, light2);
    float difu = max(difu1, difu2);
    float shininess=10.0;
    float h = max(dot(R, rd), 0.);
    float specular    =  pow(max(dot(R, rd), 0.), shininess);

    //float amb = 0.5;// + 0.5 * max(dot(nor, light1), 0.);
    col = col*clamp(difu, 0.3, 1.0) + 0.5*specular*specular;
    col = clamp(col, vec3(0.), vec3(1.0));
    return col;   
}

//===================https://www.shadertoy.com/view/wsXGWS======================
float sgn(float x) {
  return x < 0.0? -1.0: 1.0; // Return 1 for x == 0
}

int quadratic(float A, float B, float C, out vec2 res) {
  float x1,x2;
  float b = -0.5*B;
  float q = b*b - A*C;
  if (q < 0.0) return 0;
  float r = b + sgn(b)*sqrt(q);
  if (r == 0.0) {
    x1 = C/A; x2 = -x1;
  } else {
    x1 = C/r; x2 = r/A;
  }
  res = vec2(x1,x2);
  return 2;
}

int quadratic(vec3 coeffs, out vec2 res) {
  return quadratic(coeffs[0],coeffs[1],coeffs[2],res);
}

void eval(float X, float A, float B, float C, float D,
          out float Q, out float Q1, out float B1,out float C2) {
  float q0 = A*X;
  B1 = q0+B;
  C2 = B1*X+C;
  Q1 = (q0+B1)*X + C2;
  Q = C2*X + D;
}

// Solve: Ax^3 + Bx^2 + Cx + D == 0
// Find one real root, then reduce to quadratic.
int cubic(float A, float B, float C, float D, out vec3 res) {
  float X,b1,c2;
  if (A == 0.0) {
    X = 1e8; A = B; b1 = C; c2 = D;
  } else if (D == 0.0) {
    X = 0.0; b1 = B; c2 = C;
  } else {
    X = -(B/A)/3.0;
    float t,r,s,q,dq,x0;
    eval(X,A,B,C,D,q,dq,b1,c2);
    t = q/A; r = pow(abs(t),1.0/3.0); s = sgn(t);
    t = -dq/A; if (t > 0.0) r = 1.324718*max(r,sqrt(t));
    x0 = X - s*r;
    if (x0 != X) {
      for (int i = 0; i < 6; i++) {
        X = x0;
        eval(X,A,B,C,D,q,dq,b1,c2);
        if (dq == 0.0) break;
        x0 -= (q/dq);
      }
      if (abs(A)*X*X > abs(D/X)) {
        c2 = -D/X; b1 = (c2 - C)/X;
      }
    }
  }
  res.x = X;
  return 1 + quadratic(A,b1,c2,res.yz);
}

int cubic(vec4 coeffs, out vec3 res) {
  float A = coeffs[0], B = coeffs[1], C = coeffs[2], D = coeffs[3];
  return cubic(A,B,C,D,res);
}

//https://www.shadertoy.com/view/wsXGWS



HIT giper3D(vec3 ro, vec3 rd, float t, float r)
{
    float a = ro.x;
    float b = rd.x;
    float c = ro.y;
    float d = rd.y;
    float e = ro.z;
    float f = rd.z;

    float k = 2.;
    //https://github.com/flightphone/shaderty/blob/master/staples_polynomial.py
    //for generate this expression used python script staples_polynomial.py
    float a0 = 1.*a*a*c + 3.*a*c*e + 1.*a*a*e + 1.*a*c*c + 1.*c*c*e + 1.*c*e*e + 1.*a*e*e-1.*a*c*t-1.*c*e*t-1.*a*e*t-1.*a*c*e*k;
    float a1 = 1.*a*a*d + 2.*a*b*c + 3.*a*c*f + 3.*a*d*e + 2.*a*b*e + 1.*a*a*f + 2.*a*c*d + 1.*b*c*c + 1.*c*c*f + 2.*c*d*e + 3.*b*c*e + 2.*c*e*f + 1.*d*e*e + 1.*b*e*e + 2.*a*e*f-1.*a*d*t-1.*b*c*t-1.*c*f*t-1.*d*e*t-1.*b*e*t-1.*a*f*t-1.*a*c*f*k-1.*a*d*e*k-1.*b*c*e*k;
    float a2 = 2.*a*b*d + 3.*a*d*f + 2.*a*b*f + 2.*b*c*d + 2.*c*d*f + 3.*b*c*f + 3.*b*d*e + 2.*d*e*f + 2.*b*e*f-1.*b*d*t-1.*d*f*t-1.*b*f*t + 1.*b*b*c + 1.*b*b*e + 1.*a*d*d + 1.*d*d*e + 1.*c*f*f + 1.*a*f*f-1.*a*d*f*k-1.*b*c*f*k-1.*b*d*e*k;
    float a3 = 1.*b*b*d + 3.*b*d*f + 1.*b*b*f + 1.*b*d*d + 1.*d*d*f + 1.*d*f*f + 1.*b*f*f-1.*b*d*f*k;
    //https://github.com/flightphone/shaderty/blob/master/staples_polynomial.py

    vec3 roots = vec3(dist_infin);
    int nroots = cubic(a3, a2, a1, a0, roots);
    
    
    float dist = dist_infin;
    vec3 pos = vec3(0.0);
    vec3 nor = vec3(0.0);
    for (int i = 0; i < 3; i++)
    {
        if (i >= nroots)
            break;
        if (roots[i] < 0.0)
            continue;
        vec3 p = ro + roots[i]*rd;
        if (length(p) > r)    
            continue;
        if (roots[i] < dist)    
        {
            dist = roots[i];
            pos = p;
        }

    }
    if (dist < dist_infin)
    {
        nor = vec3(0.+2.*pos.x*pos.y+3.*pos.y*pos.z+2.*pos.x*pos.z+1.*pos.y*pos.y+1.*pos.z*pos.z-1.*t*pos.y-1.*t*pos.z-1.*k*pos.y*pos.z, 0.+1.*pos.x*pos.x+3.*pos.x*pos.z+2.*pos.x*pos.y+2.*pos.y*pos.z+1.*pos.z*pos.z-1.*t*pos.x-1.*t*pos.z-1.*k*pos.x*pos.z, 0.+3.*pos.x*pos.y+1.*pos.x*pos.x+1.*pos.y*pos.y+2.*pos.y*pos.z+2.*pos.x*pos.z-1.*t*pos.y-1.*t*pos.x-1.*k*pos.x*pos.y);
        nor = normalize(nor);
    }
    return HIT(dist, nor, pos);
    
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

//https://www.shadertoy.com/view/dlKBRc
vec3 lightingv3(vec3 lightColor, vec3 backcol, vec3 rd, vec3 L, vec3 normal) 
{   
    float d = dot(rd, normal);
    if (d < 0.0)
    {
        lightColor = backcol;
        normal = -normal;
    }

    vec3 V = rd;
    vec3 N = normal;
    vec3 R = reflect (-L, N);
    float shadow = 1.;
    float occ = 1.;
    float Ka = 0.5;
    //vec3 ambient = 0.0;//Ka + Ka * lightColor; //dot(normal, vec3(0., 1., 0.))
    //ambient*=0.5;
    float diff= clamp(dot(N, L), 0.3, 1.0);
    vec3 diffuse =  lightColor * diff;
    float shininess=10.0;
    float specular    =  pow(max(dot(R, V), 0.), shininess);
    vec3 back = 0.5 * lightColor * clamp(dot(N, -L), 0.0, 1.0); // back
    vec3 colOut = occ*lightColor*(diffuse*shadow+.25 +back) + vec3(.5)*specular*specular;
    return clamp(colOut, vec3(0.0), vec3(1.0));
}


//https://www.shadertoy.com/view/Xds3zN , 555 row
vec3 iqcolor(vec3 col, vec3 backcol, vec3 rd, vec3 light, vec3 light2, vec3 nor)
{
    if (dot(rd, nor) < 0.)
    {
        col = backcol;
        nor = -nor;
    }
    float ks = 1.;
    // lighting
    float occ = .8;//calcAO( pos, nor );
    vec3 lin = vec3(0.0);
    vec3 ref = reflect( rd, nor );


        // sun
        {
            vec3  lig = light;
            vec3  hal = normalize( lig-rd );
            float dif = clamp( dot( nor, lig ), 0.3, 1.0 );
            //if( dif>0.0001 )
        	//dif *= calcSoftshadow( pos, lig, 0.02, 2.5 );
			float spe = pow( clamp( dot( nor, hal ), 0.0, 1.0 ),16.0);
                  spe *= dif;
                  spe *= 0.04+0.96*pow(clamp(1.0-dot(hal,lig),0.0,1.0),5.0);
                //spe *= 0.04+0.96*pow(clamp(1.0-sqrt(0.5*(1.0-dot(rd,lig))),0.0,1.0),5.0);
            lin += col*2.20*dif*vec3(1.30,1.00,0.70);
            lin +=     5.00*spe*vec3(1.30,1.00,0.70)*ks;
        }

        // sun2
        {
            vec3  lig = light2;
            vec3  hal = normalize( lig-rd );
            float dif = clamp( dot( nor, lig ), 0.3, 1.0 );
            //if( dif>0.0001 )
        	//dif *= calcSoftshadow( pos, lig, 0.02, 2.5 );
			float spe = pow( clamp( dot( nor, hal ), 0.0, 1.0 ),16.0);
                  spe *= dif;
                  spe *= 0.04+0.96*pow(clamp(1.0-dot(hal,lig),0.0,1.0),5.0);
                //spe *= 0.04+0.96*pow(clamp(1.0-sqrt(0.5*(1.0-dot(rd,lig))),0.0,1.0),5.0);
            lin += col*2.20*dif*vec3(1.30,1.00,0.70);
            lin +=     5.00*spe*vec3(1.30,1.00,0.70)*ks;
        }

        lin /= 2.;


        
        // sky
        
        {
            float dif = sqrt(clamp( 0.5+0.5*nor.z, 0.0, 1.0 ));
                  dif *= occ;
            float spe = smoothstep( -0.2, 0.2, ref.z );
                  spe *= dif;
                  spe *= 0.04+0.96*pow(clamp(1.0+dot(nor,rd),0.0,1.0), 5.0 );
          //if( spe>0.001 )
                  //spe *= calcSoftshadow( pos, ref, 0.02, 2.5 );
            lin += col*0.60*dif*vec3(0.40,0.60,1.15);
            lin +=     2.00*spe*vec3(0.40,0.60,1.30)*ks;
        }
        

        // back
        /*
        {
        	float dif = clamp( dot( nor, normalize(vec3(0.5,0.0,0.6))), 0.0, 1.0 )*clamp( 1.0-pos.y,0.0,1.0);
                  dif *= occ;
        	lin += col*0.55*dif*vec3(0.25,0.25,0.25);
        }
        */
        // sss
        {
            float dif = pow(clamp(1.0+dot(nor,rd),0.0,1.0),2.0);
                  dif *= occ;
        	lin += col*0.25*dif*vec3(1.00,1.00,1.00);
        }
        return lin;
}
void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    //surface (x+y+z-a)(xy+yz+zx) - kxyz = 0
    vec3 light = normalize(vec3(.0, 1.0, 1.0)); //light
    vec3 light2 = normalize(vec3(.0, -1.0, -1.0)); //light

    float ra = 3.0;
    float g = 1.0;

    float t = iTime/2.0;
    vec2 m = vec2(0.0, 0.0);
    //if  (iMouse.z > 0.0)
    {
        m = (-iResolution.xy + 2.0*(iMouse.xy))/iResolution.y;
        t = 0.;
    }
    vec3 ro = vec3(0.0, 0.0, 6.); // camera
    ro = rotateY(-m.x*TAU)*rotateX(-m.y*PI)*ro; //camera rotation
    
    
    const float fl = 1.5; // focal length
    float dist = dist_infin;
    float fi = PI/4.5;
    
    mat3 rota  = rotateX(fi)*rotateY(-fi)*rotateX(-PI/2.)*rotateY(t);
    mat3 rota_1  = rotateY(-t)*rotateX(PI/2.)*rotateY(fi)*rotateX(-fi);
    //mat3 sky = rotateZ(0.0)*rotateY(PI/2.0);
    
    vec2 torus = vec2(1.0,0.3);
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
        HIT giper = giper3D(rota*ro, rota*rd, g, ra);
        if (giper.dist < dist)
        {
            col = vec3(0.1, 0.1, .9);
            vec3 backcol = vec3(1.0, 0.2, 0.2);
            vec3 nor = rota_1*giper.nor;
            col = culccolor(col, backcol, -rd, light, light2, nor);
            //col = iqcolor(col, backcol, rd, light, light2, nor);
            //col = lightingv3(col, backcol, -rd, light, nor);
            // gamma
            //col = pow( col, vec3(0.7) ); 
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