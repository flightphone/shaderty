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

float aafi(vec2 p) {
    float fi = atan(p.y, p.x);
    fi += step(p.y, 0.0)*TAU;
    return fi;
}


//converts a vector on a sphere to longitude and latitude
vec2 lonlat (vec3 p)
{
    float lon = aafi(p.xy)/TAU;
    float lat = aafi(vec2(p.z, length(p.xy)))/PI;
    return vec2(lon, lat);
}

const float dist_infin = 6.0;
const HIT hit_inf = HIT(dist_infin, vec3(0.0), vec3(0.0));
#define nn 128
const float eps = 0.001;

vec3 culccolor(vec3 col_in, vec3 backcol, vec3 rd, vec3 light, vec3 nor)
{
    vec3 col = col_in;
    float d = dot(rd, nor);
    if (d < 0.0)
        col = backcol;
    
    nor *= -sign(d);
    float difu = dot(nor, light);
    //float difu2 = dot(nor, light2);
    //float difu = max(difu1, difu2);
        col *= clamp(difu, 0.3, 1.0);
    return col;   
}

//https://iquilezles.org/articles/smin/
// polynomial smooth min 1 (k=0.1)
float smin( float a, float b, float k )
{
    float h = clamp( 0.5+0.5*(b-a)/k, 0.0, 1.0 );
    return mix( b, a, h ) - k*h*(1.0-h);
}

float dot2( in vec3 v ) { return dot(v,v); }

float sdSphere( vec3 p, float s )
{
  return length(p)-s;
}

float sdBoxFrame( vec3 p, vec3 b, float e)
{
       p = abs(p  )-b;
  vec3 q = abs(p+e)-e;

  return sqrt(min(min(dot2(max(vec3(p.x,q.y,q.z),0.0)),
                      dot2(max(vec3(q.x,p.y,q.z),0.0))),
                      dot2(max(vec3(q.x,q.y,p.z),0.0)))) 
         +min(0.0,min(min( max(p.x,max(q.y,q.z)),
                           max(p.y,max(q.z,q.x))),
                           max(p.z,max(q.x,q.y))));
}
float sdRoundBox( vec3 p, vec3 b, float r )
{
  vec3 q = abs(p) - b;
  return length(max(q,0.0)) + min(max(q.x,max(q.y,q.z)),0.0) - r;
}

float sdTorus( vec3 p, vec2 t )
{
  vec2 q = vec2(length(p.xz)-t.x,p.y);
  return length(q)-t.y;
}

float sdHexagram( in vec2 p, in float r )
{
    const vec4 k = vec4(-0.5,0.86602540378,0.57735026919,1.73205080757);
    
    p = abs(p);
    p -= 2.0*min(dot(k.xy,p),0.0)*k.xy;
    p -= 2.0*min(dot(k.yx,p),0.0)*k.yx;
    p -= vec2(clamp(p.x,r*k.z,r*k.w),r);
    return length(p)*sign(p.y);
}

float sdHexagram3( in vec3 p, in float h, in float r )
{
    float d = sdHexagram(p.xy, r);
    vec2 w = vec2( d, abs(p.z) - h );
    return min(max(w.x,w.y),0.0) + length(max(w,0.0)) - 0.05;
}
float sdLine(in vec2 p, vec2 ro, vec2 rd)
{
    vec2 sec = p-ro;
    float d = dot(sec, rd)/dot(rd, rd);
    return length(sec - rd*d);
}

float sdSegment(in vec3 p, vec3 a, vec3 b)
{
    vec3 sec = p-a;
    vec3 ba = b-a;
    float d = clamp(dot(sec, ba)/dot(ba, ba), 0.0, 1.0);
    return length (sec - ba*d);
}

float sdSegmentR(in vec3 p, vec3 a, vec3 b, float r)
{
    float h = (p.x - a.x)/(b.x-a.x);
    float dx = (r + 0.5*r*cos(30.*h));
    float d = sdSegment(p, a, b);
    //if (d < 2.0 * dx); 
        d = d-dx;
    return d;
}


float desp(vec3 p)
{
    //return  0.01 * sin(20.*p.x)*sin(20.*p.y)*sin(20.*p.z);
    //return sdTorus(p, vec2(0.1, 0.6));
    //return 0.;
    vec2 fon = lonlat(p); //get longitude and latitude
    float f3 = texture(iChannel0, fon).w;
    f3 = 0.1*pow(f3, 1.0);
    return  f3;
}

float desp2(vec3 p)
{
    vec2 fon = lonlat(p); //get longitude and latitude
    if (fon.y < 0.05 || fon.y > 0.95)        
        return 0.;
    float k = 150.0;
    float f3 = sin(k*fon.x) + sin(k*fon.y);
    return  f3*0.02;
}

float k = 15.;

float ff (float x)
{
    return 0.1*x + 0.05*sin(k*x + iTime) + 0.2*cos(k*x - iTime);
    //return 0.5*x*x;
    
    //return sqrt(1.5 - x*x);
}

vec3 ffp (float x)
{

    const float h = 0.0001;
    vec2 p = normalize(vec2(2.*h, ff(x+h)-ff(x-h)));
    return vec3(p.x, p.y, 0.);
}

float sdFunP(in vec3 p, float a, float b)
{
    float x = p.x; //clamp(p.x, a, b);
    vec3 val = vec3(x, ff(x), 0.0);
    vec3 sec = vec3(p.x - val.x, p.y - val.y, 0.0);
    float l = length(sec);
    
    float d = 0.5*length(cross(sec, ffp(x)));
    d = sqrt(p.z*p.z * d * d / l / l   + d*d);
    return d;
}



float sdFun(in vec3 p, float a, float b)
{
    float x = clamp(p.x, a, b);
    vec2 val = vec2(x, ff(x));
    float d = length(p.xy - val)/4.;
    d = sqrt(p.z*p.z * 0.25 * 0.25 + d*d);
    return d;
}

float sdFun2(in vec2 p, float a, float b)
{
    float x = clamp(p.x, a, b);
    vec2 val = vec2(x, ff(x));
    float d = length(p.xy - val)/4.;
    return d;
}

float sdFun3(in vec3 p, float a, float b)
{
    float h = 0.4;// + 0.2 * sin(3.*p.x);
    float r = 0.005;
    float d = sdFun2(p.xy, a, b);
    vec2 w = vec2( d, abs(p.z) - h );
    return min(max(w.x,w.y),0.0) + length(max(w,0.0)) - r;
}

float ff3(float x, float y)
{
    float k = 15.;
    return 0.1*(sin(k*x) + sin(k*y));
    
}

float sdEggBox(in vec3 p, float a, float b)
{
    float x = clamp(p.x, a, b);
    float y = clamp(p.y, a, b);
    vec3 val = vec3(x, y, ff3(x, y));
    float d = length(p - val)/4.;
    return d;
}

float sdTexture(in vec3 p, float a, float b)
{
    float x = clamp(p.x, a, b);
    float y = clamp(p.y, a, b);

    float x1 = (x-a)/(b-a);
    float y1 = (y-a)/(b-a);

    float f3 = texture(iChannel0, vec2(x1,y1)).w;
    f3 = pow(f3, 0.1);
    

    vec3 val = vec3(x, y, f3);
    float d = length(p - val)/4.;
    return d;
}    

float udTriangle( vec3 p, vec3 a, vec3 b, vec3 c )
{
  vec3 ba = b - a; vec3 pa = p - a;
  vec3 cb = c - b; vec3 pb = p - b;
  vec3 ac = a - c; vec3 pc = p - c;
  vec3 nor = cross( ba, ac );

  return sqrt(
    (sign(dot(cross(ba,nor),pa)) +
     sign(dot(cross(cb,nor),pb)) +
     sign(dot(cross(ac,nor),pc))<2.0)
     ?
     min( min(
     dot2(ba*clamp(dot(ba,pa)/dot2(ba),0.0,1.0)-pa),
     dot2(cb*clamp(dot(cb,pb)/dot2(cb),0.0,1.0)-pb) ),
     dot2(ac*clamp(dot(ac,pc)/dot2(ac),0.0,1.0)-pc) )
     :
     dot(nor,pa)*dot(nor,pa)/dot2(nor) );
}

float udQuad( vec3 p, vec3 a, vec3 b, vec3 c, vec3 d )
{
  vec3 ba = b - a; vec3 pa = p - a;
  vec3 cb = c - b; vec3 pb = p - b;
  vec3 dc = d - c; vec3 pc = p - c;
  vec3 ad = a - d; vec3 pd = p - d;
  vec3 nor = cross( ba, ad );

  return sqrt(
    (sign(dot(cross(ba,nor),pa)) +
     sign(dot(cross(cb,nor),pb)) +
     sign(dot(cross(dc,nor),pc)) +
     sign(dot(cross(ad,nor),pd))<3.0)
     ?
     min( min( min(
     dot2(ba*clamp(dot(ba,pa)/dot2(ba),0.0,1.0)-pa),
     dot2(cb*clamp(dot(cb,pb)/dot2(cb),0.0,1.0)-pb) ),
     dot2(dc*clamp(dot(dc,pc)/dot2(dc),0.0,1.0)-pc) ),
     dot2(ad*clamp(dot(ad,pd)/dot2(ad),0.0,1.0)-pd) )
     :
     dot(nor,pa)*dot(nor,pa)/dot2(nor) );
}

vec3 point(vec2 ll, float r)
{
    float f1 = ll.x * TAU;
    float f2 = ll.y * PI;
    float z = r*cos(f2);
    float d = abs(r*sin(f2));
    float x = d*cos(f1);
    float y = d*sin(f1);
    return vec3(x, y, z);
}

float sdDiscoBall(vec3 pos, float r)
{
    vec2 ll = lonlat(pos);
    float n = 15.;
    float n2 = 30.;
    ll.x = floor(ll.x*n2);
    ll.y = floor(ll.y*n);
    vec3 a = point(vec2(ll.x/n2, ll.y/n), r);
    vec3 b = point(vec2(ll.x/n2, (ll.y+1.)/n), r);
    vec3 c = point(vec2((ll.x + 1.)/n2, (ll.y+1.)/n), r);
    float d = dot(normalize(cross(b-a, c-a)), (pos - a));
    return abs(d*0.9);
}

float map_old( in vec3 pos )
{
    //return sdBoxFrame(pos, vec3(0.5,0.3,0.5), 0.025 );
    //return sdRoundBox(pos, vec3(0.5, 0.2, 0.15), 0.03);
    //return sdTorus(pos, vec2(0.5, 0.2));
    float d1 = sdRoundBox(pos, vec3(0.5, 0.2, 0.15), 0.03);
    float d2 = sdTorus(pos, vec2(0.5, 0.2));
    //return min(d1,d2);
    return max(d1,-d2);
}

float map( in vec3 pos )
{
    //return sdBoxFrame(pos, vec3(0.5,0.3,0.5), 0.025 );
    //return sdRoundBox(pos, vec3(0.5, 0.2, 0.15), 0.03);
    //return sdHexagram3(pos, 0.2, 0.5);
    //return sdSegmentR(pos, vec3(-.8, 0., 0.), vec3(.8, 0., 0.), 0.1);
    //return sdSphere(pos, 1.3) + desp2(pos);
    //return sdSphere(pos, 1.3) + desp(pos);
    //return sdSegmentR(pos, vec3(-1., 0., 0.), vec3(1., 0., 0.), 0.2);
    //return sdFunP(pos, -1.5, 1.5) - 0.01;//  + desp(pos);
    //return sdFun(pos, -1.5, 1.5) - 0.01;
    //return sdFun3(pos, -1., 1.);
    //return sdEggBox(pos, -1., 1.) - 0.001;
    //return sdTexture(pos, -1.5, 1.5) - 0.001;
    return sdDiscoBall(pos, 1.0);
    
}

float fpolar(float fi)
{
    //float r = 1. * cos(fi) * cos(fi);
    float r = sin(2.0*fi);
    return r;
}

float ffpl(float fi)
{
    return -2.0 * cos(fi)*sin(fi);
}

vec3 fpolarp(float fi)
{
    float x = ffpl(fi)*cos(fi) - fpolar(fi)*sin(fi);
    float y = ffpl(fi)*sin(fi) + fpolar(fi)*cos(fi);
    vec3 res = normalize(vec3(x, y, 0.));
    return res;
}

float sdPolar(vec3 p)
{
    float fi = aafi(p.xy);
    float r  = fpolar(fi);
    float e = .3;

    float d = abs(length(p.xy) - r) * e;
    d = sqrt(p.z*p.z * e * e + d*d);
    d -= 0.02;
    return d;
}

float sdArc(vec3 p)
{
    float k = 0.1;
    float e = 0.5;
    float L = length(p.xy);
    float fi = aafi(p.xy);
    float d = dist_infin;
    for (float i = 0. ; i < 6.; i++)
    {
        float r = abs(L - k*(fi + TAU*i)) * 0.8;
        d = min(d, r);
    }
    d = sqrt(p.z*p.z * e * e  + d*d);
    d -= 0.05;
    return d;
}

float map3( in vec3 pos )
{

/*
    float d =  sdPolar(pos);
    float d2 = sdSphere(pos, 0.3);
    return smin(d, d2, 0.1);
*/    


    float d =  sdArc(pos);
    float d2 = sdSphere(pos, 0.2);
    //return smin(d, d2, 0.1);
    return min(d, d2);

//return sdTexture(pos, -1.5, 1.5) - 0.001;    
}
vec3 calcNormal( in vec3 pos )
{
    const float h = 0.0001; // replace by an appropriate value
    const vec2 k = vec2(1,-1);
    return normalize( k.xyy*map(pos + k.xyy*h ) + 
                      k.yyx*map(pos + k.yyx*h ) + 
                      k.yxy*map(pos + k.yxy*h ) + 
                      k.xxx*map(pos + k.xxx*h ) );
}


HIT giper3D(vec3 ro, vec3 rd)
{
    float t  = 0.;
    for (int i = 0; i < nn; i++)
    {
        vec3 pos = ro + rd*t;
        float h = map(pos);
        if (h < eps || t >= dist_infin)
            break;
        t += h;  
    }    

    if (t >= dist_infin)
        return hit_inf;
      
    vec3 pos = ro + t*rd;
    vec3 nor = calcNormal(pos);
    return HIT(t, nor, pos);
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
    vec3 light2 = normalize(vec3(0.0, 1.0, .0)); //light

    

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
    t = 0.;
    mat3 rota  = rotateZ(t)*rotateY(-t);
    mat3 rota_1  = rotateY(t)*rotateZ(-t);
    
    vec3 tot = vec3(0.0);
    
    #define AA 2
    //antiblick
    for( int m=0; m<AA; m++ )
    for( int n=0; n<AA; n++ )
    {
        vec2 o = vec2(float(m),float(n)) / float(AA) - 0.5;
        vec2 p = (-iResolution.xy + 2.0*(fragCoord+o))/iResolution.y;
        //vec2 p = (2.0*fragCoord-iResolution.xy)/iResolution.y;
        //vec3 rd = normalize( vec3(p,fl) ); // ray direction
        vec3 rd = GetRayDir(p, ro, vec3(0,0.,0), fl); //ray direction
        vec3 col = vec3(0.0);
        
        HIT giper = giper3D(rota*ro, rota*rd);
        if (giper.dist < dist)
        {
            vec3 nor = rota_1*giper.nor;
            float dif = clamp( dot(nor,light), 0.0, 1.0 );
            float amb = 0.5 + 0.5*dot(nor,light2);
            col = vec3(0.5,0.6,0.7)*amb + vec3(0.85,0.75,0.65)*dif;
        }
        // gamma        
        col = sqrt( col );
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