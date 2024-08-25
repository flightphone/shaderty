//https://www.johndcook.com/blog/2017/02/16/simulating-seashells/

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
#define TAU 6.28318530718

// The MIT License
// Copyright © 2018 Inigo Quilez
// Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions: The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software. THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

// Distance to a quadratic bezier segment, which can be solved analyically with a cubic.

// List of some other 2D distances: https://www.shadertoy.com/playlist/MXdSRf
//
// and iquilezles.org/articles/distfunctions2d


// 0: exact, using a cubic colver
// 1: approximated
//#define METHOD 1

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


float dot2( in vec2 v ) { return dot(v,v); }
float cro( in vec2 a, in vec2 b ) { return a.x*b.y - a.y*b.x; }

//#if METHOD==0
// signed distance to a quadratic bezier
float sdBezier( in vec2 pos, in vec2 A, in vec2 B, in vec2 C )
{    
    vec2 a = B - A;
    vec2 b = A - 2.0*B + C;
    vec2 c = a * 2.0;
    vec2 d = A - pos;

    float kk = 1.0/dot(b,b);
    float kx = kk * dot(a,b);
    float ky = kk * (2.0*dot(a,a)+dot(d,b))/3.0;
    float kz = kk * dot(d,a);      

    float res = 0.0;
    float sgn = 0.0;

    float p  = ky - kx*kx;
    float q  = kx*(2.0*kx*kx - 3.0*ky) + kz;
    float p3 = p*p*p;
    float q2 = q*q;
    float h  = q2 + 4.0*p3;

    if( h>=0.0 ) 
    {   // 1 root
        h = sqrt(h);
        vec2 x = (vec2(h,-h)-q)/2.0;

        #if 0
        // When p≈0 and p<0, h-q has catastrophic cancelation. So, we do
        // h=√(q²+4p³)=q·√(1+4p³/q²)=q·√(1+w) instead. Now we approximate
        // √ by a linear Taylor expansion into h≈q(1+½w) so that the q's
        // cancel each other in h-q. Expanding and simplifying further we
        // get x=vec2(p³/q,-p³/q-q). And using a second degree Taylor
        // expansion instead: x=vec2(k,-k-q) with k=(1-p³/q²)·p³/q
        if( abs(p)<0.001 )
        {
            float k = p3/q;              // linear approx
          //float k = (1.0-p3/q2)*p3/q;  // quadratic approx 
            x = vec2(k,-k-q);  
        }
        #endif

        vec2 uv = sign(x)*pow(abs(x), vec2(1.0/3.0));
        float t = clamp( uv.x+uv.y-kx, 0.0, 1.0 );
        vec2  q = d+(c+b*t)*t;
        res = dot2(q);
    	sgn = cro(c+2.0*b*t,q);
    }
    else 
    {   // 3 roots
        float z = sqrt(-p);
        float v = acos(q/(p*z*2.0))/3.0;
        float m = cos(v);
        float n = sin(v)*1.732050808;
        vec3  t = clamp( vec3(m+m,-n-m,n-m)*z-kx, 0.0, 1.0 );
        vec2  qx=d+(c+b*t.x)*t.x; float dx=dot2(qx), sx = cro(c+2.0*b*t.x,qx);
        vec2  qy=d+(c+b*t.y)*t.y; float dy=dot2(qy), sy = cro(c+2.0*b*t.y,qy);
        if( dx<dy ) { res=dx; sgn=sx; } else {res=dy; sgn=sy; }
    }
    
    return sqrt( res )*sign(sgn);
}
//#else

// This method provides just an approximation, and is only usable in
// the very close neighborhood of the curve. Taken and adapted from
// http://research.microsoft.com/en-us/um/people/hoppe/ravg.pdf
//#endif

float udSegment( in vec2 p, in vec2 a, in vec2 b )
{
	vec2 pa = p - a;
	vec2 ba = b - a;
	float h = clamp( dot(pa,ba)/dot(ba,ba), 0.0, 1.0 );
	return length( pa - ba*h );
}

float rand(float t, float n)
{
    float res =  fract(sin(t)*n);
    return res;
}
float getfi(float t, float n)
{
    if (t == 0.0)   
        return PI/3.0;
    float r = rand(t, n);
    return  r*PI/3.0 + PI/20.0;
}

float yBezier(vec2 v0, vec2 v1, vec2 v2, float x)
{
    //t^2(P0-2P1+P2) + t(2P1-2P0) + P0
    vec2 a2 = v0 - 2.0*v1 + v2;
    vec2 a1 = (2.0*v1 - 2.0*v0);
    vec2 a0 = v0;
    /*
    vec2 res = vec2(0.0);
    int nroot = quadratic(a2.x, a1.x, a0.x - x, res);
    float t = 0.;
    for (int i = 0; i < 2; i++)
    {
        if (res[i] >= 0.0 && res[i] <= 1.0)
        {
            t = res[i];
            break;
        }
    }
    */
    float t = x;
    float yy = t*t*a2.y + t*a1.y + a0.y;
    return yy;
}

float tileCurve(float p_x, float n)
{
    float t0 = floor(p_x);
    float t1 = t0 + 1.0;
    float px = fract(p_x);
    float a = getfi(t1, n) ;
    float b = getfi(t0, n) ;
    float x = tan(b)/(tan(a) + tan(b));
    float h = x*tan(a);
    x = 1.0 - x;
    //x = 0.5;
    //h = 0.2;
    float y = h;
    if (mod(t1, 2.0) == 0.0)
        y = - h;

    vec2 v0 = vec2(0.,0.);
    vec2 v1 = vec2(x,y);
    vec2 v2 = vec2(1.0,0.);

    float val = yBezier(v0, v1, v2, px);
    return val;
}
void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
	vec2 p = (2.0*fragCoord-iResolution.xy)/iResolution.y;
    float nc = 5.0;
    //vec2 p = fragCoord/iResolution.xy;
    
    //float f = smoothstep(-0.2,0.2,cos(2.0*iTime));
    //vec3 col = vec3(1.0) - vec3(0.1,0.4,0.7)*mix(sign(d),1.0,f);
    vec3 col = vec3(0.0);
    float z = tileCurve(p.x*nc, 2873.56);
    float psd = smoothstep(0.05,0.0,abs(z-p.y) );
	col = mix( col, vec3(1.0), psd);
	fragColor = vec4(col,1.0);
}

/////=====================================================================================
void main() {
    vec4 fragColor = vec4(0);
    mainImage(fragColor, gl_FragCoord.xy);
    gl_FragColor = fragColor;
}