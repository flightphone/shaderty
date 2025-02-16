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

//noise

// Precision-adjusted variations of https://www.shadertoy.com/view/4djSRW

//float hash(float p) {    p = fract(p * 0.011);    p *= p + 7.5;    p *= p + p;    return fract(p);}

float hash(float n) {
    return fract(sin(n) * 437558.5453123);
}
//float hash(vec2 p) {vec3 p3 = fract(vec3(p.xyx) * 0.13); p3 += dot(p3, p3.yzx + 3.333); return fract((p3.x + p3.y) * p3.z); }
float hash(vec2 p) {
    return fract(1e4 * sin(17.0 * p.x + p.y * 0.1) * (0.1 + abs(sin(p.y * 13.0 + p.x))));
}  //!!!Best
//	<https://www.shadertoy.com/view/4dS3Wd>
//	By Morgan McGuire @morgan3d, http://graphicscodex.com

float noise(float x) {
    float i = floor(x);
    float f = fract(x);
    float u = f * f * (3.0 - 2.0 * f);
    return mix(hash(i), hash(i + 1.0), u);
}

float noise(vec2 x) {
    vec2 i = floor(x);
    vec2 f = fract(x);
	// Four corners in 2D of a tile
    float a = hash(i);
    float b = hash(i + vec2(1.0, 0.0));
    float c = hash(i + vec2(0.0, 1.0));
    float d = hash(i + vec2(1.0, 1.0));

	// Simple 2D lerp using smoothstep envelope between the values.
	// return mix(mix(a, b, smoothstep(0.0, 1.0, u.x)),
	//			mix(c, d, smoothstep(0.0, 1.0, u.x)),
	//			smoothstep(0.0, 1.0, u.y));

	// Same code, with the clamps in smoothstep and common subexpressions
	// optimized away.
    vec2 u = f * f * (3.0 - 2.0 * f);
    return mix(a, b, u.x) + (c - a) * u.y * (1.0 - u.x) + (d - b) * u.x * u.y;
}

// This one has non-ideal tiling properties that I'm still tuning
float noise(vec3 x) {
    const vec3 step = vec3(110, 241, 171);

    vec3 i = floor(x);
    vec3 f = fract(x);

	// For performance, compute the base input to a 1D hash from the integer part of the argument and the 
	// incremental change to the 1D based on the 3D -> 1D wrapping
    float n = dot(i, step);

    vec3 u = f * f * (3.0 - 2.0 * f);
    return mix(mix(mix(hash(n + dot(step, vec3(0, 0, 0))), hash(n + dot(step, vec3(1, 0, 0))), u.x), mix(hash(n + dot(step, vec3(0, 1, 0))), hash(n + dot(step, vec3(1, 1, 0))), u.x), u.y), mix(mix(hash(n + dot(step, vec3(0, 0, 1))), hash(n + dot(step, vec3(1, 0, 1))), u.x), mix(hash(n + dot(step, vec3(0, 1, 1))), hash(n + dot(step, vec3(1, 1, 1))), u.x), u.y), u.z);
}

float mapB(vec3 p) {
	//wave
	p += 15.;
	//p.x-= iTime*0.5;
	float r = noise(1.5*p.xy); //+ iTime * 0.2;
	r =  noise(r) + 0.5 * noise(4.*r) +  .3* noise(r * 10.) + 0.1* noise(r*20.);
	r = length(p) - r;
	r = r + 0.5*noise(r) + .2 * noise(2.*r);// +  0.1 * noise(r * 20.) + 0.1*noise(r*40.);
	return r;
	
}

vec3 calcNormalB(in vec3 p) {
	const float eps = 0.0001;
	vec2 q = vec2(0.0, eps);
	vec3 res = vec3(mapB(p + q.yxx) - mapB(p - q.yxx), mapB(p + q.xyx) - mapB(p - q.xyx), mapB(p + q.xxy) - mapB(p - q.xxy));
	return normalize(res);
}

vec3 wave(vec2 p) {
	vec3 pos = vec3(p, 3.);
	vec3 col = vec3(0.08, 0.42, 0.87);
	vec3 nor = calcNormalB(pos);
	vec3 light = normalize(vec3(0.0, 1.0, 2.5)); //light
	vec3 rd = normalize(vec3(1., 1., -1.));
	vec3 R = reflect(light, nor);
	float specular = pow(max(abs(dot(R, rd)), 0.), 16.);
	float difu = abs(dot(nor, light));
	col = col * (clamp(difu, 0., 1.0) + 0.5) + vec3(1., .7, .4) * specular;
	//float h = clamp(difu, 0., 1.0);
	//col = mix(col, pink, h*h);
	//col += vec3(1., .7, .4) * specular;
	float fre = pow(clamp(dot(nor, rd) + 1., .0, 1.), 3.); // Fresnel, for some mild glow.
	col += vec3(.1, .1, 0.1) * fre; //?
	col = sqrt(col);
	return col;
}

/////=====================================================================================

#define PI  3.14159265359
#define TAU 6.28318530718
#define rot(f) mat2(cos(f), -sin(f), sin(f), cos(f))

const float dist_infin = 10.0;
#define nn 128
const float eps = 0.001;

vec3 col1 = vec3(0.73, 0.59, 0.3);
vec3 col2 = vec3(0.8117, 0.1764, 0.8078);
vec3 resColor = vec3(0.73, 0.59, 0.3);

float ring(vec3 p, float h, float R, float r) {
    float d = abs(length(p.xy) - R) - r;
    vec2 w = vec2(d, abs(p.z) - h);
    return min(max(w.x, w.y), 0.0) + length(max(w, 0.0));
}

float quin(vec3 p) {
    float k = 1., x = p.x * k, y = p.y * k, z = p.z * k;
    float res = (x * x + y * y) - z * z * .1 - 0.5 + pow(min(z + 2., 0.), 2.) + pow(max(z - 4., 0.), 6.) - 2. * max(z - 4., 0.) + 1. * min(z + 2., 0.);
    z += 3.;
    float res2 = ring(vec3(x, y, z), 0.5, 1.2, 0.1) - 0.1;
    return min(res2, res);
}
float lines(vec3 p) {
    float a = 1., n = 5.;
    p.xy = clamp(p.xy, vec2(-a), vec2(a));
    float h = 0.75, r = 0.05;

    p.xy = (p.xy + a);
    p.xy /= ((2. * a) / n);
    p.xy = fract(p.xy) - 0.5;
    p.xy *= 2. * a / n;
    p.z -= clamp(p.z, -h, h);
    return length(p) - r;
}

vec3 colx = vec3(0.7, 0., 0.);
vec3 coly = vec3(0., 0.7, 0.);
vec3 colz = vec3(0., 0., 0.7);

float grid(vec3 p) {
    float t2 = lines(p);
    float t = t2;
    resColor = colx;

    float t3 = lines(p.yzx);
    if(t3 < t) {
        t = t3;
        resColor = coly;
    }
    float t4 = lines(p.xzy);
    if(t4 < t) {
        t = t4;
        resColor = colz;
    }
    return t;
}
float rose(vec3 p) {
    float lon = atan(p.y, p.x);
    float lat = atan(p.z, length(p.xy));
    float ll = sin(4. * lon) + cos(7. * lat) + iTime * 0.3;
    float r = 2. + 0.5 * noise(ll) + 0.1 * noise(ll * 5.);
    return length(p) - r;
}

float noiseSph(vec3 p) {
    float r = 1. + 0.7 * noise(p) + 0.3 * noise(p * 3.) + 0.05 * noise(p * 15.);
    return length(p) - r;
}

float cylinder(vec3 p, float R, float h) {

    float d = length(p.xy) - R;
    vec2 w = vec2(d, abs(p.z) - h);
    return min(max(w.x, w.y), 0.0) + length(max(w, 0.0));
}

float smin(float d1, float d2, float k) {
    float h = clamp(0.5 + 0.5 * (d2 - d1) / k, 0.0, 1.0);
    return mix(d2, d1, h) - k * h * (1.0 - h);
}

float smin_out(float d1, float d2, float k) {
    float h = clamp(0.5 - 0.5 * (d2 + d1) / k, 0.0, 1.0);
    return mix(d2, -d1, h) + k * h * (1.0 - h);
}

float cylinder2(vec3 p) {
    float t = cylinder(p, 0.7, 2.);
    t = smin(t, cylinder(p.yzx, 0.7, 2.), 0.3);
    t = smin(t, cylinder(p.xzy, 0.7, 2.), 0.3);   
    /*
    t = max(t, -cylinder(p, 0.3, 3.) );
    t = max(t, -cylinder(p.yzx, 0.3, 3.) );
    t = max(t, -cylinder(p.xzy, 0.3, 3.) );
    */
    t = smin_out(cylinder(p, 0.3, 3.), t, 0.05);
    t = smin_out(cylinder(p.yzx, 0.3, 3.), t, 0.05);
    t = smin_out(cylinder(p.xzy, 0.3, 3.), t, 0.05);

    return t - 0.05;
}

vec3 getSg(vec3 p, float nseg) {
    float fi = mod(atan(p.y, p.x), TAU);
    fi = mod(fi + PI / nseg, TAU);
    float n = floor(fi / TAU * nseg);
    p.xy *= rot(-n * TAU / nseg);
    return p;
}
float sdSegment(in vec2 p, in vec2 a, in vec2 b) {
    vec2 pa = p - a, ba = b - a;
    float h = clamp(dot(pa, ba) / dot(ba, ba), 0.0, 1.0);
    return length(pa - ba * h);
}

float umb(vec3 p) {
    p.yz*=rot(PI/2.);
    float n1 = 14., R = 2., h = 2.1;
    //============================handle=======================
    float d3 = length(vec2(length(vec2(p.xy)), max(abs(p.z) - h, 0.))) - 0.07;
    float d1 = length(vec2(length(vec2(p.xy)), 
    p.z - clamp(p.z, -h, -0.8*h))) - 0.12;
    
    p = getSg(p, n1);
    //============================spokes=====================
    float dz = 0.5, dr = sqrt(R*R - dz*dz);
    vec2 a0 = vec2(dr * cos(PI / n1), dr * sin(PI / n1)*sign(p.y));
    float d2 = sdSegment(p.xy, a0, vec2(0.));
    d2 = length(vec2(d2, p.z-dz)) - 0.02;
    //=======================================================

    //===================tent=========================
    vec3 p0 = normalize(p);
    vec2 a = vec2(R * cos(PI / n1), R * sin(PI / n1)), b = vec2(R * cos(PI / n1), -R * sin(PI / n1));
    float t = 0.;
    if(p.z > 0.) {
        t = sqrt(R * R / (p0.z * p0.z + p0.x * p0.x / cos(PI / n1) / cos(PI / n1)));
        float l = p0.x * t / cos(PI / n1);
        a = vec2(l * cos(PI / n1), l * sin(PI / n1));
        b = vec2(l * cos(PI / n1), -l * sin(PI / n1));
    }
    float d0 = sdSegment(p.xy, a, b);
    d0 = 0.95*length(vec2(d0, p.z - p0.z * t))-0.02;
    //=================================================

    //=====================edges========================
    p.xy *= rot(-PI / n1 * sign(p.y));
    float d = abs(length(p.xz) - R);
    if(p.z < 0.)
        d = length(p.xz - vec2(R, 0));
    d = length(vec2(d, p.y)) - 0.03;
    return min(min(min(min(d, d0), d3), d2), d1);

}

float sdBox( vec3 p, vec3 b )
{
  vec3 q = abs(p) - b;
  return length(max(q,0.0)) + min(max(q.x,max(q.y,q.z)),0.0);
}

float halloween(vec3 p)
{
    p.yz*=rot(PI/2.);
    float n1 = 14., R = 1.5;
    float d0 = length(p) - 1.45;

    float fi1 = TAU/n1+PI/n1, li = PI/4., r = R*cos(li);//PI/n1
    vec3 e1p = vec3(r*cos(fi1), r*sin(fi1), R*sin(li));
    float e1 = length(p - e1p) - .45;
    //fi1+= 3.*TAU/n1;
    fi1 = -fi1;
    e1p = vec3(r*cos(fi1), r*sin(fi1), R*sin(li));
    float e2 = length(p - e1p) - .45;

    fi1 = 0., li = -PI/10.;
    vec3 lip = vec3(r*cos(fi1), r*sin(fi1), R*sin(li));
    vec3 lv = p - lip;
    
    float lips = sdBox(lv, vec3(0.9, 1., 0.2));
    
    
    p = getSg(p, n1);
    p.xy *= rot(-PI / n1 * sign(p.y));
    float d = abs(length(p.xz) - R);
    d = length(vec2(d, p.y)) - 0.35;

    d = max(d, -d0);
    d = smin_out(e1, d, 0.1);
    d = smin_out(e2, d, 0.1);
    d = smin_out(lips, d, 0.2);
    //d = max(d, -lips);
    //d = max(d, -e2);
    return d;
}

float insole2 (vec3 p)
{
    float k = 0.5, x = p.x * k, y = p.y * k, z = p.z * k;
    float w = 1.8, r = 0.3, h = 0.01;
    x += w/2.;
    z += .5;
    float z0 = 1.0 - smoothstep(0.1, w * 0.75, x);
    z -= z0;
    float dx = clamp(x, 0., w);
    //r = r + dx * 0.12;
    r = r + 0.1*sin((w-dx)/w*TAU);

    float d = length(vec2(x - clamp(x, 0., w), y)) - r;
   
    if (x > w)
    {
        //d = length(vec2(x-w, max(abs(y)-r, 0.)));
        d = length(vec2(x-w+r, y)) - r*sqrt(2.)+0.007;
    }
    
    vec2 w2 = vec2(d, abs(z) - h);
    float res =  min(max(w2.x, w2.y), 0.0) + length(max(w2, 0.0));
    return res*0.95 - 0.02;
}

float heel2(vec3 p)
{
    float k = 0.45, x = p.x * k, y = p.y * k, z = p.z * k;
    float r = 0.05, w = 0.05, d = 0.;
    z += .5;
    x += .75;
    r = r * (2.2 * z * z + 1.);
    w = w * (2.2 * z * z + 1.);
    d = length(vec2(max(abs(x) - w, 0.), y)) - r;
    if (x > 0.)
    {
        d = length(vec2(x, max(abs(y)-r, 0.)));
    }
    vec2 w2 = vec2(d, abs(z) - clamp(z, 0.,  0.95));
    float res =  min(max(w2.x, w2.y), 0.0) + length(max(w2, 0.0));
    return res - 0.02;
}

float strap2(vec3 p)
{
    p.x += 1.;
    p.z -= 1.;
    float r = 0.7, h = 0.5;
    float d = abs(length(p.xy) - r);
    if (p.x < 0.)
    {
        d = length(vec2(p.x, max(abs(p.y)-r, 0.)));
    }
    vec2 w2 = vec2(d, abs(p.z) - clamp(p.z, 0.,  h));
    float res =  min(max(w2.x, w2.y), 0.0) + length(max(w2, 0.0));
    return res - 0.01;
}

float sh(vec3 p)
{
    float d = smin(heel2(p), insole2 (p), 0.01);
    d = smin(d, strap2(p.zyx), 0.01);
    return d;
}

float sdEllipse( in vec2 p, in vec2 ab )
{
    p = abs(p); if( p.x > p.y ) {p=p.yx;ab=ab.yx;}
    float l = ab.y*ab.y - ab.x*ab.x;
    float m = ab.x*p.x/l;      float m2 = m*m; 
    float n = ab.y*p.y/l;      float n2 = n*n; 
    float c = (m2+n2-1.0)/3.0; float c3 = c*c*c;
    float q = c3 + m2*n2*2.0;
    float d = c3 + m2*n2;
    float g = m + m*n2;
    float co;
    if( d<0.0 )
    {
        float h = acos(q/c3)/3.0;
        float s = cos(h);
        float t = sin(h)*sqrt(3.0);
        float rx = sqrt( -c*(s + t + 2.0) + m2 );
        float ry = sqrt( -c*(s - t + 2.0) + m2 );
        co = (ry+sign(l)*rx+abs(g)/(rx*ry)- m)/2.0;
    }
    else
    {
        float h = 2.0*m*n*sqrt( d );
        float s = sign(q+h)*pow(abs(q+h), 1.0/3.0);
        float u = sign(q-h)*pow(abs(q-h), 1.0/3.0);
        float rx = -s - u - c*4.0 + 2.0*m2;
        float ry = (s - u)*sqrt(3.0);
        float rm = sqrt( rx*rx + ry*ry );
        co = (ry/sqrt(rm-rx)+2.0*g/rm-m)/2.0;
    }
    vec2 r = ab * vec2(co, sqrt(1.0-co*co));
    return length(r-p) * sign(p.y-r.y);
}

float tube(vec3 p)
{
    
    float h = 4., r = .8;
    p.z += h/2.;
    if (p.z < 0.)
        return length(vec3(p.z, max(abs(p.x)-r, 0.), p.y));
    if (p.z > h)
        //return length(vec2(p.z-h, max(length(p.xy)-r, 0.)));
        return length(p - vec3(0, 0, h-r)) - (r-0.01)*sqrt(2.);
    
    return sdEllipse(p.xy, vec2(r, r*p.z/h));    

   
}

float map(vec3 p) {
    //return umb(p);
    //return halloween(p);
    //return quin(p);
    //return ring(p, 0.5, 1.2, 0.001) - 0.01;
    //return grid(p);
    //return rose(p);
    //return noiseSph(p);
    //return cylinder2(p);
    //return insole (p);
    //return strep(p);
    //return min(heel(p), insole (p));
    //return sh(p.xzy); //heels
    return tube(p)-0.03;
    

}

// https://iquilezles.org/articles/normalsSDF
vec3 calcNormal(in vec3 pos) {
    const float h = 0.0001; // replace by an appropriate value
    const vec2 k = vec2(1, -1);
    return normalize(k.xyy * map(pos + k.xyy * h) +
        k.yyx * map(pos + k.yyx * h) +
        k.yxy * map(pos + k.yxy * h) +
        k.xxx * map(pos + k.xxx * h));
}

vec3 GetRayDir(vec2 uv, vec3 p, vec3 l, float z) {
    vec3 f = normalize(l - p), r = normalize(vec3(f.z, 0, -f.x)), u = cross(f, r), c = f * z, i = c + uv.x * r + uv.y * u;
    return normalize(i);
}
/*
#if HW_PERFORMANCE==0
#define AA 1
#else
#define AA 2
#endif
*/
#define AA 1

void mainImage(out vec4 fragColor, in vec2 fragCoord) {
    vec3 light = normalize(vec3(0.0, .0, 1.)); //light
    vec3 light2 = normalize(vec3(0.0, 0.0, -1.)); //light
    vec2 mo = 1.5 * cos(0.5 * iTime + vec2(0, 11));
    //if  (iMouse.z > 0.0)
    {
        mo = (-iResolution.xy + 2.0 * (iMouse.xy)) / iResolution.y;
    }
    vec3 ro = vec3(0.0, 0.0, 5.); // camera
    //camera rotation
    ro.yz *= rot(mo.y);
    ro.xz *= rot(-mo.x - 1.57);

    const float fl = 1.5; // focal length
    float dist = dist_infin;

    vec3 b1 = vec3(0.0509, 0.2980, 0.4705), b2 = vec3(0.3764, 0.7529, 0.8784), bg = mix(b1, b2, vec3((1.0 - abs(fragCoord.x - iResolution.x / 2.) / iResolution.y * 2.) * fragCoord.y / iResolution.x));   
    //antialiasing
    vec3 tot = vec3(0.0);
    for(int m = 0; m < AA; m++) for(int n = 0; n < AA; n++) {
            vec2 o = vec2(float(m), float(n)) / float(AA) - 0.5;
            vec2 p = (-iResolution.xy + 2.0 * (fragCoord + o)) / iResolution.y;
            vec3 rd = GetRayDir(p, ro, vec3(0, 0., 0), fl); //ray direction
            //vec3 col = bg * bg; // background  
            vec3 col = b2*b2;
            //vec3 col = wave(p + 5.);
            //==========================raymatch=============================
            float td = 0.;
            vec3 pos = vec3(0.);
            for(int i = 0; i < nn; i++) {
                pos = ro + rd * td;
                float h = map(pos);
                if(h < eps || td >= dist_infin)
                    break;
                td += h;
            }
            if(td < dist_infin) {
                col = resColor;
                vec3 nor = calcNormal(pos);
                vec3 R = reflect(light, nor);
                float specular = pow(max(abs(dot(R, rd)), 0.), 16.);
                float difu = abs(dot(nor, light));
                col = col * (clamp(difu, 0., 1.0) + 0.5) + vec3(1., .7, .4) * specular;
                float fre = pow(clamp(dot(nor, rd) + 1., .0, 1.), 3.); // Fresnel, for some mild glow.
                col += vec3(.1, .1, 0.1) * fre; //?
                col = sqrt(col);
            }
            //==========================raymatch=============================
            tot += col;
        }
    tot = tot / float(AA) / float(AA);
    //tot = tot / float(AA);
    //antialiasing
    fragColor = vec4(tot, 1.0);
}
/////=====================================================================================
void main() {
    vec4 fragColor = vec4(0);
    mainImage(fragColor, gl_FragCoord.xy);
    gl_FragColor = fragColor;
}