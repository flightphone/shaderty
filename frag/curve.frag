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

#define PI  3.14159265359
#define TAU 6.28318530718

const float dist_infin = 100.0;
#define nn 128
const float eps = 0.001;

float sdPolar(vec3 p, float a) {
    float fi = atan(p.y, p.x);//aafi(p.xy);
    float r = a * cos(fi) * cos(fi);

    float L = length(p.xy);
    float d = abs(L - r);
    float fi2 = acos(sqrt(L / a));
    float fi3 = PI - fi2;
    d = min(2.0 * abs(sin((fi - fi2) / 2.0)) * L, d);
    d = min(2.0 * abs(sin((fi - fi3) / 2.0)) * L, d);
    d = min(2.0 * abs(sin((fi + fi2) / 2.0)) * L, d);
    d = min(2.0 * abs(sin((fi + fi3) / 2.0)) * L, d);

    d = sqrt(p.z * p.z + d * d);
    d *= .5;
    d -= 0.03;
    return d;
}

float sdBacket(vec3 p, float a, float b, float m, float n) {
    float fi = atan(p.y, p.x); //aafi(p.xy)
    float w = 20.;
    for(float i = 0.; i < 10.; i++) {
        if(mod(i, m) == 0.0 && i > 0.)
            break;
        float t = (fi + TAU * i) / m;
        float wt = abs(p.z - b * sin(n * t));
        w = min(w, wt);
    }
    float r = length(vec2(length(p.xy) - a, w)) / 2.0;
    return r - 0.03;
}
float smin2(float d1, float d2) {

    return d2;
    //return d1*d2 / length(vec2(d1, d2));

    //return min(d2, d1);
}

float sdCosNp(vec3 p, float a) {

    float fi = atan(p.y, p.x);
    float L = length(p.xy);
    float d1 = dist_infin;
    float d2 = dist_infin;
    float n = 2.;

    float r = a * cos(n * fi);
    if(p.x < 0.)
        r = 0.;
    d1 = min(abs(L - r), d1);

    float f = acos(L / a) / n;

    d2 = min(2.0 * abs(sin((fi - f) / 2.0)) * L, d2);
    d2 = min(2.0 * abs(sin((fi + f) / 2.0)) * L, d2);

    float d = min(d1, d2);
    d = sqrt(p.z * p.z + d * d);
    d *= .4;
    d -= 0.01;
    return d;
}

float sdCosN(vec3 p, float a, float n) {

    float fi = atan(p.y, p.x);
    float L = length(p.xy);
    float d1 = dist_infin;
    float d2 = dist_infin;

    for(float i = 0.; i < 10.; i++) {
        if(fract(n * i) == 0. && i > 0.)
            break;
        float r = a * cos(n * fi + n * i * TAU);
        d1 = min(abs(L - r), d1);
    }
    float f = acos(L / a) / n;
    for(float i = 0.; i < 10.; i++) {

        d2 = min(2.0 * abs(sin((fi - (f + i * TAU / n)) / 2.0)) * L, d2);
        d2 = min(2.0 * abs(sin((fi + (f + i * TAU / n)) / 2.0)) * L, d2);
        d2 = min(2.0 * abs(sin((fi - (f - i * TAU / n)) / 2.0)) * L, d2);
        d2 = min(2.0 * abs(sin((fi + (f - i * TAU / n)) / 2.0)) * L, d2);

        /*
        d2 = min(abs(sin((fi - (f + i * TAU / n)))) * L, d2);
        d2 = min(abs(sin((fi + (f + i * TAU / n)))) * L, d2);
        d2 = min(abs(sin((fi - (f - i * TAU / n)))) * L, d2);
        d2 = min(abs(sin((fi + (f - i * TAU / n)))) * L, d2);
        */
    }
    //float d = smin2(d1, d2);
    float d = min(d1, d2);
    d = sqrt(p.z * p.z + d * d);
    d *= .5;
    d -= 0.02;
    return d;
}
//https://iquilezles.org/articles/smin/
// polynomial smooth min 1 (k=0.1)
float smin(float a, float b, float k) {
    float h = clamp(0.5 + 0.5 * (b - a) / k, 0.0, 1.0);
    return mix(b, a, h) - k * h * (1.0 - h);
}

float sdCosN3(vec3 p, float a, float n) {

    float fi = atan(p.y, p.x);
    float L = length(p.xy);
    float ro = dist_infin;
    float d = dist_infin;
    float d2 = dist_infin;
    fi += step(p.y, 0.0) * TAU;

    float f = acos(L / a) / n;
    for(float i = 0.; i < 6.; i++) {
        d = min(2. * abs(sin((fi - (f + i * TAU / n)) / 2.0)) * L, d);
        d = min(2. * abs(sin((fi + (f + i * TAU / n)) / 2.0)) * L, d);
        d = min(2. * abs(sin((fi - (f - i * TAU / n)) / 2.0)) * L, d);
        d = min(2. * abs(sin((fi + (f - i * TAU / n)) / 2.0)) * L, d);
    }

    for(float i = 0.; i < 5.; i++) {
        if(fract(n * i) == 0. && i > 0.)
            break;
        float r = a * cos(n * fi + n * i * TAU);
        float l = abs(L - r);
        if(l < d2) {
            d2 = l;
            ro = r;
        }

    }
    /*
    float f2 = atan(d, d2);
    d = d2 * sin(f2);
    ro = d2*cos(f2)*cos(f2) + ro;
    */

    //d = smin(d, d2, 0.2);
    d = min(d, d2);
    if(d < d2)
        ro = L;

    float z = p.z;// - 0.4 * ro * ro;
    d = sqrt(z * z + d * d);
    d *= .3;
    d -= 0.015;
    return d;
}
float aafi(vec2 p) {
    float fi = atan(p.y, p.x);
    fi += step(p.y, 0.0) * TAU;
    return fi;
}
mat3 rotateX(float f) {
    return mat3(vec3(1.0, 0.0, 0.0), vec3(0.0, cos(f), -sin(f)), vec3(.0, sin(f), cos(f)));
}

mat3 rotateZ(float f) {
    return mat3(vec3(cos(f), -sin(f), 0.0), vec3(sin(f), cos(f), 0.0), vec3(0.0, 0.0, 1.0));

}

mat3 rotateY(float f) {
    return mat3(vec3(cos(f), 0.0, sin(f)), vec3(0.0, 1.0, 0.0), vec3(-sin(f), 0.0, cos(f)));
}
float fround(float i) {
    float r = floor(i);
    if(ceil(i) >= 0.5)
        r += 1.;
    //r += step(0.5, ceil(i));
    return r;
}

float sdRound(vec3 p, float r, float f) {
    float d = abs(length(p.xy) - r * cos(f));
    d = length(vec2(p.z - r * sin(f), d));
    return d;
}
float sdLonLat(vec3 p, float r) {
        #define ll 20.
    float fi = atan(p.x, p.y);
        //fi += step(p.y, 0.0)*TAU;
    fi = mod(fi, TAU);
    float ln = floor(fi / TAU * ll);
    float l1 = ln * TAU / ll;
    float l2 = l1 + TAU / ll;
    float d = min(sdRound(rotateX(l1) * rotateY(PI / 2.) * p, r, 0.), sdRound(rotateX(l2) * rotateY(PI / 2.) * p, r, 0.));

    fi = atan(p.z, length(p.xy));
    float mm = ll / 4.0;
    ln = floor(abs(fi) / PI * 2.0 * mm);
    l1 = ln * PI / 2.0 / mm;
    l2 = l1 + PI / 2.0 / mm;
    float d2 = min(sdRound(p, r, l1 * sign(p.z)), sdRound(p, r, l2 * sign(p.z)));
    d = min(d2, d);
    return d - 0.03;
}

float sdCoil(vec3 p) {
    float k = 0.01;

    float r = 0.1;
    float f = atan(p.y, p.x);
    f = mod(f, TAU);
    float z = dist_infin;
    float d = dist_infin;

    float nc = 38.;

    float ii = abs(p.z - f * k) / (TAU * k);
    float i = floor(ii);
    i = clamp(i, 0., nc - 2.);
    z = min(abs(p.z - f * k - TAU * k * i), abs(p.z - f * k - TAU * k * (i + 1.)));
    d = abs(length(p.xy) - r);
    d = length(vec2(d, z));
    d = min(d, length(p - vec3(r, 0., 0.)));
    d = min(d, length(p - vec3(r, 0., TAU * k * nc)));
    d *= 0.5;
    return d - 0.015;

}

float sdArc(vec3 p) {
    float ra = 0.02;
    float k = 0.01;
    float e = 0.9;
    float L = length(p.xy);
    float fi = aafi(p.xy);
    float d = dist_infin;
    float nc = 20.0;

    float ii = abs(L - k * fi) / (TAU * k);
    float i = floor(ii);
    i = clamp(i, 1., nc - 2.);
    d = min(abs(L - k * (fi + TAU * i)) * e, abs(L - k * (fi + TAU * (i + 1.))) * e);

    float r = length(p.xy - vec2(k * TAU * nc, 0.0)) * e;
    d = min(d, r);
    r = length(p.xy - vec2(k * TAU, 0.0)) * e;
    d = min(d, r);

    d = sqrt(p.z * p.z * e * e + d * d);
    d -= ra;
    return d;
}

float sdCeli(vec3 p) {
    float ra = 1.;
    float n = 15.;
    float l = length(p.xy);
    float lon = mod(atan(p.y, p.x), TAU);
    float lat = atan(p.z, l);
    //lat += PI/2.;
    lat += PI;

    float x = fract(lon / PI * n);
    float f = floor(lon / PI * n);
    float sg = mod(f, 2.0);
    float lat1 = 0.;
    float lon1 = 0.;
    float e = 1.;
    float d = dist_infin;
    float d2 = dist_infin;
    if(sg == 0.) {
        lat1 = x * PI;
        lon1 = lat / n + f * PI / n;
    } else {
        lat1 = PI * (1. - x);
        lon1 = (PI - lat) / n + f * PI / n;
    }
    d = smin2(ra * sin(abs(lat - lat1) / 2.) * 2., (l * sin(abs(lon - lon1) / 2.) * 2.));

    float r = length(vec2(abs(length(p) - ra), d));
    return r * 0.6 - 0.02;
}

float sdCeli2(vec3 p) {
    float ra = 1.0;
    float n = 2.;
    float l = length(p.xy);
    float lon = mod(atan(p.y, p.x), TAU);
    float lat = atan(p.z, l) + PI / 2.;

    float x = fract(lon / PI * n);
    float f = floor(lon / PI * n);
    float sg = mod(f, 2.0);
    float lat1 = 0.;
    float lon1 = 0.;

    if(sg == 0.) {
        lat1 = x * PI;
        lon1 = lat / n + f * PI / n;
    } else {
        lat1 = PI * (1. - x);
        lon1 = (PI - lat) / n + f * PI / n;
    }

    float d = smin2(ra * sin(abs(lat - lat1) / 2.) * 2., l * sin((abs(lon - lon1)) / 2.) * 2.);

    float r = length(vec2(abs(length(p) - ra), d));
    return r * 0.5 - 0.03;

}

float sdCeli1(vec3 p) {
    float ra = .8;
    float n = 7.;
    float l = length(p.xy);
    float lon = mod(atan(p.y, p.x), TAU);
    float lat = atan(p.z, l) + PI / 2.;
    //lat += PI;

    /*
    float d2 = abs(lat + PI/4.);
    d2 = smoothstep(0.01, 0., d2)*0.01;

    float d3 = abs(lon + PI/2.);
    d3 = smoothstep(0.01, 0., d3)*0.01;
    
    float d4 = abs(lon - PI/2.);
    d4 = smoothstep(0.01, 0., d4)*0.01;
    */
    float r = (length(p) - ra);
    /*
    float d = abs(lon);
    d = smoothstep(0.02, 0., d)*0.02;
    r += d;
    d = abs(lon-PI);
    d = smoothstep(0.02, 0., d)*0.02;
    r+=d;
    d = abs(lat-PI/2.);
    d = smoothstep(0.02, 0., d)*0.02;
    r-=d;
    */
    float d = abs(lat - mod(n * lon, PI));
    d = smoothstep(0.02, 0., d) * 0.02;
    r -= d;
    return r;
}
float tani(float d1, float d2, float t) {
    if(cos(t) * sin(t) > 0.) {
        if(cos(t) > 0.)
            return d1;
        else
            return d2;
    } else {
        if(cos(t) > 0.)
            return d2 + (d1 - d2) * cos(t);
        else
            return d1 + (d1 - d2) * cos(t);

    }
}
float map(in vec3 pos) {
    //return sdCeli2(pos);
    //return sdCosNp(pos, 1.0);
    //float d1 =   sdBacket(pos, .9, .4, 3.5, 9.); 
    //float d2 = sdCeli2(pos);
    //float d2 = sdLonLat(pos, 1.);
    //float d = tani(d1, d2, iTime/2.);
    //return d;
    //return sdCeli(pos);
    //return sdBacket(pos, .9, .4, 3.5, 9.); 

    //float d2 = sdCeli(-pos);
    //return min(d, d2);
    //return sdPolar(pos, 1.0);

    //return sdCosN3(pos, 1.0, .2);
    //eturn sdCosN3(pos, 1.0, 0.25);
    //return sdCosN3(pos, 1.0, 0.4);
    //return sdCosN3(pos, 1.0, 0.5);
    //return sdCosN3(pos, 1.0, 1.5);
    //return sdCosN3(pos, 1.0, 2.);
    //return sdCosN3(pos, 1.0, 2.2);
    //return sdCosN(pos, 1.0, 2.2);

    //return sdLonLat(pos, 1.0);

    /*
    float d = sdLonLat(pos, 1.0);
    float d2 =  sdCoil(vec3(pos.x, pos.y, pos.z+1.2));
    return min(d, d2);
    */

    /*    
    float d = sdArc(pos);
    d = min(d, sdCoil(vec3(pos.x, pos.y, pos.z + 1.2)));
    d = min(d, sdLonLat(pos, .9));
    return d;
    */
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

struct HIT {
    float dist;
    vec3 nor;
    vec3 pos;
};
const HIT hit_inf = HIT(dist_infin, vec3(0.0), vec3(0.0));

HIT giper3D(vec3 ro, vec3 rd) {
    float t = 0.;
    for(int i = 0; i < nn; i++) {
        vec3 pos = ro + rd * t;
        float h = map(pos);
        if(h < eps || t >= dist_infin)
            break;
        t += h;
    }

    if(t >= dist_infin)
        return hit_inf;

    vec3 pos = ro + t * rd;
    vec3 nor = calcNormal(pos);
    return HIT(t, nor, pos);
}

vec3 GetRayDir(vec2 uv, vec3 p, vec3 l, float z) {
    vec3 f = normalize(l - p), r = normalize(cross(vec3(0, 1, 0), f)), u = cross(f, r), c = f * z, i = c + uv.x * r + uv.y * u;
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
// IQ's vec2 to float hash.
float hash21(vec3 p){  
    return fract(sin(mod(dot(p, vec3(27.609, 57.583, 11.2345)), 6.2831853))*43758.5453); 
}

vec3 point(vec3 p) {
    float f = atan(p.y, p.x);
    float d = atan(length(p.xy), p.z);
    float n = 200.; 
    float m = 200.;
    f = floor(f*n)/n;
    d = floor(d*m)/m;
    return vec3(sin(d)*cos(f), sin(d)*sin(f), cos(d));
}

void mainImage(out vec4 fragColor, in vec2 fragCoord) {
    vec3 light = normalize(vec3(0.0, 1.0, 1.0)); //light
    vec3 light2 = normalize(vec3(0.0, 1.0, -1.0)); //light
    float t = iTime / 3.;
    vec2 m = vec2(0.0, 0.0);
    //if  (iMouse.z > 0.0)
    {
        m = (-iResolution.xy + 2.0 * (iMouse.xy)) / iResolution.y;
        t = 0.;
    }
    vec3 ro = vec3(0.0, 0.0, 2.0); // camera
    ro = rotateY(-m.x * TAU) * rotateX(-m.y * PI) * ro; //camera rotation

    const float fl = 1.5; // focal length
    float dist = dist_infin;
    mat3 rota = rotateZ(t) * rotateX(-t);
    mat3 rota_1 = rotateX(t) * rotateZ(-t);

    vec3 tot = vec3(0.0);

    //antiblick
    for(int m = 0; m < AA; m++) for(int n = 0; n < AA; n++) {
            vec2 o = vec2(float(m), float(n)) / float(AA) - 0.5;
            vec2 p = (-iResolution.xy + 2.0 * (fragCoord + o)) / iResolution.y;
            vec3 rd = GetRayDir(p, ro, vec3(0, 0., 0), fl); //ray direction
            vec3 col = vec3(0.0);
            HIT giper = giper3D(rota * ro, rota * rd);
             if(giper.dist < dist) {
                    
                  //float  ff = hash21(point(giper.pos));
                  //if (ff > 0.9)
                  {
                    //col = vec3(1.);

                    
                    vec3 nor = rota_1 * giper.nor;
                    float dif = clamp(dot(nor, light), 0.2, 1.0);
                    float amb = 0.5 + 0.5 * dot(nor, light2);
                    col = vec3(0.2, 0.3, 0.4) * amb + vec3(0.85, 0.75, 0.65) * dif;
                    
                  }   
                    
                    
                }
            
        // gamma        
            col = sqrt(col);
            tot += col;
        }
    //antiblick
    tot /= float(AA * AA);
    fragColor = vec4(tot, 1.0);
}
/////=====================================================================================
void main() {
    vec4 fragColor = vec4(0);
    mainImage(fragColor, gl_FragCoord.xy);
    gl_FragColor = fragColor;
}