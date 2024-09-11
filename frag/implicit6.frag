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
float hash(float p) {
    p = fract(p * 0.011);
    p *= p + 7.5;
    p *= p + p;
    return fract(p);
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

float n3D(vec3 p) {
    return noise(p);
}
/////=====================================================================================

#define PI  3.14159265359
#define TAU 6.28318530718
#define rot(f) mat2(cos(f), -sin(f), sin(f), cos(f))
#define nn 128.
#define newton 5

float sdRound(vec2 p, float r) {
    return length(p) - r;
}
float sdPentagon(in vec2 p, in float r) {
    const vec3 k = vec3(0.809016994, 0.587785252, 0.726542528);
    p.x = abs(p.x);
    p -= 2.0 * min(dot(vec2(-k.x, k.y), p), 0.0) * vec2(-k.x, k.y);
    p -= 2.0 * min(dot(vec2(k.x, k.y), p), 0.0) * vec2(k.x, k.y);
    p -= vec2(clamp(p.x, -r * k.z, r * k.z), r);
    return length(p) * sign(p.y);
}

float gr(vec2 p) {
    float t = sdPentagon(p, 1.);
    t = smoothstep(0., 0.05, t + 0.15) * smoothstep(0., 0.05, 0.15 - t);
    t *= 2.;

    float t2 = sdRound(p, 0.7);
    t2 = smoothstep(0., 0.05, t2 + 0.1) * smoothstep(0., 0.05, 0.1 - t2);
    t2 *= 2.;
    t += t2;
    return -t;
}

float sdRound3(vec3 p) {
    float R0 = 2.0, d2 = dot(p.xy, p.xy) - R0 * R0;
    if(p.z > 0.)
        d2 += 8. * p.z * p.z * p.z + gr(p.xy);
    else
        d2 += 50. * p.z * p.z;

    //float x = p.x, y = p.z, z = p.y, 
    float r = 0.1, R = 0.6, v = (dot(p, p) + R * R - r * r);
    v = v * v - 4. * R * R * (dot(p.xz, p.xz));
    if(p.z > .0)
        v = 10.;
    return d2 * v - 0.01;

}

float sdHexagon(in vec2 p, in float r) {
    const vec3 k = vec3(-0.866025404, 0.5, 0.577350269);
    p = abs(p);
    p -= 2.0 * min(dot(k.xy, p), 0.0) * k.xy;
    p -= vec2(clamp(p.x, -k.z * r, k.z * r), r);
    return length(p) * sign(p.y);
}

float sdBox(in vec2 p, in vec2 b) {
    vec2 d = abs(p) - b;
    return length(max(d, 0.0)) + min(max(d.x, d.y), 0.0);
}

float sdEquilateralTriangle(in vec2 p, in float r) {
    const float k = sqrt(3.0);
    p.x = abs(p.x) - r;
    p.y = p.y + r / k;
    if(p.x + k * p.y > 0.0)
        p = vec2(p.x - k * p.y, -k * p.x - p.y) / 2.0;
    p.x -= clamp(p.x, -2.0 * r, 0.0);
    return -length(p) * sign(p.y);
}

float ringdisp3(vec3 p) {
    float x = mod(atan(p.y, p.x), TAU), r = 0.3, y = 1.5 * p.z / 0.3, f = 10.;
    float t = y - sin(x * f), t2 = y - sin((x + PI / 2.) * f);
    t = smoothstep(0., 0.05, 0.1 - t) * smoothstep(0., 0.05, 0.1 + t);
    t2 = smoothstep(0., 0.05, 0.1 - t2) * smoothstep(0., 0.05, 0.1 + t2);
    return (t + t2) * 0.1;
}

float ringdisp(vec3 p) {

    float n = 5., r = 0.3, R = 2., x = p.z / r; 
    //y = p.y / r;
    float fi = mod(atan(p.y, p.x), TAU), dlon = TAU / n, i = floor(fi / dlon), fi1 = i * dlon + dlon / 2.;
    float y = (fi1 - fi) * 10.;
    float t = 0.;

    if(i == 0.) {
        t = sdPentagon(vec2(x, y), 1.);
    }

    if(i == 1.) {
        t = sdHexagon(vec2(x, y), 1.);
    }

    if(i == 2.) {
        t = sdBox(vec2(x, y), vec2(1.0));
    }

    if(i == 3.) {
        t = sdEquilateralTriangle(vec2(x, y), 1.0);
    }

    if(i == 4.) {
        t = sdRound(vec2(x, y), 1.);
    }

    t = smoothstep(0., 0.06, 0.2 - t);//*smoothstep(0., 0.06, 0.2 + t);

    return -t * 0.8;
}

float ring(vec3 p) {
    float R = 2., r = 0.3, x = p.z, y = length(p.xy) - R;
    float res = -1.;
    if(y > 0.) {
        res = x * x + y * y * 4. - r + ringdisp(p) - 0.05;
    } else
        res = x * x + y * y * 20. - r + ringdisp3(p) - 0.05;

    return res;

}

float dish(vec3 p) {

    float res = 1., R = 3., x = p.x, y = p.y, z = p.z, r = 0.2;
    if(p.z < 0.)
        res = 10.;
    else {
        res = (x * x + y * y + z * z * 2. - (R) * (R)) * (x * x + y * y + z * z * 4. - (R) * (R));
    }
    
    return res + 1.;
}

float map(vec3 p) {
    //return sdRound3(p);
    //return ring(p);
    return dish(p);
    //return trefoil(p);

}

vec3 calcNormal(in vec3 p) {
    const float eps = 0.0001;
    vec2 q = vec2(0.0, eps);
    vec3 res = vec3(map(p + q.yxx) - map(p - q.yxx), map(p + q.xyx) - map(p - q.xyx), map(p + q.xxy) - map(p - q.xxy));
    return normalize(res);
}

vec3 getPoint(vec3 a, vec3 b, float v0, float v1) {
    vec3 m;
    //binary search with  n iterations, n = newton
    for(int i = 0; i < newton; i++) {
        m = (a + b) * 0.5;
        float v = map(m);
        if(v == 0.)
            break;

        if(sign(v) * sign(v0) <= 0.) {
            v1 = v;
            b = m;
        } else {
            v0 = v;
            a = m;
        }
    }
    return m;
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
//shane color
//https://www.shadertoy.com/view/4td3zj
// I keep a collection of occlusion routines... OK, that sounded really nerdy. :)
// Anyway, I like this one. I'm assuming it's based on IQ's original.
float calculateAO(in vec3 p, in vec3 n) {
    float sca = 2., occ = 0.;
    for(float i = 0.; i < 5.; i++) {

        float hr = .01 + i * .5 / 4.;
        float dd = map(n * hr + p);
        occ += (hr - dd) * sca;
        sca *= 0.7;
    }
    return clamp(1.0 - occ, 0., 1.);
}

// Simple environment mapping. Pass the reflected vector in and create some
// colored noise with it. The normal is redundant here, but it can be used
// to pass into a 3D texture mapping function to produce some interesting
// environmental reflections.
vec3 envMap(vec3 rd, vec3 sn) {

    vec3 sRd = rd; // Save rd, just for some mixing at the end.

    // Add a time component, scale, then pass into the noise function.
    rd.xy -= iTime * .25;
    rd *= 3.;

    float c = n3D(rd) * .57 + n3D(rd * 2.) * .28 + n3D(rd * 4.) * .15; // Noise value.
    c = smoothstep(.4, 1., c); // Darken and add contast for more of a spotlight look.

    vec3 col = vec3(c, c * c, c * c * c * c); // Simple, warm coloring.
    //vec3 col = vec3(min(c*1.5, 1.), pow(c, 2.5), pow(c, 12.)); // More color.

    // Mix in some more red to tone it down and return.
    return mix(col, col.yzx, sRd * .25 + .25);

}

#define AA 1
void mainImage(out vec4 fragColor, in vec2 fragCoord) {

    //float dist_infin = 1.8;
    //float hh = 3.5;

    float dist_infin = 2.8;
    float hh = 4.8;

    vec3 light = normalize(vec3(.0, 1.0, -2.5)); //light
    //vec2 mo = 1.5*cos(0.5*rottime() + vec2(0,11));
    vec2 mo = 1.5 * cos(0.5 * iTime + vec2(0, 11));
    vec3 ro = vec3(0.0, 0.0, hh); // camera
    //if  (iMouse.z > 0.0)
    {
        mo = (-iResolution.xy + 2.0 * (iMouse.xy)) / iResolution.y;
        ro.yz *= rot(mo.y * 2.);
        ro.xz *= rot(-mo.x * 1.5);
    }
    /*
    else
    {
        ro.yz *= rot(mo.y);
        ro.xz *= rot(-mo.x - 1.57);
    }
    */
    //camera rotation

    // Light position, hovering around behind the camera.
    // https://www.shadertoy.com/view/4td3zj
    float tm = iTime / 2.;
    vec3 lp = ro + vec3(cos(tm / 2.) * .5, sin(tm / 2.) * .5, -.5);

    const float fl = 1.5; // focal length
    float dist = dist_infin;

    vec3 bg =  vec3(0.7, 0.7, 0.9) * 0.6; //vec3(0.); //
    vec3 col1 = vec3(1, .8, .4);//vec3(0.73, 0.59, 0.3);
    

    //antialiasing
    vec3 tot = vec3(0.0);
    for(int m = 0; m < AA; m++) for(int n = 0; n < AA; n++) {
            vec2 o = vec2(float(m), float(n)) / float(AA) - 0.5;
            vec2 p = (-iResolution.xy + 2.0 * (fragCoord + o)) / iResolution.y;
            vec3 rd = GetRayDir(p, ro, vec3(0, 0., 0), fl); //ray direction
            vec3 col = bg; // background  

            //STEP 1. Calculating bounding sphere
            float d = length(cross(ro, rd));
            if(d >= dist) {
                tot += col;
                continue;
            }
            /*
            STEP 2.
            ray tracing inside the bounding sphere, 
            searching for a segment with different signs of the function value 
            at the ends of the segment
            */
            float td = abs(dot(ro, rd));
            d = sqrt(dist * dist - d * d);
            vec3 pos0 = ro + rd * (td - d);
            vec3 pos1 = ro + rd * (td + d);
            vec3 rd0 = pos1 - pos0;
            vec3 pos = pos0;
            float val0 = map(pos0);
            for(float i = 1.; i < nn; i++) {
                pos1 = pos0 + rd0 * i / (nn - 1.);
                float val1 = map(pos1);
                if(sign(val0) * sign(val1) <= 0.) {
                    //different signs of the function value  at the ends of the segment
                    //STEP 3. binary search to clarify the intersection of a ray with a surface.
                    pos = getPoint(pos, pos1, val0, val1);
                    vec3 nor = calcNormal(pos);
                    //shane color
                    //https://www.shadertoy.com/view/4td3zj
                    // Surface postion, surface normal and light direction.
                    vec3 sp = pos;
                    vec3 sn = nor;
                    vec3 ld = lp - sp;
                    vec3 oC = col1;
                    
                    float tx = noise(pos*2.);
                    tx = fract(tx*5.);
                    tx = smoothstep(0., 0.01, tx-0.5);
                    oC*=tx;
                     

                    float lDist = max(length(ld), 0.001); // Light distance.
                    float atten = 1. / (1. + lDist * .125); // Light attenuation.

                    ld /= lDist; // Normalizing the light direction vector.

                    float diff = max(dot(ld, sn), 0.); // Diffuse.
                    float spec = pow(max(dot(reflect(-ld, sn), -rd), 0.0), 16.); // Specular.
                    float fre = pow(clamp(dot(sn, rd) + 1., .0, 1.), 3.); // Fresnel, for some mild glow.

                    // Shading. Note, there are no actual shadows. The camera is front on, so the following
                    // two functions are enough to give a shadowy appearance.
                    //crv = crv * .9 + .1; // Curvature value, to darken the crevices.
                    float crv = 1., edge = 0.0;
                    float ao = calculateAO(sp, sn); // Ambient occlusion, for self shadowing.

                    // Combining the terms above to light the texel.

                    col = oC * (diff + .5) + vec3(1., .7, .4) * spec * 2. + vec3(.4, .7, 1) * fre;

                    //col += (oC * .5 + .5) * envMap(reflect(rd, sn), sn) * 6.; // Fake environment mapping.

                    // Edges.
                    col *= 1. - edge * .85; // Darker edges.   

                    // Applying the shades.
                    col *= (atten * crv * ao);

                    // Rough gamma correction, then present to the screen.
                    col = sqrt(clamp(col, 0., 1.));

                    break;
                }
                
                val0 = val1;
                pos = pos1;
            }
            tot += col;
        }
    tot = tot / float(AA) / float(AA);
    //antialiasing
    fragColor = vec4(tot, 1.0);
}
/////=====================================================================================
void main() {
    vec4 fragColor = vec4(0);
    mainImage(fragColor, gl_FragCoord.xy);
    gl_FragColor = fragColor;
}