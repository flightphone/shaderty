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

//Fork of @Efim's Christmas Tree https://www.shadertoy.com/view/lcfGW8

// SNOW background from @g1mishr's "Snow Simple " https://www.shadertoy.com/view/DlGczD
#define TILES 10.0

//2D random from https://www.shadertoy.com/view/WstGDj
float random (vec2 uv) {
    return fract(sin(dot(uv, vec2(135., 263.))) * 103.214532);
}

vec4 drawSnow(vec2 curid, vec2 uv, vec4 fragColor, float r, float c)
{
    float maxoff = 2.0 / TILES; //calculate the max offset a particle can have (two tiles)

    //loop through neighboring tiles
    for(int x=-2; x<=1; x++)
    {
        for(int y=-2; y<=0; y++)
        {
            float rad = (1.0 / (TILES * 5.0)) * r; //set default radius
            vec2 id = curid + vec2(x, y); //get the id of the tile we're visiting
            vec2 pos = id / TILES; //calculate position
            float xmod = mod(random(pos), maxoff);
            pos.x += xmod; //add a random x-offset
            pos.y += mod(random(pos+vec2(4,3)), maxoff); //add a random y-offset
            rad *= mod(random(pos), 1.0); //vary the radius by multiplying by a random val
            pos.x += 0.5*(maxoff-xmod)*sin(iTime*r + random(pos)*100.0); //dynamic sin wave x-offset
            
            float len = length(uv - pos); //calculate distance from tile's particle

            //if we're inside the particle, draw it
            float v = smoothstep(0.0, 1.0, (rad - len) / rad*0.75);
            fragColor = mix(fragColor, vec4(c), v);      
        }
    }
    
    return fragColor;
}


vec4 snowBackground( vec2 fragCoord )
{
    vec4 fragColor = vec4(0.0);
    vec2 uv = (2.0*fragCoord - iResolution.xy)/iResolution.x;
    uv.y -= 0.3;
    
    //uv.x -= 0.6;

    
    vec3 col = mix(vec3(0.0, 0.45, 0.85), vec3(1), -0.3-uv.y);

    // Output to screen
    fragColor = vec4(col,1.0);
    
    vec4 bg = vec4(.529, .808, .922, 1) * 0.25;
    vec2 uvNorm = fragCoord.xy / iResolution.xy; //normalized UV coordinate [0, 1]
    vec2 uvog = fragCoord.xy / iResolution.y; //UV coordinate (will remain static)
    uv = fragCoord.xy / iResolution.y; //UV coordinate (we'll modify this one)
    
    //draw the closest snow layer
    uv += 0.2*vec2(-iTime, iTime); //move the UV coords based on time
    vec2 curid = floor(uv * TILES); //calculate the ID associated with the current UV
    curid += vec2(0.5); //center the ID
    
    //if(curid.y > 10.0)
    {
    fragColor = drawSnow(curid, uv, fragColor, 1.0, 0.9); //draw closest snow layer
    
    //draw the middle snow layer, calculate new UV and ID
    uv = uvog + 0.1*vec2(-iTime - 100.0, iTime + 100.0);
    curid = floor(uv * TILES);
    curid += vec2(0.5);
    fragColor += drawSnow(curid, uv, vec4(0), 0.75, 0.45); 
    
    //draw the far snow layer, calculate new UV and ID
    uv = uvog + 0.05*vec2(-iTime - 150.0, iTime + 150.0);
    curid = floor(uv * TILES);
    curid += vec2(0.5);
    fragColor += drawSnow(curid, uv, vec4(0), 0.5, 0.225);
    
    //fragColor = smoothstep(0.0, 3.0, iTime)*fragColor;
    }
    return fragColor;
}

// END Snow Simple https://www.shadertoy.com/view/DlGczD



#define PI  3.14159265359
#define TAU 6.28318530718
#define rot(f) mat2(cos(f), -sin(f), sin(f), cos(f))

#define COLOR_RANGE 128.   //360.



vec3 hsb2rgb( in vec3 c )
{
    vec3 rgb = clamp(abs(mod(c.x*6.0+vec3(0.0,4.0,2.0),
                             6.0)-3.0)-1.0,
                     0.0,
                     1.0 );
    rgb = rgb*rgb*(3.0-2.0*rgb);
    return (c.z * mix( vec3(1.0), rgb, c.y));
}



const float dist_infin = 10.0;
#define nn 128
const float eps = 0.001;

vec3 sdfColor;
vec3 resColor;

vec3 col1 = vec3(0.024,0.357,0.153);  // base green
//vec3 col2 = vec3(0.412,0.827,0.439);  //snow highlight
vec3 col2 = vec3(0.7686274509803922, 0.8235294117647058, 0.8745098039215686);
vec3 col3 = vec3(1., 0.8431, 0.);
float sdfReflect = 0.5;
float resReflect = 0.5;

vec3 csky(vec3 p) {
    float n = 5., m = 5., dlat = PI / n, dlon = TAU / m;
    float lon = mod(atan(p.y, p.x), TAU), lat = atan(length(p.xy), p.z);
    float fo = fract(lon / dlon), fa = fract(lat / dlat);

    float pst = fo * fa * (1. - fo) * (1. - fa);
    pst = smoothstep(0.0, 0.0625, pst);
    pst = clamp(pst, 0.1, 1.0);
    return vec3(pst);
}

float sdSolidAngle(vec3 p, vec2 c, float ra) {
  // c is the sin/cos of the angle
    vec2 q = vec2(length(p.xz), p.y);
    float l = length(q) - ra;
    float m = length(q - c * clamp(dot(q, c), 0.0, ra));
    return max(l, m * sign(c.y * q.x - c.x * q.y));
}

float heigthBranch(vec2 p) {
    float n = 2.5;
    float df = PI / n / 2.5;
    float fi = atan(p.y, p.x);
    float L = length(p.xy);
    float r = cos(n * fi);
    if(abs(fi) > df)
        r = 0.;
    float d = r - L;
    float h = smoothstep(0., 0.3, d * L * L);
    if (h > 0.)
    {
        
        sdfColor = col1;
        float pst = smoothstep(0.2, 0., abs(L-0.6));
        sdfColor = mix(col1, col2, pst);
        sdfReflect = mix(0.2, 0., pst);
    }
    return h;
}

float getlon(float lon, float n, float shift) {
    lon = lon - shift;
    float dlon = TAU / n, lon1 = floor(lon / dlon) * dlon;
    if((lon - lon1) >= dlon / 2.)
        lon1 +=  dlon;
    return lon1 + shift; ////mod(lon1 + shift, TAU);
}

float sdTree(vec3 p, float l, float r) {
    float mfi = PI / 8.;
    float d = sdSolidAngle(p, vec2(sin(mfi), cos(mfi)), l) - r;
    if(p.y < 0. || p.y > l * cos(mfi)) {
        sdfColor = col2;
        sdfReflect = 0.;
        return d;
    }
    sdfColor = col1;
    sdfReflect = 0.1;

    float n = 8., m = 5., nc = 6., dnc = l/nc;
    float lss = l/2./m,  ls = 2.*lss;
    float z = clamp(p.y, 0., l);
    float lon = mod(atan(p.z, p.x), TAU), dlon = TAU / n;

    
    float j = floor(z / lss);
    float h1 = j * lss, shift1 = mod(j, 2.) * dlon / 2.;//,h2 = h1 + lss, shift2 = mod((j + 1.), 2.) * dlon / 2.;
    float h3 = h1 - lss, shift3 = mod(j - 1., 2.) * dlon / 2.;//h4 = h1 - 2.*lss, shift4 = mod((j - 2.), 2.) * dlon / 2.;

    float lon1 = getlon(lon, n, shift1);//, lon2 = getlon(lon, n, shift2);
    float lon3 = getlon(lon, n, shift3);//, lon4 = getlon(lon, n, shift4);
    
    float h = 0.;
    if (j < n && h1 > 0.)
        h = max(heigthBranch(vec2((p.y - h1)/ls, (lon-lon1)/dlon*0.5))*l, h);
    if (h3 > 0. && h3 + ls < l)
        h = max(heigthBranch(vec2((p.y - h3)/ls, (lon-lon3)/dlon*0.5))*l, h);
    

    float blink=1.0-cos(5.0*2.0*iTime);
    float glowFact = 0.;//blink*0.01;
    
	vec3 glowColor = col1;
	
	//@Efim's correct logic to avoid loop
	float level=j;
	if(h1+lss-z<z-h1) 
		level = j+1.;
	if(level>2.0) {
    
    ///@Efim's logic computes if there are any iterations of the loop that do anything - and avoids the loop otherwise
    ///  very nice!
    ///for(float level=2.0;level<10.0;level+=1.) {
        //shpere
        float hp = /*6.*/level*lss;
        float offset = mod(level,2.0) * mfi;
        float lonsp = getlon(lon, n, offset);
        float dx = hp*tan(mfi)*(lon - lonsp);
        float dy = (z - hp)/cos(mfi), dr = l/15.;
        float ra = length(vec2(dx, dy)/dr);
        lonsp=abs(lonsp);
        if(lonsp<mfi) lonsp+=mfi;
        float colorfact=lonsp;
        if (ra < 0.4-(0.2/level))
        {
            h = max(sqrt(0.16-ra*ra), h);
            //sdfColor = vec3(0.698*lonsp,0.098,0.176);
            sdfColor = hsb2rgb(vec3(offset+COLOR_RANGE/(level*colorfact),1.,blink+0.9));
            sdfReflect = 0.4;
//        } else if (ra < 0.4+glowFact) {
//           glowColor = hsb2rgb(vec3(330./(level*colorfact*1.5),1.,blink+0.9));
//           sdfReflect = 0.4;
        }
    }

    float pst = smoothstep(0.1, 0., fract(z/dnc));
    sdfColor = mix(sdfColor, col3, pst);
    sdfReflect = mix(sdfReflect, 0.8, pst);
    
    return d * 0.3 - h*0.06*sqrt(z/l);

}

float map(vec3 p) {
    float l = 2.3;
    p.xy *= rot(PI);
    p += vec3(0., l / 2., 0.);
    p.xz *= rot(iTime/2.);
    float d = sdTree(p, l, 0.05);
    resColor = sdfColor;
    resReflect = sdfReflect;
    return d;
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

vec3 calccolor(vec3 col_in, vec3 backcol, vec3 rd, vec3 light1, vec3 light2, vec3 nor) {
    vec3 col = col_in;
    float d = dot(rd, nor);
    if(d < 0.0)
        col = backcol;

    nor *= -sign(d);
    float difu1 = dot(nor, light1);
    float difu2 = dot(nor, light2);
    float difu = max(difu1, difu2);
    col *= clamp(difu, 0.3, 1.0);
    return col;
}

void mainImage(out vec4 fragColor, in vec2 fragCoord) {

    vec3 snowBgcol = snowBackground( fragCoord ).rgb;
    
    vec3 light = normalize(vec3(1.0, .0, -2.5)); //light
    vec3 light2 = normalize(vec3(-1.0, -.0, 2.5)); //light
    vec2 mo = vec2(0.0, 0.0);
    //if  (iMouse.z > 0.0)
    {
        mo = (-iResolution.xy + 2.0 * (iMouse.xy)) / iResolution.y;
    }
    vec3 ro = vec3(0.0, 0.0, 2.5); // camera
    //camera rotation
    ro.yz *= rot(mo.y * PI);
    ro.xz *= rot(-mo.x * TAU);

    const float fl = 1.5; // focal length
    float dist = dist_infin;

    //vec3 b1 = vec3(0.23529411764705882, 0.4235294117647059, 0.7725490196078432), b2 = vec3(0.3686274509803922, 0.5725490196078431, 0.8941176470588236);
    //vec3 bg = mix(b1, b2, vec3((1.0 - abs(fragCoord.x - iResolution.x / 2.) / iResolution.y * 2.) * fragCoord.y / iResolution.x));   
    vec3 bg = snowBgcol; //mix(b2, b1, fragCoord.y / iResolution.y);   
    //antialiasing
    vec3 tot = vec3(0.0);
    for(int m = 0; m < AA; m++) for(int n = 0; n < AA; n++) {
            vec2 o = vec2(float(m), float(n)) / float(AA) - 0.5;
            vec2 p = (-iResolution.xy + 2.0 * (fragCoord + o)) / iResolution.y;
            vec3 rd = GetRayDir(p, ro, vec3(0, 0., 0), fl); //ray direction
            vec3 col = bg * bg; // background  
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

                //reflection

                vec3 psk = reflect(rd, nor);
                vec3 c2 = csky(psk);

                col = calccolor(col, col, -rd, light, light2, nor);
                col = mix(col, c2, resReflect);

                //col += c2*0.1;

            }
            //==========================raymatch=============================
            tot += col;
        }
    tot = sqrt(tot) / float(AA);
    //tot = pow(tot, vec3(0.7)) / float(AA);
    //antialiasing
    fragColor = vec4(tot, 1.0);
}


/////=====================================================================================
void main() {
    vec4 fragColor = vec4(0);
    mainImage(fragColor, gl_FragCoord.xy);
    gl_FragColor = fragColor;
}