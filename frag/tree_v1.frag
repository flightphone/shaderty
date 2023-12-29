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
#define rot(f) mat2(cos(f), -sin(f), sin(f), cos(f))

const float dist_infin = 5.0;
#define nn 128
const float eps = 0.001;

vec3 sdfColor;
vec3 resColor;
vec3 col1 = vec3(0.13725490196078433, 0.4823529411764706, 0.28627450980392155);
vec3 col2 = vec3(0.3607843137254902, 0.16470588235294117, 0.027450980392156862);

float spine(vec3 p, float l, float x, float z, float th) {
    vec3 s = vec3(sin(z) * cos(x), sin(z) * sin(x), cos(z));
    float d = clamp(dot(p, s), 0., l);
    float r1 = 0.03, r2 = 0.005, r = r1 + (r2 - r1) * d / l;
    //r *= th;
    return length(p - d * s) - r;

}

float getlon(float lon, float n, float shift) {
    lon = mod(lon - shift, TAU);
    float dlon = TAU / n, lon1 = floor(lon / dlon) * dlon;
    if((lon - lon1) > dlon / 2.)
        lon1 += dlon;
    return lon1 + shift;
}

float branch(vec3 p, float l, float r, float ls, float lss, float fi, float n, float th, int prec) {
    float z = clamp(p.z, 0., l);// r2 = r - r * z / l*0.6;
    float d = length(p - vec3(0, 0., z)) - r;
    sdfColor = col2;
    
    float lon = mod(atan(p.y, p.x), TAU), dlon = TAU / n;
    float j = floor(z / lss);
    float h1 = j * lss, h2 = h1 + lss, shift1 = mod(j, 2.) * dlon / 2., shift2 = mod((j + 1.), 2.) * dlon / 2.;
    float h3 = h1 - lss, h4 = h1 - 2.*lss, shift3 = mod(j - 1., 2.) * dlon / 2., shift4 = mod((j - 2.), 2.) * dlon / 2.;

    float lon1 = getlon(lon, n, shift1), lon2 = getlon(lon, n, shift2);
    float lon3 = getlon(lon, n, shift3), lon4 = getlon(lon, n, shift4);

    float d2 = dist_infin;
    d2 = spine(p - vec3(r * cos(lon1), r * sin(lon1), h1), ls, lon1, fi, th);
    //d2 = laspine(p, r, h1, lon1, ls, fi, th, l);
    if(d2 < d) {
        d = d2;
        sdfColor = col1;
    }
    
    //need for fi > PI/6.
    if (prec == -1)
    if(h2 < l) {
        d2 = spine(p - vec3(r * cos(lon2), r * sin(lon2), h2), ls, lon2, fi, th);
        //d2 = laspine(p, r, h2, lon2, ls, fi, th, l);
        if(d2 < d) {
            d = d2;
            sdfColor = col1;
        }
    }
    
    if (prec >= 1)
    if(h3 >= 0.) {
            d2 = spine(p - vec3(r * cos(lon3), r * sin(lon3), h3), ls, lon3, fi, th);
            //d2 = laspine(p, r, h3, lon3, ls, fi, th, l);
            if(d2 < d) {
                d = d2;
                sdfColor = col1;
            }
    }
    
    if (prec >= 2)
    if(h4 >= 0.) {
            d2 = spine(p - vec3(r * cos(lon4), r * sin(lon4), h4), ls, lon4, fi, th);
            if(d2 < d) {
                d = d2;
                sdfColor = col1;
            }
    }
        
    
    return d;
}


float iTT = 0.;

float map(vec3 p) {
    /*
    p = p.xzy;
    p += vec3(.0, .0, 1.3);
    float td = branch(p, 2., 0.1, 1., 0.3, PI/6., 20., 0.3, 2);
    resColor = sdfColor;
    return td;
    */

    float l = 4.;
    p.xy *= rot(PI/3. + iTT);
    p.yz *= rot(PI/3.+ iTT); 
    //float f = PI /4. * p.z / l;
    //p.yz *= rot(f);
    p += vec3(0.0, 0.0, 2.);
    float d = (length(p - vec3(.5, .5, 2.)) - 0.5)*0.8;
    vec3 col = vec3(1.0, 0.2901960784313726, 0.17647058823529413);

    float d2 = dist_infin;
    d2 = branch(p, l, 0.025, 0.5, 0.3, PI / 6., 20., 1., 1);
    if(d2 < d) {
        col = sdfColor;
        d = d2;
    }

    p -= vec3(0., 0., l / 2.);
    vec3 p2 = p;
    p2.xz *= rot(PI / 4.);
    d2 = branch(p2, l / 2., 0.025, 0.5, 0.3, PI / 6., 20., 1., 1);
    if(d2 < d) {
        col = sdfColor;
        d = d2;
    }

    p2 = p;
    p2.xz *= rot(PI / 6.);
    d2 = branch(p2, l / 2., 0.025, 0.5, 0.3, PI / 6., 20., 1., 1);
    if(d2 < d) {
        col = sdfColor;
        d = d2;
    }

    p2 = p;
    p2.xz *= rot(-PI / 4.);
    d2 = branch(p2, l / 2., 0.025, 0.5, 0.3, PI / 6., 20., 1., 1);
    if(d2 < d) {
        col = sdfColor;
        d = d2;
    }
    /*
    p2 = p;
    p2.xz *= rot(-PI *0.7);
    d2 = branch(p2, l / 2., 0.025, 0.5, 0.3, PI / 6., 20., 1., 1);
    if(d2 < d) {
        col = sdfColor;
        d = d2;
    }
    */
    sdfColor = col;
    resColor = sdfColor;
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
/*
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
*/

vec3 ccolor(vec3 col, vec3 rd, vec3 light, vec3 nor) {
    vec3 R = reflect (light, nor);
    float shininess=20.0;
    float specular    =  pow(max(dot(R, rd), 0.), shininess);
    
    float difu = dot(nor, light);
    //col = col*clamp(difu, 0.3, 1.0) + vec3(.5)*specular*specular;
    col = col*(col*clamp(difu, 0., 1.0) + 0.3) + vec3(.5)*specular*specular;
    return col;
}

vec3 calccolor(vec3 col, vec3 backcol, vec3 rd, vec3 light1, vec3 light2, vec3 nor) {
    float difu1 = dot(nor, light1);
    float difu2 = dot(nor, light2);
    float difu = max(difu1, difu2);
    

    vec3 R1 = reflect (light1, nor);
    vec3 R2 = reflect (light2, nor);
    float shininess=10.0;
    float specular1    =  pow(max(dot(R1, rd), 0.), shininess);
    float specular2    =  pow(max(dot(R2, rd), 0.), shininess);
    float specular = max(specular1, specular2);

    
    //col = col*clamp(difu, 0.3, 1.0) + vec3(.5)*specular*specular;
    col = col*(col*clamp(difu, 0., 1.0) + 0.3) + vec3(.5)*specular*specular;

    return col;
}

//https://www.shadertoy.com/view/dlKBRc
vec3 lightingv3(vec3 lightColor, vec3 rd, vec3 L, vec3 normal) 
{   
    vec3 V = rd;
    vec3 N = normal;
    vec3 R = reflect (-L, N);
    float shadow = 1.;
    float occ = 1.;
    float Ka = 0.5;
    vec3 ambient = Ka + Ka * dot(normal, vec3(0., 1., 0.))*lightColor;
    ambient*=0.5;
    vec3 fresnel =  lightColor *  pow(clamp(1.0 + dot(rd, N), 0.0, 1.0), 2.0);;
    float diff= clamp(dot(N, L), 0., 1.0);
    vec3 diffuse =  lightColor * diff;
    float shininess=10.0;
    float specular    = pow(max(dot(R, V), 0.0), shininess);
    vec3 back = 0.5 * lightColor * clamp(dot(N, -L), 0.0, 1.0); // back
    vec3 colOut = occ*lightColor*(ambient+diffuse*shadow+.25 +back) + vec3(.5)*specular*specular;
    return colOut;
}
// https://knarkowicz.wordpress.com/2016/01/06/aces-filmic-tone-mapping-curve/
vec3 ACESFilm(vec3 x){
    return clamp((x * (2.51 * x + 0.03)) / (x * (2.43 * x + 0.59) + 0.14), 0.0, 1.0);
}

// https://iquilezles.org/articles/nvscene2008/rwwtt.pdf
float calcAO( in vec3 pos, in vec3 nor )
{
	float occ = 0.0;
    float sca = 1.0;
    for( int i=0; i<5; i++ )
    {
        float h = 0.01 + 0.12*float(i)/4.0;
        float d = map( pos + h*nor );
        occ += (h-d)*sca;
        sca *= 0.95;
        if( occ>0.35 ) break;
    }
    return clamp( 1.0 - 3.0*occ, 0.0, 1.0 ) * (0.5+0.5*nor.y);
}

//https://www.shadertoy.com/view/Xds3zN , 555 row
vec3 iqcolor(vec3 col, vec3 rd, vec3 light, vec3 nor, vec3 pos)
{
    float ks = 1.0; //0.4;
    // lighting
    float occ = 1.;//calcAO( pos, nor );
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

        // sky
        {
            float dif = sqrt(clamp( 0.5+0.5*nor.y, 0.0, 1.0 ));
                  dif *= occ;
            float spe = smoothstep( -0.2, 0.2, ref.y );
                  spe *= dif;
                  spe *= 0.04+0.96*pow(clamp(1.0+dot(nor,rd),0.0,1.0), 5.0 );
          //if( spe>0.001 )
                  //spe *= calcSoftshadow( pos, ref, 0.02, 2.5 );
            lin += col*0.60*dif*vec3(0.40,0.60,1.15);
            lin +=     2.00*spe*vec3(0.40,0.60,1.30)*ks;
        }

        // back
        {
        	float dif = clamp( dot( nor, normalize(vec3(0.5,0.0,0.6))), 0.0, 1.0 )*clamp( 1.0-pos.y,0.0,1.0);
                  dif *= occ;
        	lin += col*0.55*dif*vec3(0.25,0.25,0.25);
        }
        // sss
        {
            float dif = pow(clamp(1.0+dot(nor,rd),0.0,1.0),2.0);
                  dif *= occ;
        	lin += col*0.25*dif*vec3(1.00,1.00,1.00);
        }
        return lin;
}
void mainImage(out vec4 fragColor, in vec2 fragCoord) {
    vec3 light = normalize(vec3(0.0, 1.0, -2.5)); //light
    vec3 light2 = normalize(vec3(0.0, -1.0, 2.5)); //light
    vec2 mo = vec2(0.0, 0.0);
    iTT = iTime;
    //if  (iMouse.z > 0.0)
    {
        mo = (-iResolution.xy + 2.0 * (iMouse.xy)) / iResolution.y;
        iTT = 0.;
    }
    vec3 ro = vec3(0.0, 0.0, 2.5); // camera
    //camera rotation
    ro.yz *= rot(mo.y * PI);
    ro.xz *= rot(-mo.x * TAU);

    const float fl = 1.5; // focal length
    float dist = dist_infin;

    vec3 b1 = vec3(0.23529411764705882, 0.4235294117647059, 0.7725490196078432), b2 = vec3(0.3686274509803922, 0.5725490196078431, 0.8941176470588236);
    vec3 bg = mix(b2, b1*b1, fragCoord.y / iResolution.y);     
    //bg = vec3(1.);
    //antialiasing
    vec3 tot = vec3(0.0);
    for(int m = 0; m < AA; m++) for(int n = 0; n < AA; n++) {
            vec2 o = vec2(float(m), float(n)) / float(AA) - 0.5;
            vec2 p = (-iResolution.xy + 2.0 * (fragCoord + o)) / iResolution.y;
            vec3 rd = GetRayDir(p, ro, vec3(0, 0., 0), fl); //ray direction
            vec3 col = bg; // background  
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
                col = calccolor(col, col, rd, light, light2, nor);
                //col = ccolor(col, rd, light, nor);
                //col = lightingv3(col, -rd, light, nor); 
                //col = iqcolor(col, rd, light, nor, pos);

            }
            //==========================raymatch=============================
            tot += col;
        }
    tot = tot / float(AA);
    //tot = pow(tot, vec3(0.7)) / float(AA);
    tot = ACESFilm(tot);
    //antialiasing
    fragColor = vec4(tot, 1.0);
}
/////=====================================================================================
void main() {
    vec4 fragColor = vec4(0);
    mainImage(fragColor, gl_FragCoord.xy);
    gl_FragColor = fragColor;
}