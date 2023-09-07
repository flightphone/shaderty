#ifdef GL_ES
precision mediump float;
#endif

uniform vec2 u_resolution;
uniform vec2 u_mouse;
uniform float u_time;
uniform float u_stop;
uniform float u_stop_time;
uniform float u_reset;

uniform float u_res0;
uniform float u_res1;
uniform float u_res2;
uniform float u_res3;

uniform sampler2D u_tex0;
uniform sampler2D u_tex1;
uniform sampler2D u_tex2;
uniform sampler2D u_tex3;
uniform sampler2D u_tex4;
uniform sampler2D u_tex5;
uniform sampler2D u_tex6;
uniform sampler2D u_tex7;
uniform sampler2D u_tex8;
uniform sampler2D u_tex9;
uniform sampler2D u_tex10;
uniform sampler2D u_tex11;
uniform sampler2D u_tex12;
uniform sampler2D u_tex13;



#define PI 3.14159265359
const float sp = 1.5;

const float num_samples = 256.0;
const float tanimate =  PI/sp;

const float hh = 1.0;
const float kstep = 0.2;
const float ksum = 10.0;
const float eatt = 10.0;


vec3 attenuation(vec3 c, float d)
{
    return c/exp(eatt*d);
}


vec3 radial_blur_filter(in vec2 origin, in vec2 point)
{
    vec4 fon = texture2D(u_tex13, point);
    vec3 col = fon.rgb * fon.a;
    vec3 light = vec3(0.0);
    vec2 v = point - origin;
    float l = length(v);
    v = normalize(v); 
    float dz =  hh / num_samples * ksum;
    float l1 = l/(1.0 + hh);
    float al = (l-l1);
    for(float s = 0.0; s < num_samples; s++)
    {
        float l2 = l1 + al*s/num_samples;
        vec2 p = origin + v*l2;
        vec3 c = texture2D(u_tex13, p).rgb;
        float pst = smoothstep(kstep, 1.0, c.r);
        c *= pst;

        float fi = atan(l2, 1.0);
        float d = (l - l2)/sin(fi);
        c = attenuation(c, d);
        light += c;
    }
    light *= dz;
    float pst = smoothstep(0., 1.0, col.r);  
	light*=(1.0 - pst);
    light*= abs(sin(u_time*sp)); 
    col += light;
    return col;
}

void main() {
    float mih = 0.325;
    float mah = 1.0 - mih; 
    float xmin = 0.203;
    float xmax = 1.0 - xmin;
    
    vec2 pt = gl_FragCoord.xy/u_resolution;
    if (u_reset == 1.0)
    {
        vec4 fon = texture2D(u_tex13, pt);
        fon.rgb = fon.rgb * fon.a;
        gl_FragColor = fon;
        return;
    }


    float t = u_time - tanimate;

    //======================================init lights=================================
    if (t < 0.0)
    {
        vec4 col = texture2D(u_tex13, pt);
        vec2 mouse = vec2(0.5 - cos(u_time*sp),0.3);
        col.rgb = radial_blur_filter(mouse, pt);
        gl_FragColor = col;    
        return;
    }
    //======================================init lights=================================


    float tstop = u_stop_time;
    vec2 pto = pt;
    vec4 fon = texture2D(u_tex10, pt);

    

    float apst = 1.0 - smoothstep(0.1, 1.0, fon.a); 
    fon.rgb = mix(fon.rgb, vec3(0.0), apst);
    
    vec4 fono = fon;
    vec4 col = fon;
    
    
    

    
    vec2 c2 = vec2(0.4253393665158371, 0.5);
    float d1 = 0.14675;
    
    float sht = 2.0;
    float sht2 = 1.5;
    float ra = 0.6;
    

    if (pt.y > mih && pt.y < mah && pt.x > xmin && pt.x < xmax)
    {
        //deformation x
        float l = length(pt - 0.5);
        float z = pow((ra*ra - l*l), 0.5);
        float fi = atan(pt.x - 0.5, z);
        float aa = ra* sin(fi);
        pt.x = 0.5 + aa;
        col.rgb = mix(col.rgb, vec3(0.0), smoothstep(xmin-0.001, xmin, pt.x) - smoothstep(xmin, xmin+0.001, pt.x));
        
    }
    //
    
    if (pt.y > mih && pt.y < mah && pt.x > xmin && pt.x < xmax)
    {
        //recalc pt

        //shift
        float pos = 4.0*(pt.x - xmin)/(xmax - xmin);
        float sh = floor(pos);
        float dx = d1*(1.0 - sh);
        
        float fonpos = (pt.y-mih)/(mah-mih)*2.0/11.0;
        //if (t >= 0.0)
        //{
                vec4 col2 =vec4(0.0);
                t -= (3.0 - sh)*sht;
                if (t < 0.0)
                    t = 0.0;

                if (u_stop > 0.0)
                    tstop += sh*sht2;

                     

                float k = clamp(log(t), 0.1, 4.0);
                //k = 1.0;
                float num = mod(t*k, 10.0);
                fonpos = fract(num*0.1 + (pt.y-mih)/(mah-mih)*0.1);
                
                
                
                
                
                float f = fract(num);
                num = floor(num);
                
                if (u_stop > 0.0 && u_time > tstop)
                {
                    if (sh == 0.0) 
                        num = u_res0;
                    if (sh == 1.0) 
                        num = u_res1;    
                    if (sh == 2.0) 
                        num = u_res2;    
                    if (sh == 3.0) 
                        num = u_res3;        
                    f = 0.0;
                    fonpos = fract(num*0.1 + (pt.y-mih)/(mah-mih)*2.0/11.0);
                }
                fon = texture2D(u_tex11, vec2(fonpos, fract(pos)));
                

                
                vec2 shift = vec2(dx, (mah-mih)*f);
                vec2 pt1 = pt + shift;

                if (num == 0.0)
                    col = texture2D(u_tex0, pt1);
                if (num == 1.0)
                    col = texture2D(u_tex1, pt1);
                if (num == 2.0)
                    col = texture2D(u_tex2, pt1);
                if (num == 3.0)
                    col = texture2D(u_tex3, pt1);
                if (num == 4.0)
                    col = texture2D(u_tex4, pt1);
                if (num == 5.0)
                    col = texture2D(u_tex5, pt1);
                if (num == 6.0)
                    col = texture2D(u_tex6, pt1);
                if (num == 7.0)
                    col = texture2D(u_tex7, pt1);
                if (num == 8.0)
                    col = texture2D(u_tex8, pt1);
                if (num == 9.0)
                    col = texture2D(u_tex9, pt1);


                shift = vec2(dx, -(mah-mih)*(1.0-f));
                pt1 = pt + shift;

                if (num == 0.0)
                    col2 = texture2D(u_tex1, pt1);
                if (num == 1.0)
                    col2 = texture2D(u_tex2, pt1);
                if (num == 2.0)
                    col2 = texture2D(u_tex3, pt1);
                if (num == 3.0)
                    col2 = texture2D(u_tex4, pt1);
                if (num == 4.0)
                    col2 = texture2D(u_tex5, pt1);
                if (num == 5.0)
                    col2 = texture2D(u_tex6, pt1);
                if (num == 6.0)
                    col2 = texture2D(u_tex7, pt1);
                if (num == 7.0)
                    col2 = texture2D(u_tex8, pt1);
                if (num == 8.0)
                    col2 = texture2D(u_tex9, pt1);
                if (num == 9.0)
                    col2 = texture2D(u_tex0, pt1);

                
                //fon = fono;

                float pst = smoothstep(0.0, 0.5, abs(pt.y - 0.5));
                col.rgb = pow(col.rgb, vec3(15.0*pst));
                col2.rgb = pow(col2.rgb, vec3(15.0*pst));
               

                float t2 = 0.0;
                if (u_stop > 0.0 && u_time > tstop)
                {
                    t2 = (u_time - tstop)*3.0;
                    if (t2 < PI)
                        col.rgb = pow(col.rgb, vec3(abs(cos(t2))));
                }

                for (float i = 0.0; i < 4.0; i++)
                {
                    dx = d1*(1.0 - i);
                    tstop = u_stop_time + i*sht2;
                    if (u_stop > 0.0 && u_time > tstop)
                        t2 = (u_time - tstop)*3.0;

                    if (u_stop > 0.0 && t2 < PI) 
                    {
                        vec2 ll =  vec2(c2.x -dx, c2.y);
                        float pst = clamp(sin(t2)/exp(15.0*length(pt-ll)), 0.0, 0.8);
                        fon.rgb = mix(fon.rgb, vec3(0.95, 0.95, 0.95), pst*fon.a);
                    }
                }
                
                col = mix(fon, col, col.a);
                col = mix(col, col2, col2.a);

                //light texture
                vec4 mask = texture2D(u_tex12, pto);
                col.rgb = mix(col.rgb, mask.rgb, 2.0*mask.a);    

                //edges
                
                pst = 1.0 - smoothstep(0.0, 0.01, pt.x - xmin);
                col.rgb = mix(col.rgb, vec3(0.0), pst);
                
                fon = texture2D(u_tex10, pto);
                pst = smoothstep(xmax-0.01, xmax, pt.x);
                col.rgb = mix(col.rgb, fon.rgb, pst);
                
                pst = smoothstep(mah-0.04, mah, pt.y);
                col.rgb = mix(col.rgb, fon.rgb, pst);
                pst = smoothstep(mih+0.04, mih, pt.y);
                col.rgb = mix(col.rgb, fon.rgb, pst);
                
        //}
        
    }
    
    gl_FragColor = col;     
}
