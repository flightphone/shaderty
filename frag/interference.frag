//1. Go fullscreen
//2. Take drugs now

//iq noise fn
float hash( float n )
{
    return fract(sin(n)*43758.5453);
}
float noise( in vec3 x )
{
    vec3 p = floor(x);
    vec3 f = fract(x);

    f = f*f*(3.0-2.0*f);
    float n = p.x + p.y*57.0 + 113.0*p.z;
    return mix(mix(mix( hash(n+  0.0), hash(n+  1.0),f.x),
                   mix( hash(n+ 57.0), hash(n+ 58.0),f.x),f.y),
               mix(mix( hash(n+113.0), hash(n+114.0),f.x),
                   mix( hash(n+170.0), hash(n+171.0),f.x),f.y),f.z);
}

//x3
vec3 noise3( in vec3 x)
{
	return vec3( noise(x+vec3(123.456,.567,.37)),
				noise(x+vec3(.11,47.43,19.17)),
				noise(x) );
}

//http://dept-info.labri.fr/~schlick/DOC/gem2.ps.gz
float bias(float x, float b) {
	return  x/((1./b-2.)*(1.-x)+1.);
}

float gain(float x, float g) {
	float t = (1./g-2.)*(1.-(2.*x));	
	return x<0.5 ? (x/(t+1.)) : (t-x)/(t-1.);
}


mat3 rotation(float angle, vec3 axis)
{
    float s = sin(-angle);
    float c = cos(-angle);
    float oc = 1.0 - c;
	vec3 sa = axis * s;
	vec3 oca = axis * oc;
    return mat3(	
		oca.x * axis + vec3(	c,	-sa.z,	sa.y),
		oca.y * axis + vec3( sa.z,	c,		-sa.x),		
		oca.z * axis + vec3(-sa.y,	sa.x,	c));	
}

vec3 fbm(vec3 x, float H, float L, int oc)
{
	vec3 v = vec3(0);
	float f = 1.;
	for (int i=0; i<10; i++)
	{
		if (i >= oc) break;
		float w = pow(f,-H);
		v += noise3(x)*w;
		x *= L;
		f *= L;
	}
	return v;
}

vec3 smf(vec3 x, float H, float L, int oc, float off)
{
	vec3 v = vec3(1);
	float f = 1.;
	for (int i=0; i<10; i++)
	{
		if (i >= oc) break;
		v *= off + f*(noise3(x)*2.-1.);
		f *= H;
		x *= L;
	}
	return v;	
}

vec2 Q(float a, float b, float c)
{
	float d = b*b-4.0*a*c;
	if (d < 0.0) return vec2(-1.,-1.);
	d=sqrt(d);	
	float oo2a = 0.5/a;
	vec2 i = vec2(-b-d,-b+d)*oo2a;
//	return vec2( min(i.x,i.y), max(i.x,i.y) );
	return i;
}

vec2 RayEllipsoid(vec3 P, vec3 V, vec3 A)
{
	A*=A;
	vec3 VP = V*P;
	P *= P;
	V *= V;
	
	vec3 S=A.yzx*A.zxy;
		
	float a = dot(V,S);
	float b = 2. * dot(VP,S);
	float c = dot(P,S) - A.x*A.y*A.z;
				
	return Q(a,b,c);
}

//x^2   y^2   z^2   
//--- + --- + --- = 1
//a^2   b^2   c^2   

//grad 
// 2x
// --- 
// a^2

vec3 NormalEllipsoid(vec3 P, vec3 A)
{
	return -normalize(P/(A*A));
}

#define pi 3.1415927

vec3 RotY(vec3 p, float t) {
	float c = cos(t); float s = sin(t);
	return vec3(p.x*c+p.z*s,
				p.y,
				-p.x*s+p.z*c);
}

void MakeViewRay(in vec2 fragCoord, out vec3 eye, out vec3 ray)
{
	vec2 ooR = 1./iResolution.xy;
    vec2 q = fragCoord.xy * ooR;
//	q.x = mod(q.x*3.,1.);	//tile it for split views!
//	q.y = mod(q.y*2.,1.);	//tile it for split views!
    vec2 p =  2.*q -1.;
    p.x *= iResolution.x * ooR.y;
	
    vec3 lookAt = vec3(0.,0.,0.);
	float t = 24.4 + iTime*0.01; //iTime; //sin(iTime*.5)*.5;
	eye = vec3(-2.5,0,2.5);
	eye = RotY(eye,t*-0.5*pi);	
	
    // camera frame
    vec3 fo = normalize(lookAt-eye);
    vec3 ri = normalize(vec3(fo.z, 0., -fo.x ));
    vec3 up = normalize(cross(fo,ri));
     
    float fov = .25;
	
    ray = normalize(fo + fov*p.x*ri + fov*p.y*up);
}

vec3 Colorize(vec3 p)
{
	float time = iTime * 1.276;
	
	float slow = time*0.02;
	
	vec3 axis = 4. * fbm(p+vec3(0.,slow,0.), 0.5, 2., 8);				//random fbm axis of rotation
	
	vec3 colorVec = 0.5 * 5. * fbm(p*0.3,0.5,2.,7);		//random base color
	p += colorVec;
	
//	float mag = 4e5;	//published, rather garish?
	float mag = 0.75e5; //still clips a bit
//	mag = mag * (1.+sin(2.*3.1415927*ts)*0.75);
	vec3 colorMod = mag * smf(p,0.7,2.,8,.2);			//multifractal saturation
	colorVec += colorMod;
	
	colorVec = rotation(3.*length(axis)+slow*0.,normalize(axis))*colorVec;
	
	colorVec *= 0.075;
			
	colorVec = colorVec / (1. + length(colorVec));	//tone it all down a bit
	
	return colorVec;
}

float shlick(vec3 N, vec3 V)
{
	float f = dot(-V,N);
	f = 1.0-f;	
	float ff = f;
	f *= f;		//2
//	f *= f;		//4
//	f *= ff;	//5
	float r0 = 0.025;
	f = r0 + (1.0-r0)*f;
	return f;
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
	vec3 mainViewP, mainViewD;
	MakeViewRay(fragCoord,mainViewP, mainViewD);
	
	vec3 bg = texture(iChannel0,mainViewD).xyz;
	bg = pow(bg,vec3(2.2));
		
	vec3 c = bg;
		
//	float s = sin(iTime)*0.025 + 1.;
//	s * t * u = 1;
//  tu = 1/s
//  t=u=sqrt(1/s);
//	float tu = sqrt(1./s);
//	A *= vec3(tu,s,tu);
	
	for (int inst=0; inst<5; inst++) { 
	vec3 viewP = mainViewP;
	vec3 viewD = mainViewD;
		
	vec3 A = vec3(.3,.4,.5)*2. * (1.-float(inst)*0.3);
		
		
	vec2 uv = fragCoord.xy /iResolution.xy;
	uv.x *= iResolution.x / iResolution.y;
	
	vec3 ax=noise3(vec3(uv*2.5,iTime*0.25+A.x));
	viewD=rotation(length(ax)*(0.05),normalize(ax))*viewD;		//wibble the view ray to change bubble shape
	
//	viewP += (noise3(vec3(iTime*A.y,A.x,A.z))*2.-1.)*vec3(1.,0.6,1.1);		//float around
		
	viewP += smf(vec3(iTime*A.y*.5,A.x,A.z), 0.51, 2.13, 2, 0.2)*vec3(1.,0.6,1.1);		//float around
			
	viewP.yz -= (mod(iTime*0.05 + float(inst)*.37,1.)*2.-1.)*2.5;		//float up
	
	vec2 tt = RayEllipsoid(viewP, viewD, A);
	
	if (tt.y > tt.x)
	{
		tt = tt.yx; //oops was blending over so go back to front.
		
		for (int i=0; i<2; i++)
		{
			vec3 hit_p = tt[i]*viewD + viewP;
			
			vec3 n = NormalEllipsoid(hit_p, A);
				
			float facing = dot(n,-viewD) > 0. ? .25 : .025; //max(,0.);
			
			float e = 1.- abs(dot(n,viewD));
			
			hit_p.y += iTime*0.1; 
			vec3 soap_col = Colorize(RotY(n*0.5,A.x+iTime*0.2)-hit_p*0.1);
			
			vec3 soap = soap_col * (e + .5*length(soap_col));
			
			
			vec3 r = reflect(-viewD,n);			//reflection
			r = texture(iChannel0,r).xyz;
			r = pow(r,vec3(2.2));
			float f = shlick(n, viewD);			//fresnel
			soap += f*r*facing;
			
			soap *= 1.-pow(e,10.);	//antialias
			//soap *= 1.-f;
			
//			c += soap;
			c = mix(c,soap,0.5*(length(soap)+0.5))+soap*.5;
		}
	}
//	c += tt.y - tt.x; 	
	
	}
	
	c = pow(c,vec3(1./2.2));
	fragColor = vec4(c,1.0);
	
}