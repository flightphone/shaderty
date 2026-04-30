//malachite texture
//malachite, 2d, noise, FBM

#define PI  3.14159265359
#define TAU 6.28318530718

float hash(vec2 p) {
	return fract(1e4 * sin(117.0 * p.x + p.y * 0.1) * (0.1 + abs(sin(p.y * 133.0 + p.x))));
} 

float noise(vec2 x) {
	vec2 i = floor(x);
	vec2 f = fract(x);
	float a = hash(i);
	float b = hash(i + vec2(1.0, 0.0));
	float c = hash(i + vec2(0.0, 1.0));
	float d = hash(i + vec2(1.0, 1.0));

	vec2 u = f * f * (3.0 - 2.0 * f);
	return mix(a, b, u.x) + (c - a) * u.y * (1.0 - u.x) + (d - b) * u.x * u.y;
}

vec3 get_palette(float t) {
    vec3 a = vec3(0.0, 0.2, 0.05);
    vec3 b = vec3(0.1, 0.3, 0.1);
    vec3 c = vec3(1.0, 1.0, 1.0);
    vec3 d = vec3(0.5, 0.2, 0.25);
    return a + b * cos(6.28318 * (c * t + d));
}

void mainImage(out vec4 fragColor, in vec2 fragCoord) {
    vec2 uv = fragCoord / iResolution.y;
    vec2 p = uv*2.; 
    p.x +=  iTime*0.3;

    float angle = PI / 7.0;
    mat2 rot = mat2(cos(angle), -sin(angle), sin(angle), cos(angle));
    vec2 shift = vec2(3.0, 2.5 );

    // Fractal Noise
    float val = 0.0;
    vec2 p_t = p;
    
    p_t = p_t * rot + shift;
    val += 1.8 * noise(2.0 * p_t);
    
    p_t = p_t * rot + shift;
    val += 0.5 * noise(5.0 * p_t);
    
    p_t = p_t * rot + shift;
    val += 0.05 * noise(20.0 * p_t);


    // sigmoid function parameters
    float H = 3.0 * TAU - 0.3;      //height
    float threshold = 0.9;          //point of sharp rise
    float delta_u = 0.5;            //width
    float k = 3.0 / delta_u;
    
    // transform
    float plateau = (H / 2.0) * (tanh(k * (val - threshold)) + 1.0);

    // layers
    float layers = (sin(plateau) + 1.0) / 2.0;
    vec3 col = get_palette(layers);

    // three white lines
    vec3 colline = vec3(0.4, 0.9, 0.4);
    float layers2 = (sin(log(plateau + 1.0)) + 1.0) / 2.0;
    
    
    float h = 0.001; 
    float eps = 0.002;
    
    float s1 = 1.0 - smoothstep(h, h + eps, abs(layers2 - 1.0));
    float s2 = 1.0 - smoothstep(4.*h, 4.*h + 4.*eps, abs(layers2 - 0.98));
    float s3 = 1.0 - smoothstep(4.*h, 4.*h + 3.*eps, abs(layers2 - 0.55));
    
    col = mix(col, colline, s1);
    col = mix(col, colline, s2);
    col = mix(col, colline, s3);

    // Corn
    float grain = hash(p_t * 10.0 ); 
    col += (grain - 0.5) * 0.05;

    fragColor = vec4(col, 1.0);
}
