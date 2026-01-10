#iChannel0 "file://../img/spree_bank.webp"
#iChannel1 "file://../img/stoun.jpg"


/*

    Rectangle Grid Cayley Transform
    -------------------------------
    
    Among other things, MLA posts a lot of really nice non-Euclidean math 
    examples. I really appreciate it when someone on Shadertoy demystifies a
    particular subject by posting useable code, and MLA does a lot of that. The 
    other day, he posted a really cool-looking tiling of the unit disc that 
    caught my eye, and I really wanted to make one, so this is a rough 
    interpretation.
    
    In essence, Cayley transforms can be used to map a Euclidean rectangular 
    grid to bands of circles inside a unit disc. It's similar to way in which 
    you'd polar map a rectangular grid, except that the tile dimensions vary in 
    such a way that the horizon never exceeds the boundary of a disc. Obviously, 
    there are tiles outside the disc in this example, but they're just 
    reflections across the boundary.
    
    The imagery produced using this process is fairly common. However, MLA took 
    the extra step of flipping the Y-values across the vertical center, which 
    produced a really interesting looking wavy pattern. I don't know if this is 
    a regular thing to do, but it was new to me.
    
    It was necessary to rewrite part of the original code to suit this particular 
    example, but it's essentially the same process with window dressing applied. 
    For better or worse, I trimmed down the transform logic and changed the way 
    the texture coordinates are calculated. I also added in scaling, equal line 
    spacing and natural log functionality.
    
    By the way, there are a few "define" options below for anyone interested. In
    fact, the options and extra lighting are responsible for a large portion of
    the character count. The pattern itself doesn't require a lot of code at all.
    
    Since horocycle lines have constant curvature and this pattern (although 
    flipped), would fall under that in a piecewise sense, it would be possible to 
    perform a raymarch traversal... I'd like to do that at some stage, just to 
    show that it can be done, but I'll wait until I have more coding energy. :)
    
    

    Based on:
    
    // One of my favorite non-Euclidean examples in a while.
    Horocyclic Tiling - MLA
    https://www.shadertoy.com/view/X3fBD4
    
*/

/////////////////////////////////

// Constant defines.

// PI and 2PI.
const float PI = 3.14159265358979;
const float TAU = PI*2.;

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

/////////////////////////////////
// Variable defines.

// Only color the inner disc, not the outside.
#define INNER_COL

// Color scheme -- Red: 0, Blue: 1, Green: 2, Greyscale: 3.
#define COLOR 1

// Add central holes.
#define HOLES

// Rectangle grid cell scale.
vec2 gSc = vec2(2., 1)*.57735; // Just a scale with no meaning.
//vec2 gSc = vec2(1, LOG(2.)); // Binary split.

bool pZoom = true; // Periodic zoom.
bool halfPlane = false; // Half-plane.
bool gFlip = true; // Y-flipping across the vertical center.


// The type of logarithm -- I wanted to show that it
// works with natural logarithms too.
#define LOGTYPE 0
#if LOGTYPE == 0
// Natural log.
#define EXP exp
#define LOG log
#else
// Powers of 2.
#define EXP exp2
#define LOG log2
#endif

/////////////////////////////////



// Standard 2D rotation formula.
mat2 rot2(in float a){ float c = cos(a), s = sin(a); return mat2(c, s, -s, c); }

// Hash without Sine -- Dave Hoskins
// https://www.shadertoy.com/view/4djSRW
// 1 out, 2 in...
float hash21(vec2 p) {
 
    p = fract(p*vec2(123.34, 456.21));
    p += dot(p, p + 45.32);
    return fract(p.x*p.y);
}

// IQ's 2D box formula, wich rounding.
float sBox(in vec2 p, in vec2 b, in float rf){
  
  p = abs(p) - b + rf;
  return min(max(p.x, p.y), 0.) + length(max(p, 0.)) - rf;
}

// The Cayley transform, which is similar to a circle inversion,
// Mobius tranform, and so forth. It's relatively easy to see that
// it maps the Euclidean plane to a disc.
vec2 cayley(vec2 z) {

    z += vec2(0, 1);
    z *= 2./dot(z, z);
    z -= vec2(0, 1);
    return z;
}

 
// Mouse transform: Similar to the transform above.
vec2 mTrans(vec2 z){
    
    vec2 c = vec2(0);
    if (iMouse.z>0.) c = (2.*iMouse.xy - iResolution.xy)/iResolution.y;
    z -= c;
    z *= (dot(c, c) - 1.)/dot(z, z);
    z += c;
    if (!halfPlane || length(c)<1.) z /= -dot(z,z);
    return z;
    
} 
 
// Global texture coordinates.
vec2 gTUV;
vec2 gP;
 
vec3 distF(vec2 p){

    
    // Animation time.
    float tm = iTime/2.;
    
    // Halfplane model.
    if(halfPlane) p = cayley(p + vec2(0, 1));
     
    // Mouse translation.
    p = mTrans(p);
 
    // Flipping the pattern on one side of the vertical axis to
    // produce the wavy pattern.
    float flip;
    if(gFlip==true) flip = p.x<0.? -1. : 1.;
    else flip = 1.;
    
    // Mapping to the unit disc, or unmapping if using the halfplane model.
    p = cayley(p*vec2(1, flip));
 
    // Periodic zoom: Exp(Log(1)) = 1.
    if(pZoom) p /= EXP(flip*mod(tm/2., LOG(1.))); 
  
 
    //////////
    
    // Rectangle scale and ID.
    vec2 sc = gSc;
    vec2 ip;
    
    // Scaling back the Y-value to account for the Cayley transform, then
    // obtaining the row ID.
    p.y = flip*LOG(abs(p.y));
    ip.y = floor(p.y/sc.y);
     
    // Column translation speed and direction.
    float speed = hash21(vec2(7, ip.y)) + .5;
    float dir = hash21(vec2(3.5, ip.y))<.5? -1. : 1.;
    
    // Scaling the X-value according to the row number -- Rows further 
    // out descrease in size exponentially. We're also adding some animation, 
    // then obtaining the column ID.
    p.x *= EXP(-(floor(flip*p.y/sc.y) + .5)*sc.y);
    p.x += dir*speed*tm;//  + hash21(vec2(2, ip.y))*4.; // Animation.
    ip.x = floor(p.x/sc.x);
 
    // Saving the unwrapped coordinates for texturing purposes.
    gTUV = p;
    
    // Local grid cell coordinates.
    // Equivalent to: p = mod(p, sc) - sc/2.;
    p -= (ip + .5)*sc; 
    
    //ip.y = mod(ip.y, 8.); // Row mod, for patterns, etc.
 
    // Rendering a box with a hole in it.    
    float d = sBox(p, sc/2., sc.x*.1);
    
    #ifdef HOLES
    // The hole was originally present for practical reasons, namely, 
    // to check the X-stretch factor -- The hole needs to have constant 
    // curvature (non-elliptical.). However, I've kept it for
    // aesthetic purposes.
    d = max(d, -(length(p) - .05)); 
    #endif
    
    gP = p;
 
    // Distance and ID.
    return vec3(d, ip);

}


void mainImage(out vec4 fragColor, vec2 fragCoord) {
 
    vec2 p = (2.*fragCoord-iResolution.xy)/iResolution.y;
    vec2 oP = p;
    
    if(!halfPlane){
        // Rotating and reducing the size of the circle.
        //p *= rot2(-iTime/4.);
        p *= 1.1;
    } 
    
    // Disc and plane distances. Used for shading and so forth.
    float disc;
    if(halfPlane) disc = (p.y + 1.)/2.;
    else disc = length(p) - 1.;

    /* 
    // Rough gradient approximation. We already have the gradient, 
    // on account of the bump map samples, so it's not needed. For 3D 
    // applications, however, it would be handy.
    float dt;
    if(halfPlane) dt = 1./(p.y + 1.);
    else dt = 1.5/abs(dot(p, p) - 1.);
    */

    // The scene distance object.
    vec3 d3 = distF(p);
    
    // Distance, cell ID and texture coordinates.
    float d = d3.x; 
    vec2 ip = d3.yz;
    vec2 tuv = gTUV;
    
    // Nearby X and Y distance object samples. The samples are 
    // used for bump mapping and to determine the numeric distance
    // gradient, which can be used for things like even spacing
    // between blocks.
    float px = 1e-4;     
    vec3 dtX = distF(p + vec2(px, 0)); // Nearby X sample.
    vec3 dtY = distF(p + vec2(0, px)); // Nearby Y sample.
    float dt = length(vec2(dtX.x, dtY.x) - d)/px/sqrt(2.); // Gradient value.

   
    
    // Scaling all samples by the gradient.
    d /= dt;
    dtX.x /= dt;  
    dtY.x /= dt; 

    
    // Smoothing factor and edge width.
    float sf = 2./450.;
    float ew = .01;
    // Tapering the edge width at the disc edges, which looks a little cleaner.
    ew *= mix(.57, 1., smoothstep(0., .25, abs(disc)));

    // Applying the edge width.
    d += ew;
    dtX.x += ew;
    dtY.x += ew; 
    
    /* 
    // Beveling, if preferred.
    d = max(d, -.05);
    dtX.x = max(dtX.x, -.05);
    dtY.x = max(dtY.x, -.05);
    */
    
    

    // Using the rectangle ID to color the rectangle tiles.
    
    // Variable ID-based color.
    float range = hash21(ip + .15);
    vec3 pCol = .5 + .45*cos(TAU*range/8. + vec3(0, PI/2., PI)*.7 + .5);
    
    // Color alternate rows, either inside the disc, or on either side.
    #ifdef INNER_COL
    if((disc<0. || halfPlane) && mod(ip.y, 2.)==1.){
    #else
        if(mod(ip.y, 2.)==1.){
    #endif
        
        #if COLOR == 0
        // Red.
        pCol = mix(pow(pCol.xzy, vec3(1.5))*2., pCol, smoothstep(0., .25, abs(disc)));
        #elif COLOR == 1
        // Blue.
        pCol = mix(pCol.yzx*1.35, pCol.zyx, smoothstep(0., .3, abs(disc)));
        #elif COLOR == 2
        // Green.
        pCol = mix(pCol.zxy*pCol.zxy, pCol.yxz, smoothstep(0., .4, abs(disc)));
        #else
        // Grey. Virtually the same shade as the surrounds.
        pCol = vec3(range*.1 + .1); 
        #endif
    }
    else{
        // Greyscale
        pCol = vec3(range*.06 + .1); 
    }
     
     


    // Texturing.
    vec3 tx = texture(iChannel0, tuv/2.).xyz; tx *= tx;
    tx = smoothstep(0., .5, tx);
    if(disc>0.) pCol *= tx*2.;
    else pCol *= tx*2. + .25;
    
    
    // Hacky lighting setup: Normal, point light and unit direction vector.
    vec3 n = normalize(vec3((dtX.x - d)/px, (dtY.x - d)/px, -1));
    vec3 ld = normalize(vec3(4, 8, -2) - vec3(oP, 0)); //vec3(.5, 1, -.125)
    vec3 rd = normalize(vec3(oP, 1));
    
    // Greyscale texture value.
    //float gr = dot(tx, vec3(.299, .587, .114));
    // Diffuse lighting.
    float diff = max(dot(ld, n), 0.);
    //diff = pow(diff, 2. + gr*8.)*2.;
  	// Specular lighting.
	float spec = pow(max(dot(reflect(ld, n), rd ), 0.), 64.); 
 
    
    // Lamest texture based bump mapping ever. :D
    n.z -= pow(dot(tx, vec3(.299, .587, .114)), 4.)*.25;
    n = normalize(n);
    
    
    // Applying the distance coloring to the canvas.
    vec3 col = mix(vec3(0), pCol, 1. - smoothstep(0., sf, d + ew));
    
    // Applying the simple lighting.
    col *= .5 + diff + spec*8.;
    
    // Cheap specular reflections.
    float speR = pow(max(dot(normalize(ld - rd), n), 0.), 8.);
    vec3 rf = reflect(rd, n); // Surface reflection.
    vec2 rf2 = lonlat(rf)
    vec3 rTx = texture(iChannel1, rf2).xyz; rTx *= rTx;
    // Hacky sunset colors.
    rTx = mix(rTx*vec3(1, 1, .25), rTx.zyx, oP.y*.5 + .5);
    col = col + (col)*speR*rTx*16.;

     // Subtle vignette. Designers use them to frame things and guide
    // the viewer's eyes toward the center... or something like that.
    vec2 w = vec2(iResolution.x/iResolution.y, 1)*2.;
    col *= 1.1 - smoothstep(0., .15, sBox(oP, w/2., .15) + .15)*.4;
    

    // Rough gamma correction and screen presentation.
    fragColor = vec4(pow(max(col, 0.), vec3(1)/2.2), 1);
    
}
