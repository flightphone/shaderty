/**
 * Shukhov Tower: The Geometry of Ruled Surfaces
 * * This shader visualizes a one-sheeted hyperboloid as a "doubly ruled surface".
 * The structure is built entirely from straight lines (generators).
 * * Core logic and mathematics provided by Gemini (Google AI).
 * this shader wrated by Gemini AI
 */

// Hyperboloid parameters
const float A = 0.7;          // Throat radius (X-axis)
const float B = 0.7;          // Throat radius (Y-axis)
const float C = 1.2;          // Curvature/Height parameter
const float N = 14.0;         // Number of beams in each family
const float R_BEAM = 0.03;    // Beam thickness
const float HEIGHT = 30.0;    // Total tower height

#define PI 3.14159265359

// Rotation matrix for Y-axis
mat2 rot(float a) {
    float s = sin(a), c = cos(a);
    return mat2(c, s, -s, c);
}

// SDF for an infinite straight line
// p: current point, p0: point on the line, d: normalized direction vector
float sdLine(vec3 p, vec3 p0, vec3 d) {
    vec3 pa = p - p0;
    float h = dot(pa, d);
    return length(pa - h * d);
}

float map(vec3 p) {
    // Height clipping
    if (abs(p.y) > HEIGHT) {
        return length(p.xz) - A > 0.0 ? p.y - HEIGHT : distance(p, vec3(p.x, HEIGHT * sign(p.y), p.z));
    }

    // Use polar coordinates to find the nearest beam sector
    float angle = atan(p.z, p.x);
    float sector = 2.0 * PI / N;
    float id = floor(angle / sector + 0.5);
    
    float minDist = 1e10;

    // Check neighboring sectors to ensure smooth SDF transitions
    for (float i = -1.0; i <= 1.0; i++) {
        float a = (id + i) * sector;
        
        // P0 point on the throat ellipse (y = 0)
        vec3 p0 = vec3(A * cos(a), 0.0, B * sin(a));
        
        // Direction vectors for the two families of generators.
        // The combination of tangent and vertical vectors creates the 
        // characteristic hyperbolic curvature from straight lines.
        vec3 tangent = vec3(-A * sin(a), 0.0, B * cos(a));
        vec3 vertical = vec3(0.0, C, 0.0); 
        
        // Family 1: clockwise twist, Family 2: counter-clockwise twist
        vec3 d1 = normalize(vertical + tangent);
        vec3 d2 = normalize(vertical - tangent);
        
        minDist = min(minDist, sdLine(p, p0, d1));
        minDist = min(minDist, sdLine(p, p0, d2));
    }

    // Add horizontal structural rings every 1.0 unit of height
    // Radius at height Y is calculated using: r = sqrt(A^2 + (y^2/C^2)*A^2)
    float towerRadiusAtY = sqrt(A*A + (p.y*p.y)/(C*C)*A*A);
    float rings = length(vec2(length(p.xz) - towerRadiusAtY, mod(min(abs(p.y), 1.5) + 0.5, 1.0) - 0.5));
    minDist = min(minDist, rings);

    return minDist - R_BEAM;
}

// Standard normal calculation for lighting
vec3 getNormal(vec3 p) {
    vec2 e = vec2(0.001, 0.0);
    return normalize(vec3(
        map(p + e.xyy) - map(p - e.xyy),
        map(p + e.yxy) - map(p - e.yxy),
        map(p + e.yyx) - map(p - e.yyx)
    ));
}

// Ray direction setup for a look-at camera
vec3 GetRayDir(vec2 uv, vec3 p, vec3 l, float z) {
    vec3 f = normalize(l-p),
         r = normalize(cross(vec3(0,1,0), f)),
         u = cross(f,r),
         c = f*z,
         i = c + uv.x*r + uv.y*u;
    return normalize(i);
}

void mainImage(out vec4 fragColor, in vec2 fragCoord) {
    vec2 uv = (fragCoord - 0.5 * iResolution.xy) / iResolution.y;
    vec3 ro = vec3(0.0, 0.0, -5.0); // Camera origin
    const float fl = 1.5;          // Focal length
    
    // Camera rotation via mouse or time
    float ang = iTime * 0.3;
    if (iMouse.z > 0.0) {
        vec2 m = (-iResolution.xy + 2.0*(iMouse.xy))/iResolution.y;
        ro.xz *= -rot(m.x * PI * 2.0);
        ro.xy *= rot(m.y * PI * 2.0);
    } else {
        ro.xz *= rot(ang);
    }

    vec3 rd = GetRayDir(uv, ro, vec3(0, 0.0, 0), fl);

    // Raymarching loop
    float t = 0.0;
    for (int i = 0; i < 100; i++) {
        vec3 p = ro + rd * t;
        float d = map(p);
        if (d < 0.001 || t > 40.0) break;
        t += d;
    }

    // Shading and Background
    vec3 col = vec3(0.05, 0.05, 0.1); // Night sky background
    if (t < 40.0) {
        vec3 p = ro + rd * t;
        vec3 n = getNormal(p);
        float diff = max(dot(n, normalize(vec3(1, 2, -1))), 0.0);
        float amb = 0.2;
        col = vec3(0.7, 0.7, 0.8) * (diff + amb); // Metallic look
        
        // Fog/depth fade
        col = mix(col, vec3(0.05, 0.05, 0.1), 1.0 - exp(-0.05 * t));
    }

    // Gamma correction
    col = pow(col, vec3(0.4545));
    fragColor = vec4(col, 1.0);
}