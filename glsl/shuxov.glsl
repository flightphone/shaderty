

// Параметры гиперболоида (башни)
const float A = .7;         // Радиус горловины (по X)
const float B = .7;         // Радиус горловины (по Y)
const float C = 1.2;         // Параметр высоты
const float N = 14.0;        // Количество балок в одном семействе
const float R_BEAM = 0.03;   // Толщина балки
const float HEIGHT = 25.0;    // Общая высота башни

#define PI 3.14159265359

// Функция поворота вокруг оси Y
mat2 rot(float a) {
    float s = sin(a), c = cos(a);
    return mat2(c, s, -s, c);
}

// SDF одной бесконечной прямой
// p - точка, p0 - точка на прямой, d - нормализованный вектор направления
float sdLine(vec3 p, vec3 p0, vec3 d) {
    vec3 pa = p - p0;
    float h = dot(pa, d);
    return length(pa - h * d);
}

float map(vec3 p) {
    // Ограничиваем башню по высоте
    if (abs(p.y) > HEIGHT) return length(p.xz) - A > 0.0 ? p.y - HEIGHT : distance(p, vec3(p.x, HEIGHT, p.z));

    float angle = atan(p.z, p.x);
    float sector = 2.0 * PI / N;
    float id = floor(angle / sector + 0.5);
    
    float minDist = 1e10;

    for (float i = -1.0; i <= 1.0; i++) {
        float a = (id + i) * sector;
        
        // Точка P0 на горловине (y = 0)
        vec3 p0 = vec3(A * cos(a), 0.0, B * sin(a));
        
        // Направляющие векторы образующих. 
        // Чтобы получить гиперболоид, вектор должен быть суммой 
        // тангенциального вектора к окружности и вертикального вектора C.
        // Это и есть ключ к линейчатой поверхности:
        vec3 tangent = vec3(-A * sin(a), 0.0, B * cos(a));
        vec3 vertical = vec3(0.0, C, 0.0); 
        
        // Два семейства: закрутка вправо и влево
        vec3 d1 = normalize(vertical + tangent);
        vec3 d2 = normalize(vertical - tangent);
        
        minDist = min(minDist, sdLine(p, p0, d1));
        minDist = min(minDist, sdLine(p, p0, d2));
    }

    // Добавляем кольца через каждые 1.0 единицы высоты
    float rings = length(vec2(length(p.xz) - sqrt(A*A + (p.y*p.y)/(C*C)*A*A), mod(min(abs(p.y), 1.5) + 0.5, 1.) - 0.5));
    minDist = min(minDist, rings);

    return minDist - R_BEAM;
}
// Простая нормаль для освещения
vec3 getNormal(vec3 p) {
    vec2 e = vec2(0.001, 0.0);
    return normalize(vec3(
        map(p + e.xyy) - map(p - e.xyy),
        map(p + e.yxy) - map(p - e.yxy),
        map(p + e.yyx) - map(p - e.yyx)
    ));
}

vec3 GetRayDir(vec2 uv, vec3 p, vec3 l, float z) {
    vec3 
        f = normalize(l-p),
        r = normalize(cross(vec3(0,1,0), f)),
        u = cross(f,r),
        c = f*z,
        i = c + uv.x*r + uv.y*u;
    return normalize(i);
}

void mainImage(out vec4 fragColor, in vec2 fragCoord) {
    vec2 uv = (fragCoord - 0.5 * iResolution.xy) / iResolution.y;
    vec3 ro = vec3(0.0, 0.0, -5.0); // Позиция камеры
    const float fl = 1.5; // focal length
    
    //vec3 rd = normalize(vec3(uv, 1.2)); // Луч
    
    // Вращение камеры мышкой или временем
    float ang = iTime * 0.3;
    vec2 m = vec2(0.0, 0.0);
    if (iMouse.z > 0.0) 
    {
        m = (-iResolution.xy + 2.0*(iMouse.xy))/iResolution.y;
        ro.xz *= -rot(m.x*PI*2.);
        ro.xy *= rot(m.y*PI*2.);
    }
    else
    {
        ro.xz *= rot(ang);
    }

    vec3 rd = GetRayDir(uv, ro, vec3(0,0.,0), fl);

    
    
    

    // Raymarching
    float t = 0.0;
    for (int i = 0; i < 100; i++) {
        vec3 p = ro + rd * t;
        float d = map(p);
        if (d < 0.001 || t > 20.0) break;
        t += d;
    }

    // Отрисовка
    vec3 col = vec3(0.05, 0.05, 0.1); // Цвет фона (ночь)
    if (t < 20.0) {
        vec3 p = ro + rd * t;
        vec3 n = getNormal(p);
        float diff = max(dot(n, normalize(vec3(1, 2, -1))), 0.0);
        float amb = 0.2;
        col = vec3(0.7, 0.7, 0.8) * (diff + amb); // Цвет металла
        
        // Туман для глубины
        col = mix(col, vec3(0.05, 0.05, 0.1), 1.0 - exp(-0.05 * t));
    }

    // Гамма-коррекция
    col = pow(col, vec3(0.4545));
    fragColor = vec4(col, 1.0);
}