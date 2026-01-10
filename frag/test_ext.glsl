precision mediump float;



uniform vec2 u_resolution;
uniform sampler2D u_texture_0; 
void main() {
    vec2 uv = gl_FragCoord.xy / u_resolution.xy;
    vec4 color = texture2D(u_texture_0, uv);
    gl_FragColor = color;
}