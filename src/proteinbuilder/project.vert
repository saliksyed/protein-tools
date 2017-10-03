// energy use for this vertex (f
varying float distToCamera;

void main()
{
    vec4 vertex = gl_Vertex;
    vec4 cs_position = gl_ModelViewMatrix * vertex;
    distToCamera = -cs_position.z;
    gl_Position = gl_ModelViewProjectionMatrix * vertex;
}