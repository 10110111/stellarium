uniform mediump mat4 projectionMatrix;

attribute vec2 vertex;
attribute vec3 viewRay;

varying vec3 view_ray;

void main()
{
    view_ray = viewRay.xyz;
	gl_Position = projectionMatrix*vec4(vertex, 0., 1.);
}

