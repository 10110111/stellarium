uniform mediump mat4 projectionMatrix;

attribute vec2 vertex;
attribute vec4 viewRayAndAddLuminance;

varying vec3 view_ray;
varying vec3 additionalLuminance;

void main()
{
    view_ray = viewRayAndAddLuminance.xyz;
	gl_Position = projectionMatrix*vec4(vertex, 0., 1.);

    // Additional luminance is airglow and light pollution. It's mostly airglow and
    // light pollution, so should be in the red and green parts of the spectrum.
	additionalLuminance = viewRayAndAddLuminance.w*vec3(1,1,0);
}

