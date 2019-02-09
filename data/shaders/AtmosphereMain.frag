#version 130

const float kLengthUnitInMeters = 1000.000000;
const float earthRadius=6.36e6/kLengthUnitInMeters;
const float sunAngularRadius=0.00935/2;

const vec2 sun_size=vec2(tan(sunAngularRadius),cos(sunAngularRadius));
const vec3 earth_center=vec3(0,0,-earthRadius);
const float exposure=0.6;

uniform vec3 camera;
uniform vec3 sun_direction;

in vec3 view_ray;
out vec4 color;
vec3 GetSolarLuminance();
vec3 GetSkyLuminance(vec3 camera, vec3 view_ray, float shadow_length,
                     vec3 sun_direction, out vec3 transmittance);
vec3 adaptColorToVisionMode(vec3 rgb);
vec3 dither(vec3);
void main()
{
    vec3 view_direction=normalize(view_ray);
    float fragment_angular_size = length(dFdx(view_direction) + dFdy(view_direction));
    vec3 transmittance;
    vec3 radiance = GetSkyLuminance(camera - earth_center, view_direction, 0, sun_direction, transmittance);
    /*
    // This 'if' case is useful when we want to check that the Sun is where
    // it's supposed to be and, ideally, has the color it's supposed to have
    if (dot(view_direction, sun_direction) > sun_size.y)
        radiance += transmittance * GetSolarLuminance();
     */
    color = vec4(dither(pow(adaptColorToVisionMode(radiance * exposure), vec3(1.0 / 2.2))),1);
}
