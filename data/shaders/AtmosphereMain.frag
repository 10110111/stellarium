#version 130

const float kLengthUnitInMeters = 1000.000000;
const float earthRadius=6.36e6/kLengthUnitInMeters;
const float sunAngularRadius=0.00935/2;

// This one is non-const because some broken drivers (e.g. nvidia-340, which is
// the latest for NVIDIA ION) think tan&cos are non-constant-expressions
vec2 sun_size=vec2(tan(sunAngularRadius),cos(sunAngularRadius));
const vec3 earth_center=vec3(0,0,-earthRadius);

uniform vec3 camera;
uniform vec3 sun_direction;

in vec3 additionalLuminance;
in vec3 view_ray;
out vec4 color;
vec3 GetSolarLuminance();
vec3 GetSkyLuminance(vec3 camera, vec3 view_ray, float shadow_length,
                     vec3 sun_direction, out vec3 transmittance);
float sqr(float x) { return x*x; }

// cos(zenithAngle) for a ray tangential to horizon and originating at given altitude
float horizonMu(float altitude)
{
    float R=earthRadius;
    float h=altitude;
    return -sqrt(2*h*R+sqr(h))/(R+h);
}

void main()
{
    vec3 view_direction=normalize(view_ray);
    vec3 cameraPos=camera;

    // FIXME: we simply clamp negative altitudes to zero currently. This is not
    // quite correct: there are places with negative elevation above sea level.
    // But at least this will avoid breakage of the model.
    if(cameraPos.z<0)
        cameraPos.z=0;

    // Don't draw any type of ground. Instead reflect the ray from mathematical
    // horizon. This isn't physical, but lets Tone Reproducer work without
    // overexposures.
    vec3 p = cameraPos - earth_center;
    float p_dot_v = dot(p, view_direction);
    float p_dot_p = dot(p, p);
    float ray_earth_center_squared_distance = p_dot_p - sqr(p_dot_v);
    float distance_to_intersection = -p_dot_v - sqrt(
      sqr(earth_center.z) - ray_earth_center_squared_distance);
    // cameraPos.z==0 is a special case where distance to intersection calculation
    // is unreliable (has a lot of noise in its sign), so check it separately
    if(distance_to_intersection>0 || (cameraPos.z==0 && view_direction.z<0))
    {
        view_direction.z=2*horizonMu(cameraPos.z)-view_direction.z;
    }

    vec3 transmittance;
    vec3 radiance = GetSkyLuminance(cameraPos - earth_center, view_direction, 0, sun_direction, transmittance);
    /*
    // This 'if' case is useful when we want to check that the Sun is where
    // it's supposed to be and, ideally, has the color it's supposed to have
    if (dot(view_direction, sun_direction) > sun_size.y)
        radiance += transmittance * GetSolarLuminance();
     */

    color = vec4(radiance+additionalLuminance,1);
}
