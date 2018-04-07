R"(uniform vec3 rgbMaxValue;
uniform sampler2D bayerPattern;
vec3 dither(vec3 c)
{
    float bayer=texture2D(bayerPattern,gl_FragCoord.xy/8.).r;

    vec3 rgb=c*rgbMaxValue;
    vec3 head=floor(rgb);
    vec3 tail=fract(rgb);
    return (head+step(bayer,tail))/rgbMaxValue;
}
vec4 dither(vec4 c) { return vec4(dither(c.xyz),c.w); }
)"
