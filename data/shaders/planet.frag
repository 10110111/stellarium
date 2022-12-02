/*
 * Stellarium
 * Copyright (C) 2002-2016 Fabien Chereau and Stellarium contributors
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public License
 * as published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Suite 500, Boston, MA  02110-1335, USA.
 */

/*
  This is the fragment shader for solar system object rendering
 */

VARYING mediump vec2 texc; //texture coord
VARYING highp vec3 P; //original vertex pos in model space

const highp float PI = 3.14159265;

uniform sampler2D tex;
uniform mediump vec2 poleLat; //latitudes of pole caps, in terms of texture coordinate. x>0...north, y<1...south. 
uniform mediump vec3 ambientLight;
uniform mediump vec3 diffuseLight;
uniform highp vec4 sunInfo;
uniform mediump float skyBrightness;

uniform int shadowCount;
uniform highp mat4 shadowData;

//x = scaling, y = exponential falloff
uniform mediump vec2 outgasParameters;
//eye direction in model space, pre-normalized
uniform highp vec3 eyeDirection;

#ifdef RINGS_SUPPORT
uniform bool ring;
uniform highp float outerRadius;
uniform highp float innerRadius;
uniform sampler2D ringS;
uniform bool isRing;
#endif

#ifdef SHADOWMAP
uniform highp sampler2D shadowTex;
VARYING highp vec4 shadowCoord;
#endif

#if defined(IS_OBJ)
    #define OREN_NAYAR 1
    //light direction in model space, pre-normalized
    uniform highp vec3 lightDirection;  
    //x = A, y = B, z = scaling factor (rho/pi * E0), w roughness
    uniform mediump vec4 orenNayarParameters;
#endif

#if defined(IS_MOON)
# define HAPKE 1
//light direction in model space, pre-normalized
uniform highp vec3 lightDirection;
#endif

#ifdef IS_MOON
    uniform sampler2D earthShadow;
    uniform mediump float eclipsePush;
    uniform sampler2D normalMap;
	uniform sampler2D horizonMap;

    VARYING highp vec3 normalX;
    VARYING highp vec3 normalY;
    VARYING highp vec3 normalZ;
#else
    VARYING mediump float lambertIllum;
    VARYING mediump vec3 normalVS; //pre-calculated normals or spherical normals in model space
#endif

const highp float M_PI=3.1415926535897932384626433832795;

#ifdef SHADOWMAP
uniform highp vec2 poissonDisk[64];

lowp float offset_lookup(in highp sampler2D sTex, in highp vec4 loc, in highp vec2 offset, in highp float zbias)
{
    //the macro SM_SIZE is set to the shadowmap size
    const mediump vec2 texmapscale=vec2(1.0/float(SM_SIZE));
    //"simulates" textureProjOffset for use in GLSL < 130
    highp vec4 coords = vec4(loc.xy + (offset * texmapscale * loc.w), loc.z, loc.w);
    
    //for some reason, when not adding a LOD bias, the result is wrong in my VM (Ubuntu 16.04)
    //even if the texture has NO mipmaps and the lookup filter is GL_NEAREST???
    //I'm 99% certain this is some bug in the GL driver, 
    //because adding ANY bias here makes it work correctly
    //It should have no effect on platforms which don't have this bug
    highp float texVal = texture2DProj_3(sTex, coords, -1000.0).r;
    //perform shadow comparison
    return texVal > (loc.z-zbias)/loc.w ? 1.0 : 0.0;
}

//basic pseudo-random number generator
mediump float random(in mediump vec4 seed4)
{
    highp float dot_product = dot(seed4, vec4(12.9898,78.233,45.164,94.673));
    return fract(sin(dot_product) * 43758.5453);
}

lowp float sampleShadowMap(in highp sampler2D sTex, in highp vec4 coord, in highp float zbias)
{   
    //uncomment for a single sample
    //return offset_lookup(sTex,coord, vec2(0.0),zbias);
 
    // for some reason > 5 samples do not seem to work on my Ubuntu VM 
    // (no matter if ES2 or GL 2.1)
    // everything gets shadowed, but no errors?!
    // so to be sure we just fix the sample count at 4 for now, 
    // even though 16 would look quite a lot better
    const int SAMPLE_COUNT = 4;
    mediump float sum = 0.0;
    for(int i=0;i<SAMPLE_COUNT;++i)
    {
        //choose "random" locations out of 64 using the screen-space coordinates
        //for ES2, we have to use mod(float, float), mod(int, int) is not defined
        int index = int(mod( 64.0*random(vec4(gl_FragCoord.xyy,i)), 64.0));
        sum += offset_lookup(sTex, coord, poissonDisk[index], zbias);
    }
    
    return clamp(sum / float(SAMPLE_COUNT),0.0,1.0);
}
#endif

#ifdef OREN_NAYAR
// Calculates the Oren-Nayar reflectance (https://en.wikipedia.org/wiki/Oren%E2%80%93Nayar_reflectance_model)
// the scale parameter is actually rho/pi * E_0 here
// A and B are precalculated on the CPU side
mediump float orenNayar(in mediump vec3 normal, in highp vec3 lightDir, in highp vec3 viewDir, in mediump float A, in mediump float B, in mediump float scale, in mediump float roughSq)
{
    mediump float cosAngleLightNormal = dot(normal, lightDir);  //cos theta_i
    mediump float cosAngleEyeNormal = dot(normal, viewDir); //cos theta_r
    //acos can be quite expensive, can we avoid it?
    mediump float angleLightNormal = acos(cosAngleLightNormal); //theta_i
    mediump float angleEyeNormal = acos(cosAngleEyeNormal); //theta_r
    mediump float alpha = max(angleEyeNormal, angleLightNormal); //alpha = max(theta_i, theta_r)
    mediump float beta = min(angleEyeNormal, angleLightNormal); //beta = min(theta_i, theta_r)
    mediump float gamma = dot(viewDir - normal * cosAngleEyeNormal, lightDir - normal * cosAngleLightNormal); // cos(phi_r-phi_i)
    mediump float C = sin(alpha) * tan(beta);
	mediump float ON = max(0.0, cosAngleLightNormal) * ((A + B * max(0.0, gamma) * C) * scale); // Qualitative model done.
    // Now add third term:
    mediump float C3_2 = (4.0*alpha*beta)/(M_PI*M_PI);
    mediump float C3=0.125*roughSq/(roughSq+0.09)*C3_2*C3_2;
    mediump float third=(1.0-abs(gamma))*C3*tan(0.5*(alpha+beta))*max(0.0, cosAngleLightNormal)*scale;
    ON += third;
    // Add the intereflection term:
    mediump float betaterm = 2.0*beta/M_PI;
    mediump float ONir = max(0.0, cosAngleLightNormal) * ((1.0-gamma*betaterm*betaterm)*0.17*scale*roughSq/(roughSq+0.13));
    ON += ONir;
    // This clamp gives an ugly pseudo edge visible in rapid animation.
    //return clamp(ON, 0.0, 1.0);
    // Better use GLSL's smoothstep. However, this is N/A in GLSL1.20/GLES. Workaround from https://www.khronos.org/registry/OpenGL-Refpages/gl4/html/smoothstep.xhtml
    if (ON>0.8)
    {
        mediump float t = clamp((ON - 0.8) / (1.0 - 0.8), 0.0, 1.0);
        return 0.8 + 0.2*(t * t * (3.0 - 2.0 * t));
    }
    else
    return clamp(ON, 0.0, 1.0);
}
#endif
#ifdef HAPKE
struct HapkeParameters
{
    float theta_p;
    float phi;
    float h_C;
    float B_C0;
    float h_S;
    float B_S0;
    float w;
    float b;
    float c;
};

float hapkeBRDF(vec3 N, vec3 L, vec3 V, HapkeParameters params)
{
    float PI = 3.14159265358979323846;
    float cos_g = dot(L, V);
    float cos_i = dot(L, N); // N=L => 1
    float cos_e = dot(V, N); // N=L => cos_g
    float g = acos(clamp(cos_g, -1., 1.));
    float i = acos(clamp(cos_i, -1., 1.));
    float e = acos(clamp(cos_e, -1., 1.));
    float my_0 = cos_i;
    float my = cos_e;
    float sin_i = sin(i);
    float sin_e = sin(e);
    // Calculate Psi by projecting L and V on a plane with N as the normal and getting
    // the cosine of the angle between the projections.
    //float Psi = acos((cos_g - cos_e * cos_i) / (sin_e * sin_i));
    float cos_Psi = dot(normalize(L - cos_i * N), normalize(V - cos_e * N));
    cos_Psi = clamp(cos_Psi, -1., 1.);
	if(isnan(cos_Psi) || isinf(cos_Psi))
		cos_Psi = 1.;
    float Psi = acos(cos_Psi);
    float PsiHalf = Psi / 2;
    float f_Psi = exp(-2 * sqrt((1-cos_Psi)/(1+cos_Psi))); // = exp(-2 * tan(PsiHalf)) but preventing tan()<0 that can happen due to bad acos implementation (e.g. Intel UHD 620 @ Mesa)
    float sin_PsiHalf = sin(PsiHalf);
    float sin2_PsiHalf = sin_PsiHalf * sin_PsiHalf;
    float PsiPerPI = Psi / PI;
    float tan_theta_p = tan(params.theta_p * PI / 180);
    float tan2_theta_p = tan_theta_p * tan_theta_p;
    // The tan_theta_p is zero when theta_p is zero.
    // A zero for theta_p is somewhat unrealistic.
    float cot_theta_p = 1 / tan_theta_p;
    float cot2_theta_p = cot_theta_p * cot_theta_p;
    float tan_i = tan(i);
    float tan_e = tan(e);
    // Here, tan_i and tan_e are zero each when i and e are zero each.
    // Because 1 / 0.0 is positive infinity we should check the usages.
    // The results are only used as factors for exponential functions,
    // where the whole expressions for the arguments are negated and
    // therefore result in negative infinity. Thus overall, these
    // exponential functions simply result in zero. No fix needed.
    float cot_i = 1 / tan_i;
    float cot_e = 1 / tan_e;
    float cot2_i = cot_i * cot_i;
    float cot2_e = cot_e * cot_e;

    float E_1_i = exp(-2 / PI * cot_theta_p * cot_i);
    float E_1_e = exp(-2 / PI * cot_theta_p * cot_e);
    float E_2_i = exp(-1 / PI * cot2_theta_p * cot2_i);
    float E_2_e = exp(-1 / PI * cot2_theta_p * cot2_e);
    float chi_theta_p = 1 / sqrt(1 + PI * tan2_theta_p);
    float eta_i = chi_theta_p * (cos_i + sin_i * tan_theta_p * E_2_i / (2 - E_1_i));
    float eta_e = chi_theta_p * (cos_e + sin_e * tan_theta_p * E_2_e / (2 - E_1_e));
    float my_0e = chi_theta_p;
    float my_e = chi_theta_p;
    float S;
    // It looks worse than it is.
    if(i <= e)
    {
        my_0e *= cos_i + sin_i * tan_theta_p * (cos_Psi * E_2_e + sin2_PsiHalf * E_2_i) / (2 - E_1_e - PsiPerPI * E_1_i);
        my_e  *= cos_e + sin_e * tan_theta_p * (E_2_e - sin2_PsiHalf * E_2_i) / (2 - E_1_e - PsiPerPI * E_1_i);
        S = my_e / eta_e * my_0 / eta_i * chi_theta_p / (1 - f_Psi + f_Psi * chi_theta_p * (my_0 / eta_i));
    }
    else
    {
        my_0e *= cos_i + sin_i * tan_theta_p * (E_2_i - sin2_PsiHalf * E_2_e) / (2 - E_1_i - PsiPerPI * E_1_e);
        my_e  *= cos_e + sin_e * tan_theta_p * (cos_Psi * E_2_i + sin2_PsiHalf * E_2_e) / (2 - E_1_i - PsiPerPI * E_1_e);
        S = my_e / eta_e * my_0 / eta_i * chi_theta_p / (1 - f_Psi + f_Psi * chi_theta_p * (my / eta_e));
    }
    float KphiTerm = 1.209 * pow(params.phi, 2.0 / 3);
    // This goes into complex numbers already within [0, 1].
    float K = -log(1 - KphiTerm) / KphiTerm;
    float gHalf = g / 2;
    float tan_gHalf = tan(gHalf);
    float tan_gHalfPerh_C = tan_gHalf / params.h_C;
    float B_C = (1 + (1 - exp(-tan_gHalfPerh_C)) / tan_gHalfPerh_C) / (2 * pow(1 + tan_gHalfPerh_C, 2));
    // When g = 0, the division in the upper part causes a NaN.
    if(isnan(B_C))
    {
        B_C = 1;
    }
    float r_0Term = sqrt(1 - params.w);
    float r_0 = (1 - r_0Term) / (1 + r_0Term);
    float LS = my_0e / (my_0e + my_e);
    float b2 = params.b * params.b;
    // An approximation from the publication.
    //c = 3.29 * exp(-17.4 * b2) - 0.908;
    float oneMinusb2 = 1 - b2;
    float twobcos_g = 2 * params.b * cos_g;

    float p_g = (1 + params.c) / 2 * oneMinusb2 / pow(1 - twobcos_g + b2, 1.5) +
                (1 - params.c) / 2 * oneMinusb2 / pow(1 + twobcos_g + b2, 1.5);
    float B_S = 1 / (1 + tan_gHalf / params.h_S);
    float x_i = my_0e / K;
    float x_e = my_e / K;
    float H_i = 1 / (1 - params.w * x_i * (r_0 + (1 - 2 * r_0 * x_i) / 2 * log((1 + x_i) / x_i)));
    float H_e = 1 / (1 - params.w * x_e * (r_0 + (1 - 2 * r_0 * x_e) / 2 * log((1 + x_e) / x_e)));
    float M = H_i * H_e - 1;
    return LS * K * params.w / (4 * PI) * (p_g * (1 + params.B_S0 * B_S) + M) * (1 + params.B_C0 * B_C) * S / cos_i;
}
#endif

// calculate pseudo-outgassing effect, inspired by MeshLab's "electronic microscope" shader
lowp float outgasFactor(in mediump vec3 normal, in highp vec3 lightDir, in mediump float falloff)
{
    mediump float opac = dot(normal,lightDir);
    opac = abs(opac);
    opac = 1.0 - pow(opac, falloff);
    return opac;
}

void main()
{
    mediump float final_illumination = 1.0;
#if defined OREN_NAYAR || defined HAPKE
    mediump float lum = 1.;
#else
    mediump float lum = lambertIllum;
#endif
#ifdef RINGS_SUPPORT
    if(isRing)
        lum=1.0;
#endif
    //shadow calculation
    if(lum > 0.0)
    {
        highp vec3 sunPosition = sunInfo.xyz;
#ifdef RINGS_SUPPORT
        if(ring && !isRing)
        {
            highp vec3 ray = normalize(sunPosition);
            highp float u = - P.z / ray.z;
            if(u > 0.0 && u < 1e10)
            {
                mediump float ring_radius = length(P + u * ray);
                if(ring_radius > innerRadius && ring_radius < outerRadius)
                {
                    ring_radius = (ring_radius - innerRadius) / (outerRadius - innerRadius);
                    lowp float ringAlpha = texture2D(ringS, vec2(ring_radius, 0.5)).w;
                    final_illumination = 1.0 - ringAlpha;
                }
            }
        }
#endif

        highp float sunRadius = sunInfo.w;
        highp float L = length(sunPosition - P);
        highp float R = asin(sunRadius / L);

        for (int i = 0; i < 4; ++i)
        {
            if (shadowCount>i)
            {
                highp vec3 satellitePosition = shadowData[i].xyz;
                highp float satelliteRadius = shadowData[i].w;
                highp float l = length(satellitePosition - P);
                highp float r = asin(satelliteRadius / l);
                // The min(1.0, .) here is required to avoid spurious bright pixels close to the shadow center.
                highp float d = acos(min(1.0, dot(normalize(sunPosition - P), normalize(satellitePosition - P))));

                mediump float illumination = 1.0;
                if(d >= R + r) // distance too far
                {
                    // illumination = 1.0; // NOP
                }
                else if( d <= r - R ) // fully inside umbra
                {
#ifdef IS_MOON
                    illumination = (d / (r - R)) * 0.594; // prepare texture coordinate. 0.6=umbra edge. Smaller number->larger shadow.
#else
                    illumination = 0.0;
#endif
                }
                else if(d <= R - r)
                {
                    // penumbra completely inside (Annular, or a moon transit in front of the Sun.)
                    illumination = 1.0 - r * r / (R * R);
                }
                else // penumbra: partially inside
                {
#ifdef IS_MOON
                    //illumination = ((d - abs(R-r)) / (R + r - abs(R-r))) * 0.4 + 0.6;
                    illumination = ((d - r + R) / (2.0 * R )) * 0.406 + 0.594;
#else
                    mediump float x = (R * R + d * d - r * r) / (2.0 * d);
                    mediump float alpha = acos(x / R);
                    mediump float beta = acos((d - x) / r);
                    mediump float AR = R * R * (alpha - 0.5 * sin(2.0 * alpha));
                    mediump float Ar = r * r * (beta - 0.5 * sin(2.0 * beta));
                    mediump float AS = R * R * 2.0 * 1.57079633;
                    illumination = 1.0 - (AR + Ar) / AS;
#endif
                }
                final_illumination = min(illumination, final_illumination);
            }
        }
    }

#ifdef IS_MOON
    mediump vec2 moonTexCoord = vec2(atan(normalZ.x, -normalZ.y)/(2.*PI)+0.5, asin(normalize(normalZ).z)/PI+0.5);
    mediump vec3 normal = texture2D(normalMap, moonTexCoord).rgb-vec3(0.5, 0.5, 0);
    normal = normalize(normalX*normal.x+normalY*normal.y+normalZ*normal.z);
    // normal now contains the real surface normal taking normal map into account

	mediump float horizonShadowCoefficient = 1.;
	{
		// Check whether the fragment is in the shadow of surrounding mountains or the horizon
		mediump vec3 lonDir = normalX;
		mediump vec3 northDir = normalY;
		mediump vec3 zenith = normalZ;
		mediump float sunAzimuth = atan(dot(lightDirection,lonDir), dot(lightDirection,northDir));
		mediump float sinSunElevation = dot(zenith, lightDirection);
		mediump vec4 horizonElevSample = (texture2D(horizonMap, moonTexCoord) - 0.5) * 2.;
		mediump vec4 sinHorizElevs = sign(horizonElevSample) * horizonElevSample * horizonElevSample;
		mediump float sinHorizElevLeft, sinHorizElevRight;
		mediump float alpha;
		if(sunAzimuth >= PI/2.)
		{
			// Sun is between East and South
			sinHorizElevLeft = sinHorizElevs[1];
			sinHorizElevRight = sinHorizElevs[2];
			alpha = (sunAzimuth - PI/2.) / (PI/2.);
		}
		else if(sunAzimuth >= 0.)
		{
			// Sun is between North and East
			sinHorizElevLeft = sinHorizElevs[0];
			sinHorizElevRight = sinHorizElevs[1];
			alpha = sunAzimuth / (PI/2.);
		}
		else if(sunAzimuth <= -PI/2.)
		{
			// Sun is between South and West
			sinHorizElevLeft = sinHorizElevs[2];
			sinHorizElevRight = sinHorizElevs[3];
			alpha = (sunAzimuth + PI) / (PI/2.);
		}
		else
		{
			// Sun is between West and North
			sinHorizElevLeft = sinHorizElevs[3];
			sinHorizElevRight = sinHorizElevs[0];
			alpha = (sunAzimuth + PI/2.) / (PI/2.);
		}
		mediump float horizElevLeft = asin(sinHorizElevLeft);
		mediump float horizElevRight = asin(sinHorizElevRight);
		mediump float horizElev = horizElevLeft + (horizElevRight-horizElevLeft)*alpha;
		if(sinSunElevation < sin(horizElev))
			horizonShadowCoefficient = 0.;
	}
#else
    // important to normalize here again
    mediump vec3 normal = normalize(normalVS);
#endif
#ifdef OREN_NAYAR
    // Use an Oren-Nayar model for rough surfaces
    // Ref: http://content.gpwiki.org/index.php/D3DBook:(Lighting)_Oren-Nayar
    lum = orenNayar(normal, lightDirection, eyeDirection, orenNayarParameters.x, orenNayarParameters.y, orenNayarParameters.z, orenNayarParameters.w);
#endif
#ifdef HAPKE
	{
		HapkeParameters hapkeParams;
		hapkeParams.w = 0.32357;
		hapkeParams.b = 0.23955;
		hapkeParams.c = 0.30452;
		hapkeParams.B_C0 = 0.0;
		hapkeParams.h_C = 1.0;
		hapkeParams.B_S0 = 1.80238;
		hapkeParams.h_S = 0.07145;
		hapkeParams.theta_p = 23.4;
		hapkeParams.phi = 0.3;
		float f_r = hapkeBRDF(normal, lightDirection, eyeDirection, hapkeParams);
		if(isnan(f_r) || f_r < 0 || isinf(f_r))
		{
			lum = 0;
		}
		else
		{
			vec3 L = lightDirection, V = eyeDirection;
			float dotLV = dot(L, V);
			vec3 refNormal = dotLV > 0. ? L /* N=L */
										: normalize(L - dotLV * V) /* clamp to the rim, so that N⊥V */;
			float ref = hapkeBRDF(refNormal, lightDirection, eyeDirection, hapkeParams);
			float cos_ref_i = dot(refNormal, lightDirection);
			float illuminance = 1. / (ref * max(0., cos_ref_i));

			float cos_i = dot(normal, lightDirection);
			lum = max(0., cos_i) * f_r * illuminance;
		}
	}
#endif
#ifdef IS_MOON
	lum *= horizonShadowCoefficient;
#endif
    //calculate pseudo-outgassing/rim-lighting effect
    lowp float outgas = 0.0;
    if(outgasParameters.x > 0.0)
    {
        outgas = outgasParameters.x * outgasFactor(normal, eyeDirection, outgasParameters.y);
    }
//Reduce lum if sky is bright, to avoid burnt-out look in daylight sky.
    //lum *= (1.0-0.4*skyBrightness);
    lum *= clamp((0.9-0.05*log(skyBrightness)), 0.1, 0.9);
#ifdef SHADOWMAP
    //use shadowmapping
    //z-bias is modified using the angle between the surface and the light
    //gives less shadow acne
    highp float NdotL = clamp(dot(normal, lightDirection), 0.0, 1.0);
    highp float zbias = 0.01 * tan(acos(NdotL));
    zbias = clamp(zbias, 0.0, 0.03);
    lowp float shadow = sampleShadowMap(shadowTex, shadowCoord, zbias);
    lum*=shadow;
#endif

    //final lighting color
    mediump vec4 litColor = vec4(lum * final_illumination * diffuseLight + ambientLight, 1.0);

    //apply texture-colored rimlight
    //litColor.xyz = clamp( litColor.xyz + vec3(outgas), 0.0, 1.0);

#ifdef IS_MOON
    lowp vec4 texColor = texture2D(tex, moonTexCoord);
    // Undo the extraneous gamma encoded in the texture.
    // FIXME: ideally, we want all the calculations to be done in linear scale,
    // and then *in the end* apply the exact sRGB transfer function. Currently
    // though, we don't do actual physically-correct simulation here, so we
    // just apply the approximate sRGB gamma.
    texColor = pow(texColor, vec4(2.8/2.2));
#else
    lowp vec4 texColor = texture2D(tex, texc);
#endif

    mediump vec4 finalColor = texColor;
	// apply (currently only Martian) pole caps. texc.t=0 at south pole, 1 at north pole. 
	if (texc.t>poleLat.x-0.01+0.001*sin(texc.s*18.*M_PI)) {	// North pole near t=1
		mediump float mixfactor=1.;
		if (texc.t<poleLat.x+0.01+0.001*sin(texc.s*18.*M_PI))
			mixfactor=(texc.t-poleLat.x+0.01-0.001*sin(texc.s*18.*M_PI))/0.02;
		//finalColor.xyz=mix(vec3(1., 1., 1.), finalColor.xyz, 1.-mixfactor); 
		finalColor.xyz=mix(vec3(1., 1., 1.), finalColor.xyz, smoothstep(0., 1., 1.-mixfactor)); 
	}
	if (texc.t<poleLat.y+0.01+0.001*sin(texc.s*18.*M_PI)) {	// South pole near texc.t~0
		mediump float mixfactor=1.;
		if (texc.t>poleLat.y-0.01+0.001*sin(texc.s*18.*M_PI))
			mixfactor=(poleLat.y+0.01-texc.t-0.001*sin(texc.s*18.*M_PI))/0.02;
		//finalColor.xyz=mix(vec3(1., 1., 1.), finalColor.xyz, 1.-mixfactor); 
		finalColor.xyz=mix(vec3(1., 1., 1.), finalColor.xyz, smoothstep(0., 1., 1.-mixfactor)); 
	}
#ifdef IS_MOON
    if(final_illumination < 0.9999)
    {
        lowp vec4 shadowColor = texture2D(earthShadow, vec2(final_illumination, 0.5));
        finalColor =
		eclipsePush*(1.0-0.75*shadowColor.a)*
		mix(finalColor * litColor, shadowColor, clamp(shadowColor.a, 0.0, 0.7)); // clamp alpha to allow some maria detail.
    }
    else
#endif
    {
        finalColor *= litColor;
    }

    //apply white rimlight
    finalColor.xyz = clamp( finalColor.xyz + vec3(outgas), 0.0, 1.0);

    FRAG_COLOR = finalColor;
    //to debug texture issues, uncomment and reload shader
    //FRAG_COLOR = texColor;
}
