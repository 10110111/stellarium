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
uniform mediump vec3 ambientLight;
uniform mediump vec3 diffuseLight;
uniform highp vec4 sunInfo;
uniform mediump float skyBrightness;

uniform int shadowCount;
uniform highp mat4 shadowData;

//eye direction in model space, pre-normalized
uniform highp vec3 eyeDirection;

//light direction in model space, pre-normalized
uniform highp vec3 lightDirection;
//x = A, y = B, z = scaling factor (rho/pi * E0), w roughness
uniform mediump vec4 orenNayarParameters;

uniform sampler2D earthShadow;
uniform mediump float eclipsePush;
uniform sampler2D normalMap;
uniform sampler2D horizonMap;

VARYING highp vec3 normalX;
VARYING highp vec3 normalY;
VARYING highp vec3 normalZ;

const highp float M_PI=3.1415926535897932384626433832795;

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

void main()
{
    mediump float final_illumination = 1.0;
    mediump float lum = 1.;
    //shadow calculation
    if(lum > 0.0)
    {
        highp vec3 sunPosition = sunInfo.xyz;
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
                    illumination = (d / (r - R)) * 0.594; // prepare texture coordinate. 0.6=umbra edge. Smaller number->larger shadow.
                }
                else if(d <= R - r)
                {
                    // penumbra completely inside (Annular, or a moon transit in front of the Sun.)
                    illumination = 1.0 - r * r / (R * R);
                }
                else // penumbra: partially inside
                {
                    //illumination = ((d - abs(R-r)) / (R + r - abs(R-r))) * 0.4 + 0.6;
                    illumination = ((d - r + R) / (2.0 * R )) * 0.406 + 0.594;
                }
                final_illumination = min(illumination, final_illumination);
            }
        }
    }

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
    // Use an Oren-Nayar model for rough surfaces
    // Ref: http://content.gpwiki.org/index.php/D3DBook:(Lighting)_Oren-Nayar
    lum = orenNayar(normal, lightDirection, eyeDirection, orenNayarParameters.x, orenNayarParameters.y, orenNayarParameters.z, orenNayarParameters.w);
	lum *= horizonShadowCoefficient;
//Reduce lum if sky is bright, to avoid burnt-out look in daylight sky.
    //lum *= (1.0-0.4*skyBrightness);
    lum *= clamp((0.9-0.05*log(skyBrightness)), 0.1, 0.9);

    //final lighting color
    mediump vec4 litColor = vec4(lum * final_illumination * diffuseLight + ambientLight, 1.0);

    lowp vec4 texColor = texture2D(tex, moonTexCoord);
    // Undo the extraneous gamma encoded in the texture.
    // FIXME: ideally, we want all the calculations to be done in linear scale,
    // and then *in the end* apply the exact sRGB transfer function. Currently
    // though, we don't do actual physically-correct simulation here, so we
    // just apply the approximate sRGB gamma.
    texColor = pow(texColor, vec4(2.8/2.2));

    mediump vec4 finalColor = texColor;
    if(final_illumination < 0.9999)
    {
        lowp vec4 shadowColor = texture2D(earthShadow, vec2(final_illumination, 0.5));
        finalColor =
		eclipsePush*(1.0-0.75*shadowColor.a)*
		mix(finalColor * litColor, shadowColor, clamp(shadowColor.a, 0.0, 0.7)); // clamp alpha to allow some maria detail.
    }
    else
    {
        finalColor *= litColor;
    }

    FRAG_COLOR = finalColor;
}

