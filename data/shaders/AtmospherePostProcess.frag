/*
 * Stellarium
 * Copyright (C) 2002-2018 Fabien Chereau and Stellarium contributors
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

#version 130

const highp float pi = 3.1415926535897931;
const highp float ln10 = 2.3025850929940459;

// Variable for the xyYTo RGB conversion
uniform highp float alphaWaOverAlphaDa;
uniform highp float oneOverGamma;
uniform highp float term2TimesOneOverMaxdLpOneOverGamma;
uniform highp float brightnessScale;

uniform sampler2D radiance;

vec3 XYZ2xyY(vec3 XYZ)
{
	return vec3(XYZ.xy/(XYZ.x+XYZ.y+XYZ.z), XYZ.y);
}

vec3 RGB2XYZ(vec3 rgb)
{
	return rgb*mat3(0.4124564, 0.3575761, 0.1804375,
					0.2126729, 0.7151522, 0.0721750,
					0.0193339, 0.1191920, 0.9503041);
}

in vec2 texCoord;
out vec3 resultColor;

vec3 dither(vec3);
void main()
{
	vec3 rgb=texture(radiance, texCoord).rgb;
	vec3 xyY=XYZ2xyY(RGB2XYZ(rgb));
	float x=xyY.x, y=xyY.y, Y=xyY.z;
	///////////////////////////////////////////////////////////////////////////
	// Now we have the xyY components, need to convert to light-adapted RGB

	// 1. Hue conversion
	// if log10Y>0.6, photopic vision only (with the cones, colors are seen)
	// else scotopic vision if log10Y<-2 (with the rods, no colors, everything blue),
	// else mesopic vision (with rods and cones, transition state)

	if (Y <= 0.01)
	{
		// special case for s = 0 (x=0.25, y=0.25)
		Y *= 0.5121445;
		Y = pow(abs(Y*pi*0.0001), alphaWaOverAlphaDa*oneOverGamma)* term2TimesOneOverMaxdLpOneOverGamma;
		resultColor = vec3(0.787077, 0.9898434, 1.9256125)*Y*brightnessScale;
	}
	else
	{
		if (Y<3.9810717055349722)
		{
			// Compute s, ratio between scotopic and photopic vision
			float op = (log(Y)/ln10 + 2.)/2.6;
			float s = op * op *(3. - 2. * op);
			// Do the blue shift for scotopic vision simulation (night vision) [3]
			// The "night blue" is x,y(0.25, 0.25)
			x = (1. - s) * 0.25 + s * x;	// Add scotopic + photopic components
			y = (1. - s) * 0.25 + s * y;	// Add scotopic + photopic components
			// Take into account the scotopic luminance approximated by V [3] [4]
			float V = Y * (1.33 * (1. + y / x + x * (1. - x - y)) - 1.68);
			Y = 0.4468 * (1. - s) * V + s * Y;
		}

		// 2. Adapt the luminance value and scale it to fit in the RGB range [2]
		// Y = std::pow(adaptLuminanceScaled(Y), oneOverGamma);
		Y = pow(abs(Y*pi*0.0001), alphaWaOverAlphaDa*oneOverGamma)* term2TimesOneOverMaxdLpOneOverGamma;

		// Convert from xyY to XZY
		// Use a XYZ to Adobe RGB (1998) matrix which uses a D65 reference white
		mediump vec3 tmp = vec3(x * Y / y, Y, (1. - x - y) * Y / y);
		resultColor = vec3(2.04148*tmp.x-0.564977*tmp.y-0.344713*tmp.z, -0.969258*tmp.x+1.87599*tmp.y+0.0415557*tmp.z, 0.0134455*tmp.x-0.118373*tmp.y+1.01527*tmp.z);
		resultColor*=brightnessScale;
	}
    resultColor=dither(resultColor);
}
