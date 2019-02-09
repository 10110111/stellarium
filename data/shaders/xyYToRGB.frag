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

vec3 adaptColorToVisionMode(vec3 rgb)
{
    vec3 color=XYZ2xyY(RGB2XYZ(rgb));
	///////////////////////////////////////////////////////////////////////////
	// Now we have the xyY components in color, need to convert to RGB

	// 1. Hue conversion
	// if log10Y>0.6, photopic vision only (with the cones, colors are seen)
	// else scotopic vision if log10Y<-2 (with the rods, no colors, everything blue),
	// else mesopic vision (with rods and cones, transition state)

    mediump vec3 resultColor;
	if (color[2] <= 0.01)
	{
		// special case for s = 0 (x=0.25, y=0.25)
		color[2] *= 0.5121445;
		color[2] = pow(abs(color[2]*pi*0.0001), alphaWaOverAlphaDa*oneOverGamma)* term2TimesOneOverMaxdLpOneOverGamma;
		color[0] = 0.787077*color[2];
		color[1] = 0.9898434*color[2];
		color[2] *= 1.9256125;
		resultColor = color.xyz*brightnessScale;
	}
	else
	{
		if (color[2]<3.9810717055349722)
		{
			// Compute s, ratio between scotopic and photopic vision
			float op = (log(color[2])/ln10 + 2.)/2.6;
			float s = op * op *(3. - 2. * op);
			// Do the blue shift for scotopic vision simulation (night vision) [3]
			// The "night blue" is x,y(0.25, 0.25)
			color[0] = (1. - s) * 0.25 + s * color[0];	// Add scotopic + photopic components
			color[1] = (1. - s) * 0.25 + s * color[1];	// Add scotopic + photopic components
			// Take into account the scotopic luminance approximated by V [3] [4]
			float V = color[2] * (1.33 * (1. + color[1] / color[0] + color[0] * (1. - color[0] - color[1])) - 1.68);
			color[2] = 0.4468 * (1. - s) * V + s * color[2];
		}

		// 2. Adapt the luminance value and scale it to fit in the RGB range [2]
		// color[2] = std::pow(adaptLuminanceScaled(color[2]), oneOverGamma);
		color[2] = pow(abs(color[2]*pi*0.0001), alphaWaOverAlphaDa*oneOverGamma)* term2TimesOneOverMaxdLpOneOverGamma;

		// Convert from xyY to XZY
		// Use a XYZ to Adobe RGB (1998) matrix which uses a D65 reference white
		mediump vec3 tmp = vec3(color[0] * color[2] / color[1], color[2], (1. - color[0] - color[1]) * color[2] / color[1]);
		resultColor = vec3(2.04148*tmp.x-0.564977*tmp.y-0.344713*tmp.z, -0.969258*tmp.x+1.87599*tmp.y+0.0415557*tmp.z, 0.0134455*tmp.x-0.118373*tmp.y+1.01527*tmp.z);
		resultColor*=brightnessScale;
	}
    return resultColor;
}
