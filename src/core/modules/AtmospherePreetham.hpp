/*
 * Stellarium
 * Copyright (C) 2003 Fabien Chereau
 * Copyright (C) 2012 Timothy Reaves
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

#ifndef ATMOSPHERE_PREETHAM_HPP
#define ATMOSPHERE_PREETHAM_HPP

#include "Atmosphere.hpp"
#include "Skylight.hpp"
#include "VecMath.hpp"

#include "Skybright.hpp"
#include "StelFader.hpp"

#include <QOpenGLBuffer>

class StelProjector;
class StelToneReproducer;
class StelCore;

//! Compute and display the daylight sky color using openGL.
//! The sky brightness is computed with the SkyBright class, the color with the SkyLight.
//! Don't use this class directly but use it through the LandscapeMgr.
class AtmospherePreetham : public Atmosphere
{
public:
	AtmospherePreetham();
	virtual ~AtmospherePreetham();
	
	//! Compute sky brightness values and average luminance.
	void computeColor(double JD, Vec3d _sunPos, Vec3d moonPos, float moonPhase, float moonMagnitude, float nonExtinctedLunarMagnitude, StelCore* core,
		float latitude = 45.f, float altitude = 200.f,
		float temperature = 15.f, float relativeHumidity = 40.f);
	void draw(StelCore* core);
	void update(double deltaTime) {fader.update((int)(deltaTime*1000));}

private:
	Vec4i viewport;
	Skylight sky;
	Skybright skyb;
	int skyResolutionY,skyResolutionX;

	Vec2f* posGrid;
	QOpenGLBuffer posGridBuffer;
	QOpenGLBuffer indicesBuffer;
	Vec4f* colorGrid;
	QOpenGLBuffer colorGridBuffer;

	//! Vertex shader used for xyYToRGB computation
	class QOpenGLShaderProgram* atmoShaderProgram;
	struct {
		int bayerPattern;
		int rgbMaxValue;
		int alphaWaOverAlphaDa;
		int oneOverGamma;
		int term2TimesOneOverMaxdL;
		int brightnessScale;
		int sunPos;
		int term_x, Ax, Bx, Cx, Dx, Ex;
		int term_y, Ay, By, Cy, Dy, Ey;
		int projectionMatrix;
		int skyVertex;
		int skyColor;
	} shaderAttribLocations;

	GLuint bayerPatternTex=0;
};

#endif // ATMOSPHERE_PREETHAM_HPP
