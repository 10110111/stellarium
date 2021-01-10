/*
 * Stellarium
 * Copyright (C) 2020 Ruslan Kabatsayev
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

#ifndef ATMOSPHERE_SHOWMYSKY_HPP
#define ATMOSPHERE_SHOWMYSKY_HPP

#include "Atmosphere.hpp"
#include "VecMath.hpp"

#include "Skybright.hpp"

#include <array>
#include <memory>
#include <stdexcept>
#include <QLibrary>
#include <QOpenGLBuffer>

#include <ShowMySky/AtmosphereRenderer.hpp>


class StelProjector;
class StelToneReproducer;
class StelCore;
class QOpenGLFunctions;

class AtmosphereShowMySky : public Atmosphere
{
public:
	AtmosphereShowMySky();
	~AtmosphereShowMySky();
	
	void computeColor(double JD, Vec3d _sunPos, Vec3d moonPos, float moonPhase, float moonMagnitude,
					  float nonExtinctedLunarMagnitude, StelCore* core,
					  float latitude, float altitude, float temperature, float relativeHumidity) override;
	void draw(StelCore* core) override;

	struct InitFailure : std::runtime_error
	{
		using std::runtime_error::runtime_error;
		InitFailure(QString const& what) : std::runtime_error(what.toStdString()) {}
	};

private:
	QLibrary showMySkyLib;
	Vec4i viewport;
	int gridMaxY,gridMaxX;

	QVector<Vec2f> posGrid;
	QOpenGLBuffer posGridBuffer;
	QOpenGLBuffer indexBuffer;
	QVector<Vec4f> viewRayGrid;
	QOpenGLBuffer viewRayGridBuffer;

	std::unique_ptr<QOpenGLShaderProgram> luminanceToScreenProgram_;
	decltype(::ShowMySky_AtmosphereRenderer_create)* ShowMySky_AtmosphereRenderer_create=nullptr;

	struct {
		int rgbMaxValue;
		int bayerPattern;
		int oneOverGamma;
		int brightnessScale;
		int luminanceTexture;
		int alphaWaOverAlphaDa;
		int term2TimesOneOverMaxdL;
	} shaderAttribLocations;

	GLuint bayerPatternTex_=0;

	int prevWidth_=0, prevHeight_=0;
	GLuint vao_=0, vbo_=0;
	std::unique_ptr<ShowMySky::AtmosphereRenderer> renderer_;
	std::unique_ptr<ShowMySky::Settings> skySettings_;

	void resolveFunctions();
	void loadShaders();
	void setupBuffers();
	void regenerateGrid();
	void setupRenderTarget();
	// Gets average value of the pixels rendered to the FBO texture as the value of the deepest mipmap level
	Vec4f getMeanPixelValue(int texW, int texH);
	void resizeRenderTarget(int width, int height);
	void drawAtmosphere(Mat4f const& projectionMatrix, float sunAzimuth, float sunZenithAngle,
						float moonAzimuth, float moonZenithAngle, float altitude,
						float brightness, bool clearTarget);
};

#endif // ATMOSTPHERE_HPP
