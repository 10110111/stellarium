/*
 * Stellarium
 * Copyright (C) 2019 Ruslan Kabatsayev
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

#ifndef ATMOSPHERE_BRUNETON_HPP
#define ATMOSPHERE_BRUNETON_HPP

#include "Atmosphere.hpp"
#include "VecMath.hpp"

#include "Skybright.hpp"
#include "StelFader.hpp"

#include <array>
#include <memory>
#include <QOpenGLBuffer>

class StelProjector;
class StelToneReproducer;
class StelCore;
class QOpenGLFunctions;

//! Compute and display the daylight sky color using openGL.
//! The sky brightness is computed with the SkyBright class, the color with the SkyLight.
//! Don't use this class directly but use it through the LandscapeMgr.
class AtmosphereBruneton : public Atmosphere
{
public:
	AtmosphereBruneton();
	~AtmosphereBruneton();
	
	void computeColor(double JD, Vec3d _sunPos, Vec3d moonPos, float moonPhase, float moonMagnitude, StelCore* core,
					  float latitude, float altitude, float temperature, float relativeHumidity) override;
	void draw(StelCore* core) override;
	void update(double deltaTime) override {fader.update((int)(deltaTime*1000));}

	void setFadeDuration(float duration) override {fader.setDuration((int)(duration*1000.f));}
	float getFadeDuration() const override {return (float)fader.getDuration()/1000.f;}

	void setFlagShow(bool b) override {fader = b;}
	bool getFlagShow() const override  {return fader;}

	float getRealDisplayIntensityFactor() const override  {return fader.getInterstate()*eclipseFactor;}

	float getFadeIntensity() const override  {return fader.getInterstate();}

	float getAverageLuminance() const override  {return averageLuminance;}

	void setAverageLuminance(float overrideLum) override;
	void setLightPollutionLuminance(float f) override { lightPollutionLuminance = f; }
	float getLightPollutionLuminance() const override { return lightPollutionLuminance; }

private:
	Vec4i viewport;
	Skybright skyb;
	int gridMaxY,gridMaxX;

	QVector<Vec2f> posGrid;
	QOpenGLBuffer posGridBuffer;
	QOpenGLBuffer indexBuffer;
	QVector<Vec4f> viewRayGrid;
	QOpenGLBuffer viewRayGridBuffer;

	//! The average luminance of the atmosphere in cd/m2
	float averageLuminance;
	bool overrideAverageLuminance; // if true, don't compute but keep value set via setAverageLuminance(float)
	float eclipseFactor;
	LinearFader fader;
	float lightPollutionLuminance;

	//! Vertex shader used for xyYToRGB computation
	std::unique_ptr<QOpenGLShaderProgram> atmoShaderProgram;
	struct {
		int bayerPattern;
		int rgbMaxValue;
		int alphaWaOverAlphaDa;
		int oneOverGamma;
		int term2TimesOneOverMaxdLpOneOverGamma;
		int brightnessScale;
		int sunDir;
		int cameraPos;
		int projectionMatrix;
		int skyVertex;
		int viewRay;
		int transmittanceTexture;
		int scatteringTexture;
		int irradianceTexture;
		int singleMieScatteringTexture;
	} shaderAttribLocations;

	GLuint bayerPatternTex=0;

	enum
	{
		TRANSMITTANCE_TEXTURE,
		SCATTERING_TEXTURE,
		IRRADIANCE_TEXTURE,
		MIE_SCATTERING_TEXTURE,

		TEX_COUNT
	};
	GLuint textures[TEX_COUNT];
	Vec3d sunDir;
	double altitude;
	void loadTextures();
	void loadShaders();
	void regenerateGrid();
	void updateEclipseFactor(StelCore* core, Vec3d sunPos, Vec3d moonPos);
	bool separateMieTexture=false;

	struct TextureSize4D
	{
		std::array<std::int32_t,4> sizes;
		int muS_size() const { return sizes[0]; }
		int mu_size() const { return sizes[1]; }
		int nu_size() const { return sizes[2]; }
		int r_size() const { return sizes[3]; }
		int width() const { return muS_size(); }
		int height() const { return mu_size(); }
		int depth() const { return nu_size()*r_size(); }
		bool operator!=(TextureSize4D const& rhs) const { return sizes!=rhs.sizes; }
	};
	struct TextureSize2D
	{
		std::array<std::int32_t,2> sizes;
		int width() const { return sizes[0]; }
		int height() const { return sizes[1]; }
	};
	TextureSize2D transmittanceTextureSize, irradianceTextureSize;
	TextureSize4D scatteringTextureSize, mieScatteringTextureSize;

	static TextureSize4D getTextureSize4D(QVector<char> const& data);
	static TextureSize2D getTextureSize2D(QVector<char> const& data);
};

#endif // ATMOSTPHERE_HPP
