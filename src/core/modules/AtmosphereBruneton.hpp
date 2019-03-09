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
#include <stdexcept>
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
	
	void computeColor(double JD, Vec3d _sunPos, Vec3d moonPos, float moonPhase, float moonMagnitude, float nonExtinctedLunarMagnitude, StelCore* core,
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

	struct InitFailure : std::runtime_error
	{
		using std::runtime_error::runtime_error;
	};

private:
	Vec4i viewport;
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

	std::unique_ptr<QOpenGLShaderProgram> atmosphereRenderProgram;
	std::unique_ptr<QOpenGLShaderProgram> postProcessProgram;
	std::unique_ptr<QOpenGLShaderProgram> roundTextureSizeProgram; // for getMeanPixelValueWithWorkaround()
	bool needNPOTMipmapWorkaround=false;
	struct {
		int sunDir;
		int cameraPos;
		int projectionMatrix;
		int scatteringTexture;
		int irradianceTexture;
		int transmittanceTexture;
		int singleMieScatteringTexture;
		int viewRayAndAddLuminance; // View direction int xyz components, additional luminance in w
		int skyVertex;

		int rgbMaxValue;
		int bayerPattern;
		int oneOverGamma;
		int radianceImage;
		int brightnessScale;
		int alphaWaOverAlphaDa;
		int term2TimesOneOverMaxdLpOneOverGamma;

		int textureToRound=-1; // for getMeanPixelValueWithWorkaround()
	} shaderAttribLocations;

	GLuint bayerPatternTex=0;

	enum
	{
		TRANSMITTANCE_TEXTURE,
		SCATTERING_TEXTURE,
		IRRADIANCE_TEXTURE,
		MIE_SCATTERING_TEXTURE,
		FBO_TEXTURE,
		FBO_TEXTURE_POT, // for getMeanPixelValueWithWorkaround()

		TEX_COUNT
	};
	GLuint textures[TEX_COUNT];
	enum
	{
		FBO_MAIN,
		FBO_POT, // for getMeanPixelValueWithWorkaround()

		FBO_COUNT
	};
	GLuint fbos[FBO_COUNT];
	int fboPrevWidth, fboPrevHeight;
	double altitude;
	GLuint vao, vbo;

	void(QOPENGLF_APIENTRYP TexImage3D)(GLenum,GLint,GLint,GLsizei,
										GLsizei,GLsizei,GLint,GLenum,
										GLenum,const GLvoid*);
	void(QOPENGLF_APIENTRYP GetTexImage)(GLenum,GLint,GLenum,GLenum,GLvoid*);
	void(QOPENGLF_APIENTRYP GenVertexArrays)(GLsizei,GLuint*);
	void(QOPENGLF_APIENTRYP BindVertexArray)(GLuint);
	void(QOPENGLF_APIENTRYP DeleteVertexArrays)(GLsizei,const GLuint*);

	void checkNeedForNPOTMipmapWorkaround();
	void resolveFunctions();
	void loadShaders();
	void loadTextures();
	void setupBuffers();
	void regenerateGrid();
	void setupRenderTarget();
	static QVector<char> readFileData(QString const& path);
	// Gets average value of the pixels rendered to the FBO texture as the value of the deepest
	// mipmap level
	Vec4f getMeanPixelValue(int texW, int texH);
	//! \brief Functionally the same as @ref getMeanPixelValue, but works around Mesa bug 109816.
	//  It's less precise, so isn't used for all OpenGL implementations.
	Vec4f getMeanPixelValueWithWorkaround(int texW, int texH);
	void resizeRenderTarget(int width, int height);
	void drawAtmosphere(Mat4f const& projectionMatrix, Vec3d sunDir, float brightness);
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
