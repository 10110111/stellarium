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

#include "AtmosphereShowMySky.hpp"
#include "StelUtils.hpp"
#include "StelApp.hpp"
#include "StelProjector.hpp"
#include "StelToneReproducer.hpp"
#include "StelCore.hpp"
#include "StelPainter.hpp"
#include "StelFileMgr.hpp"
#include "Dithering.hpp"
#include "StelTranslator.hpp"

#include <cassert>
#include <cstring>

#include <QDir>
#include <QFile>
#include <QDebug>
#include <QVector>
#include <QSettings>
#include <QOpenGLShaderProgram>
#include <QOpenGLFunctions_3_3_Core>

namespace
{

constexpr float nonExtinctedSolarMagnitude=-26.74;

struct SkySettings : ShowMySky::Settings, QObject
{
	double altitude() override { return altitude_; }
	double sunAzimuth() override { return sunAzimuth_; }
	double sunZenithAngle() override { return sunZenithAngle_; }
	double moonAzimuth() override { return moonAzimuth_; }
	double moonZenithAngle() override { return moonZenithAngle_; }

	bool zeroOrderScatteringEnabled() override { init(); return zeroOrderScatteringEnabled_; }
	bool singleScatteringEnabled   () override { init(); return singleScatteringEnabled_; }
	bool multipleScatteringEnabled () override { init(); return multipleScatteringEnabled_; }

	bool onTheFlySingleScatteringEnabled       () override { init(); return false; }
	bool onTheFlyPrecompDoubleScatteringEnabled() override { init(); return false; }

	bool usingEclipseShader() override
	{
		init();
		return false; // TODO: make it dynamic
	}

	double altitude_ = 0;
	double sunAzimuth_ = 0;
	double sunZenithAngle_ = M_PI/4;
	double moonAzimuth_ = 0;
	double moonZenithAngle_ = -M_PI/2;

	bool zeroOrderScatteringEnabled_ = false;
	bool singleScatteringEnabled_ = true;
	bool multipleScatteringEnabled_ = true;


	// These variables are not used by AtmosphereRenderer, but it's convenient to keep them here
	QMatrix4x4 projectionMatrix_;

	static constexpr const char* zeroOrderScatteringPropName() { return "LandscapeMgr.flagAtmosphereZeroOrderScattering"; }
	static constexpr const char* singleScatteringPropName()    { return "LandscapeMgr.flagAtmosphereSingleScattering"; }
	static constexpr const char* multipleScatteringPropName()  { return "LandscapeMgr.flagAtmosphereMultipleScattering"; }
	bool inited=false; // We have to lazy initialize because at the time of construction, LandscapeMgr isn't registered with StelPropertyMgr
	void init()
	{
		if(inited) return;

		{
			const auto prop = StelApp::getInstance().getStelPropertyManager()->getProperty(zeroOrderScatteringPropName());
			QObject::connect(prop, &StelProperty::changed, this, [this](const QVariant& v) { zeroOrderScatteringEnabled_=v.toBool(); });
			zeroOrderScatteringEnabled_ = prop->getValue().toBool();
		}
		{
			const auto prop = StelApp::getInstance().getStelPropertyManager()->getProperty(singleScatteringPropName());
			QObject::connect(prop, &StelProperty::changed, this, [this](const QVariant& v) { singleScatteringEnabled_=v.toBool(); });
			singleScatteringEnabled_ = prop->getValue().toBool();
		}
		{
			const auto prop = StelApp::getInstance().getStelPropertyManager()->getProperty(multipleScatteringPropName());
			QObject::connect(prop, &StelProperty::changed, this, [this](const QVariant& v) { multipleScatteringEnabled_=v.toBool(); });
			multipleScatteringEnabled_ = prop->getValue().toBool();
		}
		inited=true;
	}
};

QOpenGLFunctions_3_3_Core& glfuncs()
{
	return *QOpenGLContext::currentContext()->versionFunctions<QOpenGLFunctions_3_3_Core>();
}

}

void AtmosphereShowMySky::loadShaders()
{
	const auto handleCompileStatus=[](bool success, QOpenGLShader const& shader, const char* what)
		{
			if(!success)
			{
				qCritical("Error while compiling %s: %s", what, shader.log().toLatin1().constData());
				throw InitFailure("Shader compilation failed");
			}
			if(!shader.log().isEmpty())
			{
				qWarning("Warnings while compiling %s: %s", what, shader.log().toLatin1().constData());
			}
		};

	// Shader program that converts XYZW texture to sRGB image
	{
		static constexpr char vShaderSrc[]=R"(
#version 130
in vec4 vertex;
out vec2 texCoord;
void main()
{
	gl_Position=vertex;
	texCoord=vertex.xy*0.5+0.5;
}
)";
	static constexpr char fShaderSrc[]=R"(
#version 330

uniform sampler2D luminanceXYZW;

in vec2 texCoord;
out vec4 color;

vec3 dither(vec3);
vec3 xyYToRGB(float x, float y, float Y);

void main()
{
	vec3 XYZ=texture(luminanceXYZW, texCoord).xyz;
	vec3 srgb=xyYToRGB(XYZ.x/(XYZ.x+XYZ.y+XYZ.z), XYZ.y/(XYZ.x+XYZ.y+XYZ.z), XYZ.y);
	color=vec4(dither(srgb),1);
}
)";

		QOpenGLShader vShader(QOpenGLShader::Vertex);
		handleCompileStatus(vShader.compileSourceCode(vShaderSrc), vShader,
							"ShowMySky atmosphere luminance-to-screen vertex shader");
		luminanceToScreenProgram_->addShader(&vShader);

		QOpenGLShader ditherShader(QOpenGLShader::Fragment);
		handleCompileStatus(ditherShader.compileSourceCode(makeDitheringShader()), ditherShader,
							"ShowMySky atmosphere dithering shader");
		luminanceToScreenProgram_->addShader(&ditherShader);

		QOpenGLShader toneReproducerShader(QOpenGLShader::Fragment);
		handleCompileStatus(toneReproducerShader.compileSourceFile(":/shaders/xyYToRGB.glsl"),
							toneReproducerShader, "ShowMySky atmosphere tone reproducer fragment shader");
		luminanceToScreenProgram_->addShader(&toneReproducerShader);

		QOpenGLShader luminanceToScreenShader(QOpenGLShader::Fragment);
		handleCompileStatus(luminanceToScreenShader.compileSourceCode(fShaderSrc), luminanceToScreenShader,
							"ShowMySky atmosphere luminance-to-screen fragment shader");
		luminanceToScreenProgram_->addShader(&luminanceToScreenShader);

		if(!StelPainter::linkProg(luminanceToScreenProgram_.get(), "atmosphere luminance-to-screen"))
			throw InitFailure("Shader program linking failed");
	}

	{
		static constexpr char viewDirVertShaderSrc[]=R"(
#version 330
uniform mat4 projectionMatrix;

in vec2 skyVertex;
in vec4 viewRay;

out vec3 viewDirUnnormalized;

void main()
{
	viewDirUnnormalized = viewRay.xyz;
	gl_Position = projectionMatrix*vec4(skyVertex, 0., 1.);
}
)";
		static constexpr char viewDirFragShaderSrc[]=R"(
#version 330
in vec3 viewDirUnnormalized;
vec3 calcViewDir()
{
	return normalize(viewDirUnnormalized);
}
)";

		renderer_->loadData(viewDirVertShaderSrc, viewDirFragShaderSrc);
	}
}

void AtmosphereShowMySky::resizeRenderTarget(int width, int height)
{
	renderer_->resizeEvent(width, height);

	prevWidth_=width;
	prevHeight_=height;
}

void AtmosphereShowMySky::setupRenderTarget()
{
	auto& gl=glfuncs();

	GLint viewport[4];
	gl.glGetIntegerv(GL_VIEWPORT, viewport);
	resizeRenderTarget(viewport[2], viewport[3]);
}

void AtmosphereShowMySky::setupBuffers()
{
	auto& gl=glfuncs();

	gl.glGenVertexArrays(1, &vao_);
	gl.glBindVertexArray(vao_);
	gl.glGenBuffers(1, &vbo_);
	gl.glBindBuffer(GL_ARRAY_BUFFER, vbo_);
	const GLfloat vertices[]=
	{
		-1, -1,
		 1, -1,
		-1,  1,
		 1,  1,
	};
	gl.glBufferData(GL_ARRAY_BUFFER, sizeof vertices, vertices, GL_STATIC_DRAW);
	constexpr GLuint attribIndex=0;
	constexpr int coordsPerVertex=2;
	gl.glVertexAttribPointer(attribIndex, coordsPerVertex, GL_FLOAT, false, 0, 0);
	gl.glEnableVertexAttribArray(attribIndex);
	gl.glBindVertexArray(0);
}

void AtmosphereShowMySky::resolveFunctions()
{
	if(!showMySkyLib.load())
		throw InitFailure(q_("Failed to load ShowMySky library"));
	ShowMySky_AtmosphereRenderer_create=reinterpret_cast<decltype(ShowMySky_AtmosphereRenderer_create)>(
												showMySkyLib.resolve("ShowMySky_AtmosphereRenderer_create"));
	if(!ShowMySky_AtmosphereRenderer_create)
		throw InitFailure(q_("Failed to resolve the function to create AtmosphereRenderer"));
}

AtmosphereShowMySky::AtmosphereShowMySky()
	: showMySkyLib("ShowMySky")
	, viewport(0,0,0,0)
	, gridMaxY(44)
	, gridMaxX(44)
	, posGridBuffer(QOpenGLBuffer::VertexBuffer)
	, indexBuffer(QOpenGLBuffer::IndexBuffer)
	, viewRayGridBuffer(QOpenGLBuffer::VertexBuffer)
	, luminanceToScreenProgram_(new QOpenGLShaderProgram())
{
	setFadeDuration(1.5);

	resolveFunctions();
	try
	{
		const auto& conf=*StelApp::getInstance().getSettings();
		const auto defaultPath = QDir::homePath() + "/cms";
		const auto pathToData = conf.value("landscape/atmosphere_model_path", defaultPath).toString();
		auto& gl=glfuncs();
		qDebug() << "Will load CalcMySky atmosphere model from" << pathToData;
		skySettings_.reset(new SkySettings);
		const std::function<void(QOpenGLShaderProgram&)> renderSurface=[this](QOpenGLShaderProgram& prog)
			{
				const auto& settings = *static_cast<SkySettings*>(skySettings_.get());
				prog.setUniformValue("projectionMatrix", settings.projectionMatrix_);

				const auto viewRayLoc = prog.attributeLocation("viewRay");
				const auto skyVertexLoc = prog.attributeLocation("skyVertex");

				viewRayGridBuffer.bind();
				prog.setAttributeBuffer(viewRayLoc, GL_FLOAT, 0, 4, 0);
				viewRayGridBuffer.release();
				prog.enableAttributeArray(viewRayLoc);

				posGridBuffer.bind();
				prog.setAttributeBuffer(skyVertexLoc, GL_FLOAT, 0, 2, 0);
				posGridBuffer.release();
				prog.enableAttributeArray(skyVertexLoc);

				auto& gl=glfuncs();

				indexBuffer.bind();
				std::size_t shift=0;
				for (int y=0;y<gridMaxY;++y)
				{
					gl.glDrawElements(GL_TRIANGLE_STRIP, (gridMaxX+1)*2, GL_UNSIGNED_SHORT,
									  reinterpret_cast<void*>(shift));
					shift += (gridMaxX+1)*2*2;
				}
				indexBuffer.release();

				prog.disableAttributeArray(skyVertexLoc);
				prog.disableAttributeArray(viewRayLoc);
				prog.release();
			};
		renderer_.reset(ShowMySky_AtmosphereRenderer_create(&gl, &pathToData, skySettings_.get(), &renderSurface));
		loadShaders();
		setupRenderTarget();
		setupBuffers();
	}
	catch(ShowMySky::Error const& error)
	{
		throw InitFailure(error.what());
	}

	{
		auto& prog=*luminanceToScreenProgram_;
		prog.bind();

		shaderAttribLocations.rgbMaxValue        = prog.uniformLocation("rgbMaxValue");
		shaderAttribLocations.bayerPattern       = prog.uniformLocation("bayerPattern");
		shaderAttribLocations.oneOverGamma       = prog.uniformLocation("oneOverGamma");
		shaderAttribLocations.brightnessScale    = prog.uniformLocation("brightnessScale");
		shaderAttribLocations.luminanceTexture   = prog.uniformLocation("luminance");
		shaderAttribLocations.alphaWaOverAlphaDa = prog.uniformLocation("alphaWaOverAlphaDa");
		shaderAttribLocations.term2TimesOneOverMaxdL
												 = prog.uniformLocation("term2TimesOneOverMaxdL");

		prog.release();
	}
}

AtmosphereShowMySky::~AtmosphereShowMySky()
{
	if(auto*const ctx=QOpenGLContext::currentContext())
	{
		auto& gl=glfuncs();
		gl.glDeleteBuffers(1, &vbo_);
		gl.glDeleteVertexArrays(1, &vao_);
	}
}

void AtmosphereShowMySky::regenerateGrid()
{
	const float width=viewport[2], height=viewport[3];
	gridMaxY = StelApp::getInstance().getSettings()->value("landscape/atmosphereybin", 44).toInt();
	gridMaxX = std::floor(0.5+gridMaxY*(0.5*std::sqrt(3.0))*width/height);
	const auto gridSize=(1+gridMaxX)*(1+gridMaxY);
	posGrid.resize(gridSize);
	viewRayGrid.resize(gridSize);
	const auto stepX = width / (gridMaxX-0.5);
	const auto stepY = height / gridMaxY;
	const float viewportLeft = viewport[0];
	const float viewportBottom = viewport[1];
	for(int y=0; y<=gridMaxY; ++y)
	{
		for(int x=0; x<=gridMaxX; ++x)
		{
			Vec2f& v=posGrid[y*(1+gridMaxX)+x];
			v[0] = viewportLeft + (x == 0 ? 0
										  : x == gridMaxX ? width
														  : (x-0.5*(y&1))*stepX);
			v[1] = viewportBottom+y*stepY;
		}
	}
	posGridBuffer.destroy();
	posGridBuffer.setUsagePattern(QOpenGLBuffer::StaticDraw);
	posGridBuffer.create();
	posGridBuffer.bind();
	posGridBuffer.allocate(&posGrid[0], posGrid.size()*sizeof posGrid[0]);
	posGridBuffer.release();

	// Generate the indices used to draw the quads
	QVector<GLushort> indices((gridMaxX+1)*gridMaxY*2);
	for(int y=0, i=0; y<gridMaxY; ++y)
	{
		auto g0 = y*(1+gridMaxX);
		auto g1 = (y+1)*(1+gridMaxX);
		for(int x=0; x<=gridMaxX; ++x)
		{
			indices[i++]=g0++;
			indices[i++]=g1++;
		}
	}
	indexBuffer.destroy();
	indexBuffer.setUsagePattern(QOpenGLBuffer::StaticDraw);
	indexBuffer.create();
	indexBuffer.bind();
	indexBuffer.allocate(&indices[0], (gridMaxX+1)*gridMaxY*2*2);
	indexBuffer.release();

	viewRayGridBuffer.destroy();
	viewRayGridBuffer.setUsagePattern(QOpenGLBuffer::DynamicDraw);
	viewRayGridBuffer.create();
	viewRayGridBuffer.bind();
	viewRayGridBuffer.allocate(&viewRayGrid[0], (1+gridMaxX)*(1+gridMaxY)*4*4);
	viewRayGridBuffer.release();
}

void AtmosphereShowMySky::drawAtmosphere(Mat4f const& projectionMatrix, const float sunAzimuth, const float sunZenithAngle,
										 const float moonAzimuth, const float moonZenithAngle, const float altitude,
										 const float brightness, const bool clearTarget)
{
	const auto& m = projectionMatrix;
	auto& settings = *static_cast<SkySettings*>(skySettings_.get());
	settings.projectionMatrix_ = QMatrix4x4(m[0], m[4], m[8] , m[12],
											m[1], m[5], m[9] , m[13],
											m[2], m[6], m[10], m[14],
											m[3], m[7], m[11], m[15]);
	settings.altitude_=altitude;
	settings.sunAzimuth_=sunAzimuth;
	settings.sunZenithAngle_=sunZenithAngle;
	settings.moonAzimuth_=moonAzimuth;
	settings.moonZenithAngle_=moonZenithAngle;

	renderer_->draw(brightness, clearTarget);
}

Vec4f AtmosphereShowMySky::getMeanPixelValue(int texW, int texH)
{
	auto& gl=glfuncs();

	gl.glActiveTexture(GL_TEXTURE0);
	gl.glBindTexture(GL_TEXTURE_2D, renderer_->getLuminanceTexture());
	gl.glGenerateMipmap(GL_TEXTURE_2D);

	using namespace std;
	// Formula from the glspec, "Mipmapping" subsection in section 3.8.11 Texture Minification
	const auto totalMipmapLevels = 1+floor(log2(max(texW,texH)));
	const auto deepestLevel=totalMipmapLevels-1;

#ifndef NDEBUG
	// Sanity check
	int deepestMipmapLevelWidth=-1, deepestMipmapLevelHeight=-1;
	gl.glGetTexLevelParameteriv(GL_TEXTURE_2D, deepestLevel, GL_TEXTURE_WIDTH, &deepestMipmapLevelWidth);
	gl.glGetTexLevelParameteriv(GL_TEXTURE_2D, deepestLevel, GL_TEXTURE_HEIGHT, &deepestMipmapLevelHeight);
	assert(deepestMipmapLevelWidth==1);
	assert(deepestMipmapLevelHeight==1);
#endif

	Vec4f pixel;
	gl.glGetTexImage(GL_TEXTURE_2D, deepestLevel, GL_RGBA, GL_FLOAT, &pixel[0]);
	return pixel;
}

void AtmosphereShowMySky::computeColor(const double JD, Vec3d sunPos, Vec3d moonPos, const float moonPhase,
									   const float moonMagnitude, const float nonExtinctedLunarMagnitude, StelCore*const core,
									   const float latitude, const float altitude, const float temperature, const float relativeHumidity)
{
	// The majority of calculations is done in fragment shader, but we still need a nontrivial
	// grid to pass view rays, corresponding to the chosen projection, to the shader. Of course,
	// best quality results would be if we did projection inside the shader, but that'd require
	// having two implementations for each projection type: one for CPU and one for GPU.
	const auto prj = core->getProjection(StelCore::FrameAltAz, StelCore::RefractionOff);
	if (viewport != prj->getViewport())
	{
		viewport = prj->getViewport();
		regenerateGrid();
	}
	const auto width=viewport[2], height=viewport[3];
	if(width!=prevWidth_ || height!=prevHeight_)
		resizeRenderTarget(width, height);

	if (std::isnan(sunPos.length()))
		sunPos.set(0.,0.,-1.*AU);
	if (std::isnan(moonPos.length()))
		moonPos.set(0.,0.,-1.*AU);

	// FIXME: ignoring eclipse factor, it remains =1

	// No need to calculate if not visible
	if (!fader.getInterstate())
	{
		// GZ 20180114: Why did we add light pollution if atmosphere was not visible?????
		// And what is the meaning of 0.001? Approximate contribution of stellar background? Then why is it 0.0001 below???
		averageLuminance = 0.001f;
		return;
	}


	// FIXME: ignoring the "additional luminance" like star background etc.; see AtmospherePreetham for all potentially needed terms
	const auto numViewRayGridPoints=(1+gridMaxX)*(1+gridMaxY);
	for (int i=0; i<numViewRayGridPoints; ++i)
	{
		Vec3d point(1, 0, 0);
		prj->unProject(posGrid[i][0],posGrid[i][1],point);

		viewRayGrid[i].set(point[0], point[1], point[2], 0);
	}

	viewRayGridBuffer.bind();
	viewRayGridBuffer.write(0, &viewRayGrid[0], viewRayGrid.size()*sizeof viewRayGrid[0]);
	viewRayGridBuffer.release();


	const auto sunDir = sunPos / sunPos.length();
	const auto moonDir = moonPos / moonPos.length();

	const auto sunAzimuth = std::atan2(sunDir[1], sunDir[0]);
	const auto sunZenithAngle = std::acos(sunDir[2]);
	const auto moonAzimuth = std::atan2(moonDir[1], moonDir[0]);
	const auto moonZenithAngle = std::acos(moonDir[2]);

	drawAtmosphere(prj->getProjectionMatrix(), sunAzimuth, sunZenithAngle, moonAzimuth, moonZenithAngle,  altitude, 1.0, true);
	const auto moonRelativeBrightness=std::pow(10.f,0.4f*(nonExtinctedSolarMagnitude-nonExtinctedLunarMagnitude));
	drawAtmosphere(prj->getProjectionMatrix(), moonAzimuth, moonZenithAngle, 0, M_PI, altitude, moonRelativeBrightness, false);

	if (!overrideAverageLuminance)
	{
		const auto meanPixelValue=getMeanPixelValue(width, height);
		const auto meanY=meanPixelValue[1];

		averageLuminance = meanY;
	}
}

void AtmosphereShowMySky::draw(StelCore* core)
{
	if (StelApp::getInstance().getVisionModeNight())
		return;

	StelToneReproducer* eye = core->getToneReproducer();

	if (!fader.getInterstate())
		return;

	const float atm_intensity = fader.getInterstate();

	luminanceToScreenProgram_->bind();
	float a, b, c;
	eye->getShadersParams(a, b, c);
	luminanceToScreenProgram_->setUniformValue(shaderAttribLocations.alphaWaOverAlphaDa, a);
	luminanceToScreenProgram_->setUniformValue(shaderAttribLocations.oneOverGamma, b);
	luminanceToScreenProgram_->setUniformValue(shaderAttribLocations.term2TimesOneOverMaxdL, c);
	luminanceToScreenProgram_->setUniformValue(shaderAttribLocations.brightnessScale, atm_intensity);

	StelPainter sPainter(core->getProjection2d());
	sPainter.setBlending(true, GL_ONE, GL_ONE);

	const auto rgbMaxValue=calcRGBMaxValue(sPainter.getDitheringMode());
	luminanceToScreenProgram_->setUniformValue(shaderAttribLocations.rgbMaxValue, rgbMaxValue[0], rgbMaxValue[1], rgbMaxValue[2]);

	auto& gl=glfuncs();
	gl.glActiveTexture(GL_TEXTURE0);
	gl.glBindTexture(GL_TEXTURE_2D, renderer_->getLuminanceTexture());
	luminanceToScreenProgram_->setUniformValue(shaderAttribLocations.luminanceTexture, 0);

	gl.glActiveTexture(GL_TEXTURE1);
	if(!bayerPatternTex_)
		bayerPatternTex_=makeBayerPatternTexture(*sPainter.glFuncs());
	gl.glBindTexture(GL_TEXTURE_2D, bayerPatternTex_);
	luminanceToScreenProgram_->setUniformValue(shaderAttribLocations.bayerPattern, 1);

	gl.glBindVertexArray(vao_);
	gl.glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);
	gl.glBindVertexArray(0);

	luminanceToScreenProgram_->release();
}
