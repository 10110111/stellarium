/*
 * Stellarium
 * Copyright (C) 2003 Fabien Chereau
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

#include "AtmosphereBruneton.hpp"
#include "StelUtils.hpp"
#include "StelApp.hpp"
#include "StelProjector.hpp"
#include "StelToneReproducer.hpp"
#include "StelCore.hpp"
#include "StelPainter.hpp"
#include "StelFileMgr.hpp"
#include "Dithering.hpp"

#include <cstring>

#include <QFile>
#include <QDebug>
#include <QVector>
#include <QSettings>
#include <QOpenGLShaderProgram>

namespace
{

constexpr float radiusOfSun=696000.f; // km
constexpr float radiusOfMoon=1738.f; // km
constexpr float kLengthUnitInMeters = 1000.000000;

template<typename T> T sqr(T x) { return x*x; }

QVector<char> readFileData(QString const& path)
{
	QFile file(path);
	if(!file.open(QIODevice::ReadOnly)) throw std::runtime_error("Failed to open file \""+path.toStdString()+"\"");
	QVector<char> data(file.size());
	if(file.read(data.data(),file.size())!=file.size())
		throw std::runtime_error("Failed to read file \""+path.toStdString()+"\"");
	return data;
}

void bindAndSetupTexture(QOpenGLFunctions& gl, GLenum target, GLuint texture)
{
	gl.glBindTexture(target, texture);
	gl.glTexParameteri(target, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	gl.glTexParameteri(target, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	gl.glTexParameteri(target, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	gl.glTexParameteri(target, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	if(target==GL_TEXTURE_3D)
		gl.glTexParameteri(target, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
}

} // namespace

auto AtmosphereBruneton::getTextureSize4D(QVector<char> const& data) -> TextureSize4D
{
	std::int32_t dim;
	std::memcpy(&dim, &data.back()+1-sizeof dim, sizeof dim);
	if(dim!=4)
	{
		throw std::runtime_error("Bad dimension for 4D texture: "+std::to_string(dim)+
								 " (texture size: "+std::to_string(data.size())+")");
	}

	TextureSize4D out;
	if(unsigned(data.size())<=sizeof dim+sizeof out.sizes)
		throw std::runtime_error("Too small texture file");

	std::memcpy(&out.sizes[0], &data.back()+1-(sizeof dim+sizeof out.sizes), sizeof out.sizes);

	if(unsigned(data.size())!=out.r_size()*out.mu_size()*out.muS_size()*out.nu_size()*sizeof(GLfloat)*4+sizeof dim+sizeof out.sizes)
		throw std::runtime_error("Bad 4D texture size "+std::to_string(data.size()));

	return out;
}

auto AtmosphereBruneton::getTextureSize2D(QVector<char> const& data) -> TextureSize2D
{
	std::int32_t dim;
	std::memcpy(&dim, &data.back()+1-sizeof dim, sizeof dim);
	if(dim!=2)
	{
		throw std::runtime_error("Bad dimension for 2D texture: "+std::to_string(dim)+
								 " (texture size: "+std::to_string(data.size())+")");
	}

	TextureSize2D out;
	if(unsigned(data.size())<=sizeof dim+sizeof out.sizes)
		throw std::runtime_error("Too small texture file");

	std::memcpy(&out.sizes[0], &data.back()+1-(sizeof dim+sizeof out.sizes), sizeof out.sizes);

	if(unsigned(data.size())!=out.width()*out.height()*sizeof(GLfloat)*4+sizeof dim+sizeof out.sizes)
		throw std::runtime_error("Bad 2D texture size "+std::to_string(data.size()));

	return out;
}

void AtmosphereBruneton::loadTextures()
{
	QOpenGLFunctions& gl=*QOpenGLContext::currentContext()->functions();

	// Clear error state to be able to later check for new errors
	while(gl.glGetError()!=GL_NO_ERROR);

	gl.glGenTextures(TEX_COUNT, textures);

	gl.glActiveTexture(GL_TEXTURE0);

	const auto TexImage3D=reinterpret_cast<void(QOPENGLF_APIENTRYP)(GLenum,GLint,GLint,GLsizei,
																	GLsizei,GLsizei,GLint,GLenum,
																	GLenum,const GLvoid*)>
		(QOpenGLContext::currentContext()->getProcAddress(QByteArrayLiteral("glTexImage3D")));
	// We assume at least OpenGL 3.0 (minimal requirements of Stellarium), and there this function must exist.
	Q_ASSERT(TexImage3D);

	{
		bindAndSetupTexture(gl, GL_TEXTURE_2D, textures[TRANSMITTANCE_TEXTURE]);
		const auto data=readFileData(StelFileMgr::findFile("textures/atmosphere/transmittance.dat"));
		transmittanceTextureSize=getTextureSize2D(data);
		gl.glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F,
						transmittanceTextureSize.width(), transmittanceTextureSize.height(),
						0, GL_RGBA, GL_FLOAT, data.data());
	}
	{
		bindAndSetupTexture(gl, GL_TEXTURE_3D, textures[SCATTERING_TEXTURE]);
		const auto data=readFileData(StelFileMgr::findFile("textures/atmosphere/scattering.dat"));
		scatteringTextureSize=getTextureSize4D(data);
		(*TexImage3D)(GL_TEXTURE_3D, 0, GL_RGBA32F,
					  scatteringTextureSize.width(), scatteringTextureSize.height(),
					  scatteringTextureSize.depth(),
					  0, GL_RGBA, GL_FLOAT, data.data());
	}
	{
		bindAndSetupTexture(gl, GL_TEXTURE_2D, textures[IRRADIANCE_TEXTURE]);
		const auto data=readFileData(StelFileMgr::findFile("textures/atmosphere/irradiance.dat"));
		irradianceTextureSize=getTextureSize2D(data);
		gl.glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F,
						irradianceTextureSize.width(), irradianceTextureSize.height(),
						0, GL_RGBA, GL_FLOAT, data.data());
	}
	try
	{
		bindAndSetupTexture(gl, GL_TEXTURE_3D, textures[MIE_SCATTERING_TEXTURE]);
		const auto data=readFileData(StelFileMgr::findFile("textures/atmosphere/mie_scattering.dat"));
		mieScatteringTextureSize=getTextureSize4D(data);
		(*TexImage3D)(GL_TEXTURE_3D, 0, GL_RGBA32F,
					  mieScatteringTextureSize.width(), mieScatteringTextureSize.height(),
					  mieScatteringTextureSize.depth(),
					  0, GL_RGBA, GL_FLOAT, data.data());
		separateMieTexture=true;
	}
	catch(std::runtime_error const&)
	{
		qDebug() << "Separate Mie scattering texture wasn't found, assuming the shader doesn't need it.";
		separateMieTexture=false;
	}
	if(separateMieTexture && mieScatteringTextureSize!=scatteringTextureSize)
		throw std::runtime_error("Mie scattering texture must match Rayleight scattering texture in size");

	gl.glBindTexture(GL_TEXTURE_2D, 0);
	Q_ASSERT(gl.glGetError()==GL_NO_ERROR);
}

void AtmosphereBruneton::loadShaders()
{
	QOpenGLShader vShader(QOpenGLShader::Vertex);
	if (!vShader.compileSourceFile(":/shaders/AtmosphereMain.vert"))
	{
		qFatal("Error while compiling atmosphere vertex shader: %s", vShader.log().toLatin1().constData());
	}
	if (!vShader.log().isEmpty())
	{
		qWarning() << "Warnings while compiling atmosphere vertex shader: " << vShader.log();
	}
	QOpenGLShader ditherShader(QOpenGLShader::Fragment);
	if (!ditherShader.compileSourceCode(makeDitheringShader()))
	{
		qFatal("Error while compiling atmosphere dithering shader: %s", ditherShader.log().toLatin1().constData());
	}
	if (!ditherShader.log().isEmpty())
	{
		qWarning() << "Warnings while compiling atmosphere dithering shader: " << ditherShader.log();
	}
	QOpenGLShader fShader(QOpenGLShader::Fragment);
	if (!fShader.compileSourceFile(":/shaders/AtmosphereMain.frag"))
	{
		qFatal("Error while compiling atmosphere fragment shader: %s", fShader.log().toLatin1().constData());
	}
	if (!fShader.log().isEmpty())
	{
		qWarning() << "Warnings while compiling atmosphere fragment shader: " << fShader.log();
	}
	QOpenGLShader toneReproducerShader(QOpenGLShader::Fragment);
	if (!toneReproducerShader.compileSourceFile(":/shaders/xyYToRGB.frag"))
	{
		qFatal("Error while compiling atmosphere tone reproducer shader: %s", toneReproducerShader.log().toLatin1().constData());
	}
	if (!toneReproducerShader.log().isEmpty())
	{
		qWarning() << "Warnings while compiling atmosphere tone reproducer shader: " << toneReproducerShader.log();
	}
	QOpenGLShader funcShader(QOpenGLShader::Fragment);
	{
		QFile file(":/shaders/AtmosphereFunctions.frag");
		if(!file.open(QIODevice::ReadOnly))
			throw std::runtime_error("Failed to read file "+file.fileName().toStdString());
		QString source=file.readAll();
		source.replace(QRegExp("(TRANSMITTANCE_TEXTURE_WIDTH = )[^;]+;"),  QString("\\1%1;").arg(transmittanceTextureSize.width()));
		source.replace(QRegExp("(TRANSMITTANCE_TEXTURE_HEIGHT = )[^;]+;"), QString("\\1%1;").arg(transmittanceTextureSize.height()));
		source.replace(QRegExp("(IRRADIANCE_TEXTURE_WIDTH = )[^;]+;"),     QString("\\1%1;").arg(irradianceTextureSize.width()));
		source.replace(QRegExp("(IRRADIANCE_TEXTURE_HEIGHT = )[^;]+;"),    QString("\\1%1;").arg(irradianceTextureSize.height()));
		source.replace(QRegExp("(SCATTERING_TEXTURE_R_SIZE = )[^;]+;"),    QString("\\1%1;").arg(scatteringTextureSize.r_size()));
		source.replace(QRegExp("(SCATTERING_TEXTURE_MU_SIZE = )[^;]+;"),   QString("\\1%1;").arg(scatteringTextureSize.mu_size()));
		source.replace(QRegExp("(SCATTERING_TEXTURE_MU_S_SIZE = )[^;]+;"), QString("\\1%1;").arg(scatteringTextureSize.muS_size()));
		source.replace(QRegExp("(SCATTERING_TEXTURE_NU_SIZE = )[^;]+;"),   QString("\\1%1;").arg(scatteringTextureSize.nu_size()));
		if(separateMieTexture) source.replace("#define COMBINED_SCATTERING_TEXTURES","");
		if (!funcShader.compileSourceCode(source))
		{
			qFatal("Error while compiling atmosphere render functions shader: %s", funcShader.log().toLatin1().constData());
		}
		if (!funcShader.log().isEmpty())
		{
			qWarning() << "Warnings while compiling atmosphere render functions shader: " << funcShader.log();
		}
	}
	atmoShaderProgram->addShader(&vShader);
	atmoShaderProgram->addShader(&ditherShader);
	atmoShaderProgram->addShader(&fShader);
	atmoShaderProgram->addShader(&toneReproducerShader);
	atmoShaderProgram->addShader(&funcShader);
	StelPainter::linkProg(atmoShaderProgram.get(), "atmosphere");
}

AtmosphereBruneton::AtmosphereBruneton()
	: viewport(0,0,0,0)
	, gridMaxY(44)
	, gridMaxX(44)
	, posGridBuffer(QOpenGLBuffer::VertexBuffer)
	, indexBuffer(QOpenGLBuffer::IndexBuffer)
	, viewRayGridBuffer(QOpenGLBuffer::VertexBuffer)
	, averageLuminance(0)
	, overrideAverageLuminance(false)
	, eclipseFactor(1)
	, lightPollutionLuminance(0)
	, atmoShaderProgram(new QOpenGLShaderProgram())
{
	setFadeDuration(1.5);

	loadTextures();
	loadShaders();

	atmoShaderProgram->bind();
	shaderAttribLocations.bayerPattern = atmoShaderProgram->uniformLocation("bayerPattern");
	shaderAttribLocations.rgbMaxValue = atmoShaderProgram->uniformLocation("rgbMaxValue");
	shaderAttribLocations.alphaWaOverAlphaDa = atmoShaderProgram->uniformLocation("alphaWaOverAlphaDa");
	shaderAttribLocations.oneOverGamma = atmoShaderProgram->uniformLocation("oneOverGamma");
	shaderAttribLocations.term2TimesOneOverMaxdLpOneOverGamma = atmoShaderProgram->uniformLocation("term2TimesOneOverMaxdLpOneOverGamma");
	shaderAttribLocations.brightnessScale = atmoShaderProgram->uniformLocation("brightnessScale");
	shaderAttribLocations.sunDir = atmoShaderProgram->uniformLocation("sun_direction");
	shaderAttribLocations.cameraPos = atmoShaderProgram->uniformLocation("camera");
	shaderAttribLocations.projectionMatrix = atmoShaderProgram->uniformLocation("projectionMatrix");
	shaderAttribLocations.skyVertex = atmoShaderProgram->attributeLocation("vertex");
	shaderAttribLocations.viewRay = atmoShaderProgram->attributeLocation("viewRay");
	shaderAttribLocations.transmittanceTexture = atmoShaderProgram->uniformLocation("transmittance_texture");
	shaderAttribLocations.scatteringTexture = atmoShaderProgram->uniformLocation("scattering_texture");
	shaderAttribLocations.irradianceTexture = atmoShaderProgram->uniformLocation("irradiance_texture");
	shaderAttribLocations.singleMieScatteringTexture = atmoShaderProgram->uniformLocation("single_mie_scattering_texture");
	atmoShaderProgram->release();
}

AtmosphereBruneton::~AtmosphereBruneton()
{
	if(auto*const ctx=QOpenGLContext::currentContext())
		ctx->functions()->glDeleteTextures(TEX_COUNT, textures);
}

void AtmosphereBruneton::regenerateGrid()
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

void AtmosphereBruneton::updateEclipseFactor(StelCore* core, Vec3d sunPos, Vec3d moonPos)
{
	// Update the eclipse intensity factor to apply on atmosphere model
	// these are for radii
	const float sunAngularSize = std::atan(radiusOfSun/AU/sunPos.length());
	const float moonAngularSize = std::atan(radiusOfMoon/AU/moonPos.length());
	const float touchAngle = sunAngularSize + moonAngularSize;

	// determine luminance falloff during solar eclipses
	sunPos.normalize();
	moonPos.normalize();
	const auto separationAngle = std::acos(sunPos.dot(moonPos));  // angle between them
	// qDebug("touch at %f\tnow at %f (%f)\n", touchAngle, separationAngle, separationAngle/touchAngle);
	// bright stars should be visible at total eclipse
	// TODO: correct for atmospheric diffusion
	// TODO: use better coverage function (non-linear)
	// because of above issues, this algorithm darkens more quickly than reality
	// Note: On Earth only, else moon would brighten other planets' atmospheres (LP:1673283)
	if ((core->getCurrentLocation().planetName=="Earth") && (separationAngle < touchAngle))
	{
		float darkAngle = moonAngularSize - sunAngularSize;
		float min = 0.0025f; // 0.005f; // 0.0001f;  // so bright stars show up at total eclipse
		if (darkAngle < 0)
		{
			// annular eclipse
			const auto asun = sqr(sunAngularSize);
			min = (asun - sqr(moonAngularSize))/asun;  // minimum proportion of sun uncovered
			darkAngle *= -1;
		}

		if (separationAngle < darkAngle)
			eclipseFactor = min;
		else
			eclipseFactor = min + (1-min)*(separationAngle-darkAngle)/(touchAngle-darkAngle);
	}
	else
		eclipseFactor = 1;
	// TODO: compute eclipse factor also for Lunar eclipses! (lp:#1471546)
}

void AtmosphereBruneton::computeColor(double JD, Vec3d sunPos, Vec3d moonPos, float moonPhase, float moonMagnitude,
							   StelCore* core, float latitude, float altitude, float temperature, float relativeHumidity)
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

	if (std::isnan(sunPos.length()))
		sunPos.set(0.,0.,-1.*AU);
	if (std::isnan(moonPos.length()))
		moonPos.set(0.,0.,-1.*AU);

	updateEclipseFactor(core, sunPos, moonPos);

	// No need to calculate if not visible
	if (!fader.getInterstate())
	{
		// GZ 20180114: Why did we add light pollution if atmosphere was not visible?????
		// And what is the meaning of 0.001? Approximate contribution of stellar background? Then why is it 0.0001 below???
		averageLuminance = 0.001f;
		return;
	}

	sunPos.normalize();
	moonPos.normalize();
	sunDir=Vec3d(sunPos[0],sunPos[1],sunPos[2]);
	this->altitude=altitude;

	// Calculate the atmosphere RGB for each point of the grid

	skyb.setLocation(latitude * M_PI/180, altitude, temperature, relativeHumidity);
	skyb.setSunMoon(moonPos[2], sunPos[2]);

	// Calculate the date from the julian day.
	int year, month, day;
	StelUtils::getDateFromJulianDay(JD, &year, &month, &day);
	skyb.setDate(year, month, moonPhase, moonMagnitude);

	// Variables used to compute the average sky luminance
	float sumLum = 0;

	Vec3d point(1, 0, 0);
	float lumi;

	// Compute the sky color for every point above the ground
	for (int i=0; i<(1+gridMaxX)*(1+gridMaxY); ++i)
	{
		const Vec2f &v(posGrid[i]);
		prj->unProject(v[0],v[1],point);

		Q_ASSERT(fabs(point.lengthSquared()-1.0) < 1e-10);

		lumi = skyb.getLuminance(moonPos[0]*point[0]+moonPos[1]*point[1]+
				 moonPos[2]*point[2], sunPos[0]*point[0]+sunPos[1]*point[1]+
				 sunPos[2]*point[2], point[2]);
		lumi *= eclipseFactor;
		// Add star background luminance
		lumi += 0.0001f;
		// Multiply by the input scale of the ToneConverter (is not done automatically by the xyYtoRGB method called later)
		//lumi*=eye->getInputScale();

		// Add the light pollution luminance AFTER the scaling to avoid scaling it because it is the cause
		// of the scaling itself
		lumi += fader.getInterstate()*lightPollutionLuminance;

		// Store for later statistics
		sumLum+=lumi;

		// Now need to compute the xy part of the color component
		// This is done in the openGL shader
		// Store the back projected position + luminance in the input color to the shader
		viewRayGrid[i].set(point[0], point[1], point[2], lumi);
	}

	viewRayGridBuffer.bind();
	viewRayGridBuffer.write(0, &viewRayGrid[0], viewRayGrid.size()*sizeof viewRayGrid[0]);
	viewRayGridBuffer.release();

	// Update average luminance
	if (!overrideAverageLuminance)
		averageLuminance = sumLum/((1+gridMaxX)*(1+gridMaxY));
}

// override computable luminance. This is for special operations only, e.g. for scripting of brightness-balanced image export.
// To return to auto-computed values, set any negative value.
void AtmosphereBruneton::setAverageLuminance(float overrideLum)
{
	if (overrideLum<0.f)
	{
		overrideAverageLuminance=false;
		averageLuminance=0.f;
	}
	else
	{
		overrideAverageLuminance=true;
		averageLuminance=overrideLum;
	}
}

// Draw the atmosphere using the precalc values stored in tab_sky
void AtmosphereBruneton::draw(StelCore* core)
{
	if (StelApp::getInstance().getVisionModeNight())
		return;

	StelToneReproducer* eye = core->getToneReproducer();

	if (!fader.getInterstate())
		return;

	StelPainter sPainter(core->getProjection2d());
	sPainter.setBlending(true, GL_ONE, GL_ONE);

	const float atm_intensity = fader.getInterstate();

	atmoShaderProgram->bind();
	float a, b, c;
	eye->getShadersParams(a, b, c);
	atmoShaderProgram->setUniformValue(shaderAttribLocations.alphaWaOverAlphaDa, a);
	atmoShaderProgram->setUniformValue(shaderAttribLocations.oneOverGamma, b);
	atmoShaderProgram->setUniformValue(shaderAttribLocations.term2TimesOneOverMaxdLpOneOverGamma, c);
	atmoShaderProgram->setUniformValue(shaderAttribLocations.brightnessScale, atm_intensity);
	atmoShaderProgram->setUniformValue(shaderAttribLocations.sunDir, sunDir[0], sunDir[1], sunDir[2]);
	atmoShaderProgram->setUniformValue(shaderAttribLocations.cameraPos, 0, 0, altitude/kLengthUnitInMeters);
	const Mat4f& m = sPainter.getProjector()->getProjectionMatrix();
	atmoShaderProgram->setUniformValue(shaderAttribLocations.projectionMatrix,
		QMatrix4x4(m[0], m[4], m[8], m[12], m[1], m[5], m[9], m[13], m[2], m[6], m[10], m[14], m[3], m[7], m[11], m[15]));

	const auto rgbMaxValue=calcRGBMaxValue(sPainter.getDitheringMode());
	atmoShaderProgram->setUniformValue(shaderAttribLocations.rgbMaxValue, rgbMaxValue[0], rgbMaxValue[1], rgbMaxValue[2]);
	auto& gl=*sPainter.glFuncs();

	gl.glActiveTexture(GL_TEXTURE0);
	gl.glBindTexture(GL_TEXTURE_2D, textures[TRANSMITTANCE_TEXTURE]);
	gl.glUniform1i(shaderAttribLocations.transmittanceTexture, 0);
	gl.glActiveTexture(GL_TEXTURE1);
	gl.glBindTexture(GL_TEXTURE_3D, textures[SCATTERING_TEXTURE]);
	gl.glUniform1i(shaderAttribLocations.scatteringTexture, 1);
	gl.glActiveTexture(GL_TEXTURE2);
	gl.glBindTexture(GL_TEXTURE_2D, textures[IRRADIANCE_TEXTURE]);
	gl.glUniform1i(shaderAttribLocations.irradianceTexture, 2);
	gl.glActiveTexture(GL_TEXTURE3);
	gl.glBindTexture(GL_TEXTURE_3D, textures[MIE_SCATTERING_TEXTURE]);
	gl.glUniform1i(shaderAttribLocations.singleMieScatteringTexture, 3);

	gl.glActiveTexture(GL_TEXTURE4);
	if(!bayerPatternTex)
		bayerPatternTex=makeBayerPatternTexture(*sPainter.glFuncs());
	gl.glBindTexture(GL_TEXTURE_2D, bayerPatternTex);
	atmoShaderProgram->setUniformValue(shaderAttribLocations.bayerPattern, 4);

	viewRayGridBuffer.bind();
	atmoShaderProgram->setAttributeBuffer(shaderAttribLocations.viewRay, GL_FLOAT, 0, 4, 0);
	viewRayGridBuffer.release();
	atmoShaderProgram->enableAttributeArray(shaderAttribLocations.viewRay);
	posGridBuffer.bind();
	atmoShaderProgram->setAttributeBuffer(shaderAttribLocations.skyVertex, GL_FLOAT, 0, 2, 0);
	posGridBuffer.release();
	atmoShaderProgram->enableAttributeArray(shaderAttribLocations.skyVertex);

	// And draw everything at once
	indexBuffer.bind();
	std::size_t shift=0;
	for (int y=0;y<gridMaxY;++y)
	{
		sPainter.glFuncs()->glDrawElements(GL_TRIANGLE_STRIP, (gridMaxX+1)*2, GL_UNSIGNED_SHORT, reinterpret_cast<void*>(shift));
		shift += (gridMaxX+1)*2*2;
	}
	indexBuffer.release();

	atmoShaderProgram->disableAttributeArray(shaderAttribLocations.skyVertex);
	atmoShaderProgram->disableAttributeArray(shaderAttribLocations.viewRay);
	atmoShaderProgram->release();
}
