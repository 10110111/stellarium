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
constexpr float nonExtinctedSolarMagnitude=-26.74;
constexpr float kLengthUnitInMeters = 1000.000000;

template<typename T> T sqr(T x) { return x*x; }

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

void checkFramebufferStatus(QOpenGLFunctions& gl)
{
	GLenum status=gl.glCheckFramebufferStatus(GL_FRAMEBUFFER);
	if(status!=GL_FRAMEBUFFER_COMPLETE)
	{
		QString errorDescription;
		switch(status)
		{
		case GL_FRAMEBUFFER_INCOMPLETE_ATTACHMENT:
			errorDescription="incomplete attachment";
			break;
		case GL_FRAMEBUFFER_INCOMPLETE_MISSING_ATTACHMENT:
			errorDescription="missing attachment";
			break;
		case GL_INVALID_FRAMEBUFFER_OPERATION:
			errorDescription="invalid framebuffer operation";
			break;
		case GL_FRAMEBUFFER_UNSUPPORTED:
			errorDescription="framebuffer unsupported";
			break;
		default:
			errorDescription=QString("Unknown error %1").arg(status);
			break;
		}
		qWarning() << "Error: framebuffer incomplete:" << errorDescription.toStdString().c_str();
	}
}

} // namespace

QVector<char> AtmosphereBruneton::readFileData(QString const& path)
{
	QFile file(path);
	if(!file.open(QIODevice::ReadOnly)) throw InitFailure("Failed to open file \""+path.toStdString()+"\"");
	QVector<char> data(file.size());
	if(file.read(data.data(),file.size())!=file.size())
		throw InitFailure("Failed to read file \""+path.toStdString()+"\"");
	return data;
}

auto AtmosphereBruneton::getTextureSize4D(QVector<char> const& data) -> TextureSize4D
{
	std::int32_t dim;
	std::memcpy(&dim, &data.back()+1-sizeof dim, sizeof dim);
	if(dim!=4)
	{
		throw InitFailure("Bad dimension for 4D texture: "+std::to_string(dim)+
								 " (texture size: "+std::to_string(data.size())+")");
	}

	TextureSize4D out;
	if(unsigned(data.size())<=sizeof dim+sizeof out.sizes)
		throw InitFailure("Too small texture file");

	std::memcpy(&out.sizes[0], &data.back()+1-(sizeof dim+sizeof out.sizes), sizeof out.sizes);

	if(unsigned(data.size())!=out.r_size()*out.mu_size()*out.muS_size()*out.nu_size()*sizeof(GLfloat)*4+sizeof dim+sizeof out.sizes)
		throw InitFailure("Bad 4D texture size "+std::to_string(data.size()));

	return out;
}

auto AtmosphereBruneton::getTextureSize2D(QVector<char> const& data) -> TextureSize2D
{
	std::int32_t dim;
	std::memcpy(&dim, &data.back()+1-sizeof dim, sizeof dim);
	if(dim!=2)
	{
		throw InitFailure("Bad dimension for 2D texture: "+std::to_string(dim)+
								 " (texture size: "+std::to_string(data.size())+")");
	}

	TextureSize2D out;
	if(unsigned(data.size())<=sizeof dim+sizeof out.sizes)
		throw InitFailure("Too small texture file");

	std::memcpy(&out.sizes[0], &data.back()+1-(sizeof dim+sizeof out.sizes), sizeof out.sizes);

	if(unsigned(data.size())!=out.width()*out.height()*sizeof(GLfloat)*4+sizeof dim+sizeof out.sizes)
		throw InitFailure("Bad 2D texture size "+std::to_string(data.size()));

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
	// We require at least OpenGL 3.0 (minimal requirements of Stellarium), and there this function must exist.
	if(!TexImage3D)
		throw InitFailure("glTexImage3D function not found");

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
		throw InitFailure("Mie scattering texture must match Rayleigh scattering texture in size");

	gl.glBindTexture(GL_TEXTURE_2D, 0);
	const auto error=gl.glGetError();
	if(error!=GL_NO_ERROR)
		throw InitFailure("Texture loading: glGetError returned "+std::to_string(error));
}

void AtmosphereBruneton::loadShaders()
{
	// Core atmosphere rendering shader program
	{
		QOpenGLShader vShader(QOpenGLShader::Vertex);
		if (!vShader.compileSourceFile(":/shaders/AtmosphereMain.vert"))
		{
			qCritical("Error while compiling atmosphere vertex shader: %s", vShader.log().toLatin1().constData());
		}
		if (!vShader.log().isEmpty())
		{
			qWarning() << "Warnings while compiling atmosphere vertex shader: " << vShader.log();
		}
		QOpenGLShader fShader(QOpenGLShader::Fragment);
		if (!fShader.compileSourceFile(":/shaders/AtmosphereMain.frag"))
		{
			qCritical("Error while compiling atmosphere fragment shader: %s", fShader.log().toLatin1().constData());
			throw InitFailure("Shader compilation failed");
		}
		if (!fShader.log().isEmpty())
		{
			qWarning() << "Warnings while compiling atmosphere fragment shader: " << fShader.log();
		}
		QOpenGLShader funcShader(QOpenGLShader::Fragment);
		{
			QFile file(":/shaders/AtmosphereFunctions.frag");
			if(!file.open(QIODevice::ReadOnly))
				throw InitFailure("Failed to read file "+file.fileName().toStdString());
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
				qCritical("Error while compiling atmosphere render functions shader: %s", funcShader.log().toLatin1().constData());
				throw InitFailure("Shader compilation failed");
			}
			if (!funcShader.log().isEmpty())
			{
				qWarning() << "Warnings while compiling atmosphere render functions shader: " << funcShader.log();
			}
		}
		atmosphereRenderProgram->addShader(&vShader);
		atmosphereRenderProgram->addShader(&fShader);
		atmosphereRenderProgram->addShader(&funcShader);
		StelPainter::linkProg(atmosphereRenderProgram.get(), "atmosphere");
	}

	// Post-processing shader program
	{
		QOpenGLShader vShader(QOpenGLShader::Vertex);
		constexpr char src[]=R"(
#version 130
in vec4 vertex;
out vec2 texCoord;
void main()
{
    gl_Position=vertex;
    texCoord=vertex.xy*0.5+0.5;
}
)";
		if (!vShader.compileSourceCode(src))
		{
			qCritical("Error while compiling atmosphere post-processing vertex shader: %s", vShader.log().toLatin1().constData());
		}
		if (!vShader.log().isEmpty())
		{
			qWarning() << "Warnings while compiling atmosphere post-processing vertex shader: " << vShader.log();
		}
		QOpenGLShader ditherShader(QOpenGLShader::Fragment);
		if (!ditherShader.compileSourceCode(makeDitheringShader()))
		{
			qCritical("Error while compiling atmosphere dithering shader: %s", ditherShader.log().toLatin1().constData());
			throw InitFailure("Shader compilation failed");
		}
		if (!ditherShader.log().isEmpty())
		{
			qWarning() << "Warnings while compiling atmosphere dithering shader: " << ditherShader.log();
		}
		QOpenGLShader toneReproducerShader(QOpenGLShader::Fragment);
		if (!toneReproducerShader.compileSourceFile(":/shaders/AtmospherePostProcess.frag"))
		{
			qCritical("Error while compiling atmosphere tone reproducer shader: %s", toneReproducerShader.log().toLatin1().constData());
			throw InitFailure("Shader compilation failed");
		}
		if (!toneReproducerShader.log().isEmpty())
		{
			qWarning() << "Warnings while compiling atmosphere tone reproducer shader: " << toneReproducerShader.log();
		}
		postProcessProgram->addShader(&vShader);
		postProcessProgram->addShader(&ditherShader);
		postProcessProgram->addShader(&toneReproducerShader);
		if(!StelPainter::linkProg(postProcessProgram.get(), "atmosphere-post-process"))
			throw InitFailure("Shader program linking failed");
	}
}

void AtmosphereBruneton::resizeRenderTarget(int width, int height)
{
	Q_ASSERT(fbo);
	QOpenGLFunctions& gl=*QOpenGLContext::currentContext()->functions();
	gl.glBindTexture(GL_TEXTURE_2D,textures[FBO_TEXTURE]);
	gl.glTexImage2D(GL_TEXTURE_2D,0,GL_RGBA32F,width,height,0,GL_RGBA,GL_UNSIGNED_BYTE,nullptr);
	gl.glBindTexture(GL_TEXTURE_2D,0);
	gl.glBindFramebuffer(GL_FRAMEBUFFER,fbo);
	gl.glFramebufferTexture2D(GL_FRAMEBUFFER,GL_COLOR_ATTACHMENT0,GL_TEXTURE_2D,textures[FBO_TEXTURE],0);
	checkFramebufferStatus(gl);
	gl.glBindFramebuffer(GL_FRAMEBUFFER,0);

	fboPrevWidth=width;
	fboPrevHeight=height;
}

void AtmosphereBruneton::setupRenderTarget()
{
	QOpenGLFunctions& gl=*QOpenGLContext::currentContext()->functions();
	gl.glGenFramebuffers(1,&fbo);
	gl.glBindTexture(GL_TEXTURE_2D,textures[FBO_TEXTURE]);
	gl.glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_NEAREST);
	gl.glBindTexture(GL_TEXTURE_2D,0);

	GLint viewport[4];
	gl.glGetIntegerv(GL_VIEWPORT, viewport);
	resizeRenderTarget(viewport[2], viewport[3]);
}

void AtmosphereBruneton::setupBuffers()
{
	QOpenGLFunctions& gl=*QOpenGLContext::currentContext()->functions();

	const auto GenVertexArrays=reinterpret_cast<void(QOPENGLF_APIENTRYP)(GLsizei,GLuint*)>
		(QOpenGLContext::currentContext()->getProcAddress(QByteArrayLiteral("glGenVertexArrays")));
	const auto BindVertexArray=reinterpret_cast<void(QOPENGLF_APIENTRYP)(GLuint)>
		(QOpenGLContext::currentContext()->getProcAddress(QByteArrayLiteral("glBindVertexArray")));
	// We require at least OpenGL 3.0 (minimal requirements of Stellarium), and there these functions must exist.
	if(!GenVertexArrays) throw InitFailure("glGenVertexArrays function not found");
	if(!BindVertexArray) throw InitFailure("glBindVertexArray function not found");

	(*GenVertexArrays)(1, &vao);
	(*BindVertexArray)(vao);
	gl.glGenBuffers(1, &vbo);
	gl.glBindBuffer(GL_ARRAY_BUFFER, vbo);
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
	(*BindVertexArray)(0);
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
	, atmosphereRenderProgram(new QOpenGLShaderProgram())
	, postProcessProgram(new QOpenGLShaderProgram())
{
	setFadeDuration(1.5);

	loadTextures();
	loadShaders();
	setupRenderTarget();
	setupBuffers();

	atmosphereRenderProgram->bind();
	 shaderAttribLocations.cameraPos                  = atmosphereRenderProgram->uniformLocation("camera");
	 shaderAttribLocations.sunDir                     = atmosphereRenderProgram->uniformLocation("sun_direction");
	 shaderAttribLocations.projectionMatrix           = atmosphereRenderProgram->uniformLocation("projectionMatrix");
	 shaderAttribLocations.transmittanceTexture       = atmosphereRenderProgram->uniformLocation("transmittance_texture");
	 shaderAttribLocations.scatteringTexture          = atmosphereRenderProgram->uniformLocation("scattering_texture");
	 shaderAttribLocations.irradianceTexture          = atmosphereRenderProgram->uniformLocation("irradiance_texture");
	 shaderAttribLocations.singleMieScatteringTexture = atmosphereRenderProgram->uniformLocation("single_mie_scattering_texture");
	 shaderAttribLocations.skyVertex                  = atmosphereRenderProgram->attributeLocation("vertex");
	 shaderAttribLocations.viewRay                    = atmosphereRenderProgram->attributeLocation("viewRay");
	atmosphereRenderProgram->release();

	postProcessProgram->bind();
	 shaderAttribLocations.rgbMaxValue                         = postProcessProgram->uniformLocation("rgbMaxValue");
	 shaderAttribLocations.bayerPattern                        = postProcessProgram->uniformLocation("bayerPattern");
	 shaderAttribLocations.radianceImage                       = postProcessProgram->uniformLocation("radiance");
	 shaderAttribLocations.alphaWaOverAlphaDa                  = postProcessProgram->uniformLocation("alphaWaOverAlphaDa");
	 shaderAttribLocations.oneOverGamma                        = postProcessProgram->uniformLocation("oneOverGamma");
	 shaderAttribLocations.term2TimesOneOverMaxdLpOneOverGamma = postProcessProgram->uniformLocation("term2TimesOneOverMaxdLpOneOverGamma");
	 shaderAttribLocations.brightnessScale                     = postProcessProgram->uniformLocation("brightnessScale");
	postProcessProgram->release();
}

AtmosphereBruneton::~AtmosphereBruneton()
{
	if(auto*const ctx=QOpenGLContext::currentContext())
	{
		auto& gl=*ctx->functions();
		gl.glDeleteTextures(TEX_COUNT, textures);
		gl.glDeleteBuffers(1, &vbo);
		const auto DeleteVertexArrays=reinterpret_cast<void(QOPENGLF_APIENTRYP)(GLsizei,const GLuint*)>
			(QOpenGLContext::currentContext()->getProcAddress(QByteArrayLiteral("glDeleteVertexArrays")));
		(*DeleteVertexArrays)(1, &vao);
	}
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

void AtmosphereBruneton::drawAtmosphere(Mat4f const& projectionMatrix, Vec3d sunDir, float brightness)
{
	atmosphereRenderProgram->bind();
	atmosphereRenderProgram->setUniformValue(shaderAttribLocations.sunDir, sunDir[0], sunDir[1], sunDir[2]);
	atmosphereRenderProgram->setUniformValue(shaderAttribLocations.cameraPos, 0, 0, altitude/kLengthUnitInMeters);
	const auto& m = projectionMatrix;
	atmosphereRenderProgram->setUniformValue(shaderAttribLocations.projectionMatrix,
											 QMatrix4x4(m[0], m[4], m[8], m[12],
														m[1], m[5], m[9], m[13],
														m[2], m[6], m[10], m[14],
														m[3], m[7], m[11], m[15]));

	viewRayGridBuffer.bind();
	atmosphereRenderProgram->setAttributeBuffer(shaderAttribLocations.viewRay, GL_FLOAT, 0, 4, 0);
	viewRayGridBuffer.release();
	atmosphereRenderProgram->enableAttributeArray(shaderAttribLocations.viewRay);
	posGridBuffer.bind();
	atmosphereRenderProgram->setAttributeBuffer(shaderAttribLocations.skyVertex, GL_FLOAT, 0, 2, 0);
	posGridBuffer.release();
	atmosphereRenderProgram->enableAttributeArray(shaderAttribLocations.skyVertex);

	QOpenGLFunctions& gl=*QOpenGLContext::currentContext()->functions();

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

	gl.glBlendFunc(GL_CONSTANT_COLOR, GL_ONE);
	gl.glBlendColor(brightness, brightness, brightness, brightness);
	gl.glEnable(GL_BLEND);

	indexBuffer.bind();
	std::size_t shift=0;
	for (int y=0;y<gridMaxY;++y)
	{
		gl.glDrawElements(GL_TRIANGLE_STRIP, (gridMaxX+1)*2, GL_UNSIGNED_SHORT, reinterpret_cast<void*>(shift));
		shift += (gridMaxX+1)*2*2;
	}
	indexBuffer.release();

	// FIXME: maybe restore blending state instead of setting to defaults?
	gl.glBlendFunc(GL_ONE, GL_ZERO);
	gl.glBlendColor(0,0,0,0);
	gl.glDisable(GL_BLEND);

	atmosphereRenderProgram->disableAttributeArray(shaderAttribLocations.skyVertex);
	atmosphereRenderProgram->disableAttributeArray(shaderAttribLocations.viewRay);
	atmosphereRenderProgram->release();
}

Vec4f AtmosphereBruneton::getMeanPixelValue(int texW, int texH)
{
	QOpenGLFunctions& gl=*QOpenGLContext::currentContext()->functions();

	gl.glActiveTexture(GL_TEXTURE0);
	gl.glBindTexture(GL_TEXTURE_2D, textures[FBO_TEXTURE]);
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
	const auto GetTexImage=reinterpret_cast<void(QOPENGLF_APIENTRYP)(GLenum,GLint,GLenum,GLenum,GLvoid*)>
							(QOpenGLContext::currentContext()->getProcAddress(QByteArrayLiteral("glGetTexImage")));
	// We require at least OpenGL 3.0 (minimal requirements of Stellarium), and there this function must exist.
	if(!GetTexImage)
		throw InitFailure("glGetTexImage function not found");
	(*GetTexImage)(GL_TEXTURE_2D, deepestLevel, GL_RGBA, GL_FLOAT, &pixel[0]);
	return pixel;
}

void AtmosphereBruneton::computeColor(double JD, Vec3d sunPos, Vec3d moonPos, float moonPhase, float moonMagnitude, float nonExtinctedLunarMagnitude,
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
	const auto width=viewport[2], height=viewport[3];
	if(width!=fboPrevWidth || height!=fboPrevHeight)
		resizeRenderTarget(width, height);

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
	const auto sunDir=sunPos;
	const auto moonDir=moonPos;
	this->altitude=altitude;

	// Calculate the atmosphere RGB for each point of the grid

	skyb.setLocation(latitude * M_PI/180, altitude, temperature, relativeHumidity);
	skyb.setSunMoon(moonPos[2], sunPos[2]);

	// Calculate the date from the julian day.
	int year, month, day;
	StelUtils::getDateFromJulianDay(JD, &year, &month, &day);
	skyb.setDate(year, month, moonPhase, moonMagnitude);

	// TODO:
	// 1. (DONE) Luminance due to the Sun and the Moon
	// 2. 		 Add airglow from Skybright (put Sun & Moon to nadir to avoid their influence)
	// 3. 		 Take eclipsed Sun into account.
	// 4. (DONE) Take eclipsed Moon into account.
	// 5. (DONE) Add star background luminance
	// 6.(P.DONE)Add light pollution luminance (done only to calculate average luminance)
	// 7.		 Draw light pollution

	for (int i=0; i<(1+gridMaxX)*(1+gridMaxY); ++i)
	{
		Vec3d point(1, 0, 0);
		prj->unProject(posGrid[i][0],posGrid[i][1],point);
		viewRayGrid[i].set(point[0], point[1], point[2], 1);
	}

	viewRayGridBuffer.bind();
	viewRayGridBuffer.write(0, &viewRayGrid[0], viewRayGrid.size()*sizeof viewRayGrid[0]);
	viewRayGridBuffer.release();

	QOpenGLFunctions& gl=*QOpenGLContext::currentContext()->functions();
	GLint prevFBO;
	gl.glGetIntegerv(GL_FRAMEBUFFER_BINDING, &prevFBO);
	gl.glBindFramebuffer(GL_FRAMEBUFFER, fbo);
	gl.glClearColor(0,0,0,0);
	gl.glClear(GL_COLOR_BUFFER_BIT);
	 drawAtmosphere(prj->getProjectionMatrix(), sunDir, 1.0);
	 const auto moonRelativeBrightness=std::pow(10.f,0.4f*(nonExtinctedSolarMagnitude-nonExtinctedLunarMagnitude));
	 drawAtmosphere(prj->getProjectionMatrix(), moonDir, moonRelativeBrightness);
	gl.glBindFramebuffer(GL_FRAMEBUFFER, prevFBO);

	if (!overrideAverageLuminance)
	{
		const auto meanPixelValue=getMeanPixelValue(width, height);
		// CIE 1931 luminance computed from linear sRGB
		const auto meanY=0.2126729*meanPixelValue[0]+0.7151522*meanPixelValue[1]+0.0721750*meanPixelValue[2];

		const auto starBackgroundLuminance=1e-4f;
		averageLuminance = meanY + starBackgroundLuminance + fader.getInterstate()*lightPollutionLuminance;
	}
}

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

void AtmosphereBruneton::draw(StelCore* core)
{
	if (StelApp::getInstance().getVisionModeNight())
		return;

	StelToneReproducer* eye = core->getToneReproducer();

	if (!fader.getInterstate())
		return;

	StelPainter sPainter(core->getProjection2d());
	sPainter.setBlending(true, GL_ONE, GL_ONE);
	auto& gl=*sPainter.glFuncs();

	const float atm_intensity = fader.getInterstate();

	postProcessProgram->bind();
	float a, b, c;
	eye->getShadersParams(a, b, c);
	postProcessProgram->setUniformValue(shaderAttribLocations.alphaWaOverAlphaDa, a);
	postProcessProgram->setUniformValue(shaderAttribLocations.oneOverGamma, b);
	postProcessProgram->setUniformValue(shaderAttribLocations.term2TimesOneOverMaxdLpOneOverGamma, c);
	postProcessProgram->setUniformValue(shaderAttribLocations.brightnessScale, atm_intensity);

	const auto rgbMaxValue=calcRGBMaxValue(sPainter.getDitheringMode());
	postProcessProgram->setUniformValue(shaderAttribLocations.rgbMaxValue, rgbMaxValue[0], rgbMaxValue[1], rgbMaxValue[2]);

	gl.glActiveTexture(GL_TEXTURE0);
	gl.glBindTexture(GL_TEXTURE_2D, textures[FBO_TEXTURE]);
	postProcessProgram->setUniformValue(shaderAttribLocations.radianceImage, 0);

	gl.glActiveTexture(GL_TEXTURE1);
	if(!bayerPatternTex)
		bayerPatternTex=makeBayerPatternTexture(*sPainter.glFuncs());
	gl.glBindTexture(GL_TEXTURE_2D, bayerPatternTex);
	postProcessProgram->setUniformValue(shaderAttribLocations.bayerPattern, 1);

	const auto BindVertexArray=reinterpret_cast<void(QOPENGLF_APIENTRYP)(GLuint)>
								(QOpenGLContext::currentContext()->getProcAddress(QByteArrayLiteral("glBindVertexArray")));
	(*BindVertexArray)(vao);
	gl.glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);
	(*BindVertexArray)(0);

	postProcessProgram->release();
}
