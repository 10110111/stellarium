/*
 * Copyright (C) 2014 Alexander Wolf
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

#include "StelOpenGL.hpp"
#include <QDebug>
#include "StelMainView.hpp"
#include <QOpenGLFunctions_3_3_Core>
#if QT_VERSION >= QT_VERSION_CHECK(6,0,0)
# include <QOpenGLVersionFunctionsFactory>
#endif

QOpenGLContext* StelOpenGL::mainContext = Q_NULLPTR;

const char* StelOpenGL::getGLErrorText(GLenum code) {
	switch (code) {
		case GL_INVALID_ENUM:
			return "GL_INVALID_ENUM";
		case GL_INVALID_FRAMEBUFFER_OPERATION:
			return "GL_INVALID_FRAMEBUFFER_OPERATION";
		case GL_INVALID_VALUE:
			return "GL_INVALID_VALUE";
		case GL_INVALID_OPERATION:
			return "GL_INVALID_OPERATION";
		case GL_OUT_OF_MEMORY:
			return "GL_OUT_OF_MEMORY";
		default:
			return "unknown GL error type";
	}
}

int StelOpenGL::checkGLErrors(const char *file, int line)
{
	int errors = 0;
	for(GLenum x; (x=mainContext->functions()->glGetError())!=GL_NO_ERROR; )
	{
		++errors;
		qCritical("%s:%d: OpenGL error: %d (%s)\n",
			  file, line, x, getGLErrorText(x));
	}
	return errors;
}

void StelOpenGL::clearGLErrors()
{
	while(mainContext->functions()->glGetError()!=GL_NO_ERROR)
	{ }
}

QByteArray StelOpenGL::smoothTexSamplingFunction()
{
	static const QByteArray src = 1+R"(
#line 1 109
uniform int textureSmoothingLevel;
vec4 smoothSampleTexture(in sampler2D tex, vec2 texc)
{
	if(textureSmoothingLevel <= 0)
	{
		return texture2D(tex, texc);
	}
	else if(textureSmoothingLevel == 1)
	{
		vec2 delta = 0.35/vec2(textureSize(tex,0));
		vec4 col0  = texture2D(tex, texc);
		vec4 col1a = texture2D(tex, texc +   vec2( delta.s,  delta.t));
		vec4 col1b = texture2D(tex, texc +   vec2(-delta.s,  delta.t));
		vec4 col1c = texture2D(tex, texc +   vec2( delta.s, -delta.t));
		vec4 col1d = texture2D(tex, texc +   vec2(-delta.s, -delta.t));
		return 1./3.*(col0+0.5*(col1a+col1b+col1c+col1d));
	}
	else
	{
		vec2 delta = 0.35/vec2(textureSize(tex,0));
		vec4 col0  = texture2D(tex, texc);
		vec4 col1a = texture2D(tex, texc +   vec2( delta.s,  delta.t));
		vec4 col1b = texture2D(tex, texc +   vec2(-delta.s,  delta.t));
		vec4 col1c = texture2D(tex, texc +   vec2( delta.s, -delta.t));
		vec4 col1d = texture2D(tex, texc +   vec2(-delta.s, -delta.t));
		vec4 col2a = texture2D(tex, texc + 2.*vec2( delta.s, 0));
		vec4 col2b = texture2D(tex, texc + 2.*vec2(       0,  delta.t));
		vec4 col2c = texture2D(tex, texc + 2.*vec2(-delta.s, 0));
		vec4 col2d = texture2D(tex, texc + 2.*vec2(       0, -delta.t));
		return 0.2*(col0 + 0.7*(col1a+col1b+col1c+col1d) + 0.3*(col2a+col2b+col2c+col2d));
	}
}
#line 1 0
)";
	return src;
}

QByteArray StelOpenGL::globalShaderPrefix(const ShaderType type)
{
	const auto& glInfo = StelMainView::getInstance().getGLInformation();

	if(glInfo.isHighGraphicsMode)
	{
		if(type == VERTEX_SHADER)
		{
			static const QByteArray prefix = 1+R"(
#version 330
#define ATTRIBUTE in
#define VARYING out
#define texture2D(a,b) texture(a,b)
#define texture2D_3(a,b,c) texture(a,b,c)
#define texture2DProj(a,b) textureProj(a,b)
#define texture2DProj_3(a,b,c) textureProj(a,b,c)
#define shadow2D_x(a,b) texture(a,b)
#define shadow2DProj_x(a,b) textureProj(a,b)
#define textureCube(a,b) texture(a,b)
#define textureCube_3(a,b,c) texture(a,b,c)
#line 1
)";
			return prefix;
		}
		else // FRAGMENT_SHADER
		{
			static const QByteArray prefix = 1+R"(
#version 330
#define VARYING in
out vec4 FRAG_COLOR;
#define texture2D(a,b) texture(a,b)
#define texture2D_3(a,b,c) texture(a,b,c)
#define texture2DProj(a,b) textureProj(a,b)
#define texture2DProj_3(a,b,c) textureProj(a,b,c)
#define shadow2D_x(a,b) texture(a,b)
#define shadow2DProj_x(a,b) textureProj(a,b)
#define textureCube(a,b) texture(a,b)
#define textureCube_3(a,b,c) texture(a,b,c)
#define textureGrad_SUPPORTED
#line 1
)";
			return prefix;
		}
	}

	if(type == VERTEX_SHADER)
	{
		static const QByteArray prefix = 1+R"(
#define ATTRIBUTE attribute
#define VARYING varying
#define texture2D_3(a,b,c) texture2D(a,b,c)
#define texture2DProj_3(a,b,c) texture2DProj(a,b,c)
#define shadow2D_x(a,b) shadow2D(a,b).x
#define shadow2DProj_x(a,b) shadow2DProj(a,b).x
#define textureCube_3(a,b,c) textureCube(a,b,c)
#line 1
)";
		return prefix;
	}
	else // FRAGMENT_SHADER
	{
		static const QByteArray prefix = 1+R"(
#define VARYING varying
#define FRAG_COLOR gl_FragColor
#define texture2D_3(a,b,c) texture2D(a,b,c)
#define texture2DProj_3(a,b,c) texture2DProj(a,b,c)
#define shadow2D_x(a,b) shadow2D(a,b).x
#define shadow2DProj_x(a,b) shadow2DProj(a,b).x
#define textureCube_3(a,b,c) textureCube(a,b,c)
#line 1
)";
		if(glInfo.isGLES)
			return "precision mediump float;\n" + prefix;
		return prefix;
	}
}

QOpenGLFunctions_3_3_Core* StelOpenGL::highGraphicsFunctions()
{
#if !QT_CONFIG(opengles2)
# if QT_VERSION >= QT_VERSION_CHECK(6,0,0)
	return QOpenGLVersionFunctionsFactory::get<QOpenGLFunctions_3_3_Core>(QOpenGLContext::currentContext());
# else
	return QOpenGLContext::currentContext()->versionFunctions<QOpenGLFunctions_3_3_Core>();
# endif
#else
	return nullptr;
#endif
}
