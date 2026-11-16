/*
 * Stellarium: Sunspots plugin
 * Copyright (C) 2026 Ruslan Kabatsayev
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

#include "SunspotsImage.hpp"
#include "StelSRGB.hpp"
#include "StelTranslator.hpp"

#include <fitsio.h>
#include <QColor>
#include <QFile>
#include <QSaveFile>
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace
{
constexpr int textureWidth = 8192;
constexpr int textureHeight = 4096;
constexpr double radPerDegree = M_PI / 180.0;
constexpr double arcsecPerRadian = 180.0 * 3600.0 / M_PI;
// Mean color for both observed and unseen photosphere.
const QColor meanSunColor(0xf2, 0xe6, 0xe2);

struct FitsImage
{
	long width = 0;
	long height = 0;
	std::vector<float> pixels;
};

using IntensityHistogram = std::array<quint32, 256>;

QString fitsError(int status)
{
	char message[FLEN_STATUS];
	fits_get_errstatus(status, message);
	return QString::fromLatin1(message);
}

QString readKeyword(fitsfile* file, const char* name)
{
	char value[FLEN_VALUE];
	int status = 0;
	if (fits_read_key(file, TSTRING, name, value, nullptr, &status))
		return {};
	return QString::fromLatin1(value).trimmed();
}

double number(const QMap<QString, QString>& metadata, const QString& key)
{
	const auto it = metadata.constFind(key);
	if (it == metadata.cend())
		throw std::runtime_error(q_("FITS image is missing the numeric keyword %1.").arg(key).toStdString());
	bool ok = false;
	const double value = it.value().toDouble(&ok);
	if (!ok || !std::isfinite(value))
		throw std::runtime_error(q_("FITS image has a malformed numeric value for %1: %2.")
		                         .arg(key, it.value()).toStdString());
	return value;
}

double sample(const std::vector<float>& pixels, long width, long height, double x, double y)
{
	if (x < 0. || y < 0. || x >= width - 1 || y >= height - 1)
		return NAN;
	const size_t ix = x, iy = y;
	const double dx = x - ix, dy = y - iy;
	const double a = pixels[iy * width + ix];
	const double b = pixels[iy * width + ix + 1];
	const double c = pixels[(iy + 1) * width + ix];
	const double d = pixels[(iy + 1) * width + ix + 1];
	return (1. - dy) * ((1. - dx) * a + dx * b) + dy * ((1. - dx) * c + dx * d);
}

bool readFits(const QString& fitsPath, FitsImage& fits, SunspotsImage::Result& result)
{
	fitsfile* file = nullptr;
	int status = 0;
	const QByteArray name = QFile::encodeName(fitsPath);
	if (fits_open_file(&file, name.constData(), READONLY, &status))
	{
		result.error = q_("Cannot open FITS image: %1").arg(fitsError(status));
		return false;
	}

	int hduCount = 0, imageDim = 0, hduType = 0;
	fits_get_num_hdus(file, &hduCount, &status);
	long dimensions[2] = {0, 0};
	for (int hdu = 1; hdu <= hduCount && !status; ++hdu)
	{
		fits_movabs_hdu(file, hdu, &hduType, &status);
		fits_get_img_dim(file, &imageDim, &status);
		if (imageDim == 2)
		{
			fits_get_img_size(file, 2, dimensions, &status);
			if (dimensions[0] > 0 && dimensions[1] > 0)
				break;
		}
	}
	if (status || dimensions[0] < 2 || dimensions[1] < 2 ||
	    dimensions[0] > 16384 || dimensions[1] > 16384)
	{
		result.error = q_("FITS file has no supported two-dimensional image: %1").arg(fitsError(status));
		fits_close_file(file, &status);
		return false;
	}

	const char* keys[] = {"T_OBS", "T_REC", "CRLT_OBS", "CRLN_OBS", "CROTA2", "CRPIX1", "CRPIX2",
	                      "CRVAL1", "CRVAL2", "CDELT1", "CDELT2", "CTYPE1", "CTYPE2", "CUNIT1", "CUNIT2",
	                      "DSUN_OBS", "RSUN_REF", "RSUN_OBS", "QUALITY"};
	for (const char* key : keys)
	{
		const QString value = readKeyword(file, key);
		if (!value.isEmpty()) result.metadata.insert(QString::fromLatin1(key), value);
	}

	fits.width = dimensions[0];
	fits.height = dimensions[1];
	fits.pixels.resize(size_t(fits.width) * fits.height);
	auto& pixels = fits.pixels;
	int anyNull = 0;
	float nullValue = std::numeric_limits<float>::quiet_NaN();
	status = 0;
	fits_read_img(file, TFLOAT, 1, static_cast<LONGLONG>(pixels.size()), &nullValue,
	              pixels.data(), &anyNull, &status);
	int closeStatus = 0;
	fits_close_file(file, &closeStatus);
	if (status)
	{
		result.error = q_("Cannot read FITS pixels: %1").arg(fitsError(status));
		return false;
	}
	return true;
}

QImage reproject(FitsImage& fits, const Mat4d& modelToEcliptic, const Vec3d& observerEcliptic,
	IntensityHistogram& histogram, SunspotsImage::Result& result,
	const std::function<void(int)>& reportProgress, const std::function<bool()>& isCancelled)
try
{
	auto& pixels = fits.pixels;
	const float nullValue = std::numeric_limits<float>::quiet_NaN();
	const auto& meta = result.metadata;
	if (meta.value("CTYPE1", "HPLN-TAN") != "HPLN-TAN" ||
	    meta.value("CTYPE2", "HPLT-TAN") != "HPLT-TAN")
	{
		result.error = q_("Unsupported FITS projection (expected HPLN-TAN/HPLT-TAN).");
		return {};
	}
	const double crpix1 = number(meta, "CRPIX1");
	const double crpix2 = number(meta, "CRPIX2");
	const double unit1 = meta.value("CUNIT1", "arcsec").contains("deg", Qt::CaseInsensitive) ? 3600. : 1.;
	const double unit2 = meta.value("CUNIT2", "arcsec").contains("deg", Qt::CaseInsensitive) ? 3600. : 1.;
	const double cdelt1 = number(meta, "CDELT1") * unit1;
	const double cdelt2 = number(meta, "CDELT2") * unit2;
	const double angle = number(meta, "CROTA2") * radPerDegree;
	const double c = cos(angle), s = sin(angle);
	// HMI supplies pixel scales in CDELT1/2 and image rotation in CROTA2.
	const double a = cdelt1 * c;
	const double b = -cdelt1 * s;
	const double d = cdelt2 * s;
	const double e = cdelt2 * c;
	const double determinant = a * e - b * d;
	if (!std::isfinite(determinant) || fabs(determinant) < 1.e-12)
	{
		result.error = q_("FITS image has invalid pixel scale or rotation.");
		return {};
	}

	const double apparentRadius = number(meta, "RSUN_OBS");
	const double distance = 1. / sin(apparentRadius / arcsecPerRadian);
	if (apparentRadius <= 0. || !std::isfinite(distance) || distance <= 1.)
	{
		result.error = q_("FITS image has invalid observer distance.");
		return {};
	}
	const Vec3d view = normalize(modelToEcliptic.inverse() * observerEcliptic);
	const Vec3d north = normalize(Vec3d(0., 0., 1.) - view * view[2]);
	const Vec3d east = normalize(north ^ view);
	const double crval1 = number(meta, "CRVAL1") * unit1;
	const double crval2 = number(meta, "CRVAL2") * unit2;
	if (!std::isfinite(crval1) || !std::isfinite(crval2))
	{
		result.error = q_("FITS image has invalid reference coordinates.");
		return {};
	}
	const double centerX = crpix1 - 1. - (e * crval1 - b * crval2) / determinant;
	const double centerY = crpix2 - 1. - (-d * crval1 + a * crval2) / determinant;
	// Radius of the apparent limb in the TAN projection plane, in squared arcseconds.
	const double limbRadiusSquared = arcsecPerRadian * arcsecPerRadian / (distance * distance - 1.);
	// Width of the excluded annulus inside the predicted solar limb, measured radially in FITS pixels.
	// Sampling the very edge produced a rim that was too dark in the north and too bright in the south.
	// Treat these pixels as missing before interpolation so affected texels receive the mean Sun color.
	constexpr double limbMarginPixels = 2.;

	size_t validCount = 0;
	double maximum = 0.;
	for (long y = 0; y < fits.height; ++y)
	{
		if (isCancelled && isCancelled())
		{
			result.cancelled = true;
			return {};
		}
		for (long x = 0; x < fits.width; ++x)
		{
			float& intensity = pixels[size_t(y) * fits.width + x];
			const double dx = x - centerX, dy = y - centerY;
			const double radius = std::hypot(dx, dy);
			// Test the margin farther out along the radius in FITS pixels. Masking source pixels also
			// rejects any bilinear footprint touching the excluded rim. Keep the WCS mapping unchanged.
			const double marginScale = radius > 0. ? 1. + limbMarginPixels / radius : 1.;
			const double worldX = (a * dx + b * dy) * marginScale;
			const double worldY = (d * dx + e * dy) * marginScale;
			if (!std::isfinite(intensity) || intensity < 0. ||
			    worldX * worldX + worldY * worldY >= limbRadiusSquared)
			{
				intensity = nullValue;
				continue;
			}
			++validCount;
			maximum = qMax(maximum, static_cast<double>(intensity));
		}
	}
	if (validCount < 16)
	{
		result.error = q_("FITS image contains too few valid pixels.");
		return {};
	}
	if (maximum <= 0.)
	{
		result.error = q_("FITS image has no usable intensity range.");
		return {};
	}

	QImage image(textureWidth, textureHeight, QImage::Format_Indexed8);
	if (image.isNull())
	{
		result.error = q_("Cannot allocate solar texture.");
		return {};
	}
	image.fill(0);
	for (int y = 0; y < textureHeight; ++y)
	{
		if (isCancelled && isCancelled())
		{
			result.cancelled = true;
			return {};
		}
		const double latitude = M_PI * (0.5 - (y + 0.5) / textureHeight);
		const double sinLat = sin(latitude), cosLat = cos(latitude);
		uchar* row = image.scanLine(y);
		for (int x = 0; x < textureWidth; ++x)
		{
			const double longitude = 2. * M_PI * (x + 0.5) / textureWidth;
			const Vec3d point(-sin(longitude) * cosLat, cos(longitude) * cosLat, sinLat);
			const double towardObserver = point * view;
			if (towardObserver <= 1. / distance) continue;
			const double denominator = distance - towardObserver;
			const double worldX = (point * east) / denominator * arcsecPerRadian - crval1;
			const double worldY = (point * north) / denominator * arcsecPerRadian - crval2;
			const double imageX = (e * worldX - b * worldY) / determinant + crpix1 - 1.;
			const double imageY = (-d * worldX + a * worldY) / determinant + crpix2 - 1.;
			const double intensity = sample(pixels, fits.width, fits.height, imageX, imageY);
			if (!std::isfinite(intensity) || intensity < 0.) continue;
			// Preserve zero exclusively for missing samples, including at the quantization floor.
			row[x] = qBound(1, qRound(255. * linearToSRGB(intensity / maximum)), 255);
			++histogram[row[x]];
		}
		if (y % 16 == 0) reportProgress(100 * y / textureHeight);
	}
	reportProgress(100);
	return image;
}
catch (const std::runtime_error& error)
{
	result.error = QString::fromUtf8(error.what());
	return {};
}

// Input indices represent sRGB intensities, with zero reserved for missing samples.
// Replace them with indices into a finished, opaque palette sorted by brightness.
bool colorizeAndIndex(QImage& image, const IntensityHistogram& histogram, SunspotsImage::Result& result,
	const std::function<void(int)>& reportProgress,
	const std::function<bool()>& isCancelled)
{
	std::array<double, 256> linearIntensities{};
	double sum = 0.;
	quint64 count = 0;
	// Exclude missing samples and calculate the mean in linear light, before tinting.
	for (size_t gray = 1; gray < histogram.size(); ++gray)
	{
		linearIntensities[gray] = srgbToLinear(gray / 255.);
		sum += histogram[gray] * linearIntensities[gray];
		count += histogram[gray];
	}
	if (count == 0 || sum <= 0.)
	{
		result.error = q_("Projected image has no usable intensity values.");
		return false;
	}
	const double mean = sum / count;
	const Vec3d linearColor(srgbToLinear(meanSunColor.redF()),
	                        srgbToLinear(meanSunColor.greenF()),
	                        srgbToLinear(meanSunColor.blueF()));
	const auto meanAtGain = [&](double gain)
	{
		double total = 0.;
		for (size_t gray = 1; gray < histogram.size(); ++gray)
			if (histogram[gray] != 0)
				total += histogram[gray] * qMin(1., gain * linearIntensities[gray] / mean);
		return total / count;
	};
	std::array<double, 3> gains;
	for (int channel = 0; channel < 3; ++channel)
	{
		const double target = linearColor[channel];
		if (target <= 0.)
		{
			gains[channel] = 0.;
			continue;
		}
		double lower = 0., upper = qMax(1., 2. * target);
		while (meanAtGain(upper) < target && upper < 1.e12) upper *= 2.;
		for (int iteration = 0; iteration < 32; ++iteration)
		{
			if (isCancelled && isCancelled())
			{
				result.cancelled = true;
				return false;
			}
			const double midpoint = (lower + upper) / 2.;
			if (meanAtGain(midpoint) < target) lower = midpoint;
			else upper = midpoint;
		}
		gains[channel] = upper;
	}
	std::array<QRgb, 256> colors;
	colors[0] = qRgb(meanSunColor.red(), meanSunColor.green(), meanSunColor.blue());
	for (size_t gray = 1; gray < colors.size(); ++gray)
	{
		const double intensity = linearIntensities[gray] / mean;
		const auto channel = [intensity](double gain)
		{
			return qBound(0., std::round(255. * linearToSRGB(qMin(1., gain * intensity))), 255.);
		};
		colors[gray] = qRgb(channel(gains[0]), channel(gains[1]), channel(gains[2]));
	}
	std::array<double, 256> luminances;
	std::array<uchar, 256> order;
	for (size_t index = 0; index < colors.size(); ++index)
	{
		const QRgb color = colors[index];
		luminances[index] = 0.2126 * srgbToLinear(qRed(color) / 255.) +
		                    0.7152 * srgbToLinear(qGreen(color) / 255.) +
		                    0.0722 * srgbToLinear(qBlue(color) / 255.);
		order[index] = index;
	}
	std::stable_sort(order.begin(), order.end(), [&luminances](uchar a, uchar b)
	{
		return luminances[a] < luminances[b];
	});
	QVector<QRgb> palette(256);
	std::array<uchar, 256> remapping;
	for (int index = 0; index < palette.size(); ++index)
	{
		palette[index] = colors[order[index]];
		remapping[order[index]] = index;
	}
	image.setColorTable(palette);
	for (int y = 0; y < image.height(); ++y)
	{
		if (isCancelled && isCancelled())
		{
			result.cancelled = true;
			return false;
		}
		uchar* row = image.scanLine(y);
		for (int x = 0; x < image.width(); ++x)
			row[x] = remapping[row[x]];
		if (y % 16 == 0) reportProgress(100 * y / image.height());
	}
	reportProgress(100);
	return true;
}

void savePng(const QImage& image, const QString& pngPath, SunspotsImage::Result& result,
	const std::function<bool()>& isCancelled)
{
	QSaveFile output(pngPath);
	if (!output.open(QIODevice::WriteOnly) || !image.save(&output, "PNG"))
		result.error = q_("Cannot save projected PNG: %1").arg(output.errorString());
	else if (isCancelled && isCancelled())
	{
		output.cancelWriting();
		result.cancelled = true;
	}
	else if (!output.commit())
		result.error = q_("Cannot save projected PNG: %1").arg(output.errorString());
}
}

SunspotsImage::Result SunspotsImage::prepareTexture(const QString& fitsPath, const QString& pngPath,
	const QMap<QString, QString>& jsocMetadata, const Mat4d& modelToEcliptic,
	const Vec3d& observerEcliptic, const std::function<void(PreparationStage, int)>& reportProgress,
	const std::function<bool()>& isCancelled)
{
	Result result;
	result.metadata = jsocMetadata;
	if (isCancelled && isCancelled())
	{
		result.cancelled = true;
		return result;
	}

	FitsImage fits;
	if (!readFits(fitsPath, fits, result)) return result;
	if (isCancelled && isCancelled())
	{
		result.cancelled = true;
		return result;
	}

	IntensityHistogram histogram{};
	constexpr int reprojectionProgressPercent = 90;
	QImage image = reproject(fits, modelToEcliptic, observerEcliptic, histogram, result,
		[&reportProgress](int value)
		{
			reportProgress(PreparationStage::Preparing, reprojectionProgressPercent * value / 100);
		},
		isCancelled);
	if (image.isNull()) return result;
	if (!colorizeAndIndex(image, histogram, result,
		[&reportProgress](int value)
		{
			reportProgress(PreparationStage::Preparing,
			               reprojectionProgressPercent + (100 - reprojectionProgressPercent) * value / 100);
		},
		isCancelled)) return result;
	if (isCancelled && isCancelled())
	{
		result.cancelled = true;
		return result;
	}

	reportProgress(PreparationStage::Saving, 50);
	savePng(image, pngPath, result, isCancelled);
	return result;
}

QImage SunspotsImage::loadTexture(const QString& pngPath, QString* error,
	const std::function<bool()>& isCancelled)
{
	if (isCancelled && isCancelled()) return {};
	QImage image(pngPath);
	if (isCancelled && isCancelled()) return {};
	if (image.isNull())
	{
		*error = q_("Cannot read cached PNG: %1").arg(pngPath);
		return {};
	}
	return image;
}
