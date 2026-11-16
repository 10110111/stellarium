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

#ifndef SUNSPOTSIMAGE_HPP
#define SUNSPOTSIMAGE_HPP

#include <QImage>
#include <QMap>
#include <QString>
#include <functional>

namespace SunspotsImage
{
	//! Stages reported while creating the cached solar map.
	enum class PreparationStage { Preparing, Saving };

	struct Result
	{
		QString error;
		bool cancelled = false;
		QMap<QString, QString> metadata;
	};

	//! Read image/WCS keywords from a FITS header or the FITS XML embedded in a PHI JP2.
	//! Throws std::runtime_error on source errors. PHI sources contain their own metadata.
	QMap<QString, QString> readMetadata(const QString& sourcePath);

	//! Expected quiet-Sun continuum intensity at 617.3 nm, normalized at disk center.
	//! mu is the cosine of the angle between the surface normal and direction to the spacecraft.
	double phiLimbDarkening(double mu);

	//! Prepare an equirectangular, opaque indexed-color PNG from an HMI/PHI FITS image or PHI JPEG2000 quicklook.
	//! Read source data, undo JPEG2000 intensity encoding when needed, reproject it, normalize colors,
	//! construct the indexed palette, and save the PNG.
	//! The sRGB palette normalizes observed intensities to the mean Sun color and fills missing samples with it.
	//! Exclude a source-image margin at the geometric solar limb before interpolating.
	//! The viewing direction comes from the FITS Carrington observer coordinates.
	//! Preparation uses 0--90% for geometry and 90--100% for color preparation and palette remapping.
	//! reportProgress reports the stage and percentage.
	//! Saving reports 50 because PNG writing has no progress callback.
	Result prepareTexture(const QString& sourcePath, const QString& pngPath,
	                      const QMap<QString, QString>& metadata,
	                      const std::function<void(PreparationStage, int)>& reportProgress,
	                      const std::function<bool()>& isCancelled = {});

	//! Decode a finished cached PNG, preserving its palette for the texture loader to expand.
	//! Cancellation is checked before and after Qt's non-interruptible decoding operation.
	QImage loadTexture(const QString& pngPath, QString* error,
	                   const std::function<bool()>& isCancelled = {});
}

#endif
