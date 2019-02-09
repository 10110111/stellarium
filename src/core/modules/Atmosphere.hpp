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

#ifndef ATMOSPHERE_HPP
#define ATMOSPHERE_HPP

#include "StelCore.hpp"
#include "VecMath.hpp"

class Atmosphere
{
public:
	virtual ~Atmosphere() = default;
	//! Compute sky brightness values and average luminance.
	virtual void computeColor(double JD, Vec3d _sunPos, Vec3d moonPos, float moonPhase, float moonMagnitude, StelCore* core,
	                          float latitude = 45.f, float altitude = 200.f, float temperature = 15.f, float relativeHumidity = 40.f) = 0;
	virtual void draw(StelCore* core) = 0;
	virtual void update(double deltaTime) = 0;

	//! Set fade in/out duration in seconds
	virtual void setFadeDuration(float duration) = 0;
	//! Get fade in/out duration in seconds
	virtual float getFadeDuration() const = 0;

	//! Define whether to display atmosphere
	virtual void setFlagShow(bool b) = 0;
	//! Get whether atmosphere is displayed
	virtual bool getFlagShow() const = 0;

	//! Get the actual atmosphere intensity due to eclipses + fader
	//! @return the display intensity ranging from 0 to 1
	virtual float getRealDisplayIntensityFactor() const = 0;

	// lets you know how far faded in or out the atmosphere is (0..1)
	virtual float getFadeIntensity() const = 0;

	/*!
     * Get the average luminance of the atmosphere in cd/m2
	 * If atmosphere is off, the luminance equals the background starlight (0.001cd/m2).
	 * Otherwise it includes the (atmosphere + background starlight (0.0001cd/m2) * eclipse factor + light pollution.
	 * @return the last computed average luminance of the atmosphere in cd/m2.
	 */
	virtual float getAverageLuminance() const = 0;
	//! override computable luminance. This is for special operations only, e.g. for scripting of brightness-balanced image export.
	//! To return to auto-computed values, set any negative value at the end of the script.
	virtual void setAverageLuminance(float overrideLum) = 0;
	//! Set the light pollution luminance in cd/m^2
	virtual void setLightPollutionLuminance(float f) = 0;
	//! Get the light pollution luminance in cd/m^2
	virtual float getLightPollutionLuminance() const = 0;
};

#endif
