/*
 * Stellarium
 * Copyright (C) 2017 Alexander Wolf
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

#ifndef ASTERISM_HPP
#define ASTERISM_HPP

#include "StelObject.hpp"
#include "StelTranslator.hpp"
#include "StelFader.hpp"
#include "StelSphereGeometry.hpp"
#include "AsterismMgr.hpp"

#include <vector>
#include <QString>
#include <QFont>

class StarMgr;
class StelPainter;
class QJsonObject;

//! @class Asterism
//! The Asterism class models a grouping of stars in a Sky Culture.
//! Each Asterism consists of a list of stars identified by their
//! abbreviation and Hipparcos catalogue numbers (taken from file: asterismship.fab),
//! another entry in file asterism_names.eng.fab with the defining abbreviated name
//! and translatable englishName (translation goes into nameI18).
class Asterism : public StelObject
{
	friend class AsterismMgr;
public:
	static const QString ASTERISM_TYPE;
private:
	Asterism();
	~Asterism() override;

	// StelObject method to override
	//! Get a string with data about the Asterism.
	//! Constellations support the following InfoStringGroup flags:
	//! - Name
	//! @param core the StelCore object
	//! @param flags a set of InfoStringGroup items to include in the return value.
	//! @return a QString a description of the constellation.
	QString getInfoString(const StelCore*, const InfoStringGroup& flags) const override;

	//! Get the module/object type string.
	//! @return "Asterism"
	QString getType() const override {return ASTERISM_TYPE;}
	QString getObjectType() const override {return N_("asterism"); }
	QString getObjectTypeI18n() const override {return q_(getObjectType()); }
	QString getID() const override { return abbreviation; }

	//! observer centered J2000 coordinates.
	//! These are either automatically computed from all stars forming the lines,
	//! or from the manually defined label point(s).
	Vec3d getJ2000EquatorialPos(const StelCore*) const override;

	//! @param record string containing the following whitespace
	//! separated fields: abbreviation - a three character abbreviation
	//! for the asterism, a number of lines (pairs), and a list of Hipparcos
	//! catalogue numbers which, when connected pairwise, form the lines of the
	//! asterism.
	//! @param starMgr a pointer to the StarManager object.
	//! @param excludeRefs a QSet of reference exclusions (user-defined preference in config.ini)
	//! @return false if can't parse record, else true.
	bool read(const QJsonObject& data, StarMgr *starMgr, const QSet<int> &excludeRefs);

	//! Draw the asterism name
	void drawName(const Vec3d &xyName, StelPainter& sPainter) const;

	//! Get the translated name for the Asterism.
	QString getNameI18n() const override {return culturalName.translatedI18n;}
	//! Get the English name for the Asterism.
	QString getEnglishName() const override {return culturalName.translated;}
	//! Get the native name for the Asterism
	QString getNameNative() const override {return culturalName.native;}
	//! Get (translated) pronouncement of the native name for the Asterism
	QString getNamePronounce() const override {return (culturalName.pronounceI18n.isEmpty() ? culturalName.native : culturalName.pronounceI18n);}
	//! Get the short name for the Asterism (returns the translated version of abbreviation).
	QString getShortName() const {return abbreviationI18n;}
	//! Combine screen label from various components, depending on settings in SkyCultureMgr
	QString getScreenLabel() const override;
	//! Combine InfoString label from various components, depending on settings in SkyCultureMgr
	QString getInfoLabel() const override;
	//! Underlying worker
	QString getCultureLabel(StelObject::CulturalDisplayStyle style) const;
	//! Draw the lines for the Asterism.
	//! This method uses the coords of the stars (optimized for use through
	//! the class AsterismMgr only).
	void drawOptim(StelPainter& sPainter, const StelCore* core, const SphericalCap& viewportHalfspace) const;
	//! Update fade levels according to time since various events.
	void update(int deltaTime);
	//! Turn on and off Asterism line rendering.
	//! @param b new state for line drawing.
	void setFlagLines(const bool b) { lineFader=b; }
	//! Turn on and off ray helper rendering.
	//! @param b new state for ray helper drawing.
	void setFlagRayHelpers(const bool b) { rayHelperFader=b; }
	//! Turn on and off Asterism name label rendering.
	//! @param b new state for name label drawing.
	void setFlagLabels(const bool b) { nameFader=b; }
	//! Get the current state of Asterism line rendering.
	//! @return true if Asterism line rendering it turned on, else false.
	bool getFlagLines() const { return lineFader; }
	//! Get the current state of ray helper rendering.
	//! @return true if ray helper rendering it turned on, else false.
	bool getFlagRayHelpers() const { return rayHelperFader; }
	//! Get the current state of Asterism name label rendering.
	//! @return true if Asterism name label rendering it turned on, else false.
	bool getFlagLabels() const { return nameFader; }

	//! @return true if a (real) named asterism, false for a ray helper
	bool isAsterism() const { return flagAsterism; }

	//! Asterism name. This is a culture-dependent thing, and in each skyculture an asterism has one name entry only.
	//! Given multiple aspects of naming, we need all the components and more.
	CulturalName culturalName;
	//! Abbreviation
	//! A skyculture designer must invent it. (usually 2-5 letters)
	//! This MUST be filled and be unique within a sky culture.
	//! @note Given their possible screen use, using numerical labels as abbreviation is not recommended.
	QString abbreviation;
	//! Translated version of abbreviation (the short name or designation of asterism)
	//! Latin-based languages should not translate it, but it may be useful to translate for other glyph systems.
	QString abbreviationI18n;
	//! Context for name
	QString context;
	//! Direction vectors pointing on constellation name drawing position (J2000.0 coordinates)
	//! Usually a single position is computed from averaging star positions forming the constellation, but we can override with an entry in index.json,
	//! and even give more positions (e.g. for long or split-up constellations like Serpens.
	QList<Vec3d> XYZname;
	enum class Type
	{
		RayHelper,          //!< Ray helper
		Asterism,           //!< An asterism with lines between HIP/Gaia stars
		TelescopicAsterism, //!< An asterism with lines defined by J2000.0 coordinates
	};
	//! Type of asterism
	Type typeOfAsterism = Type::Asterism;
	bool flagAsterism; //!< True for genuine asterisms, false for ray helpers
	//! List of stars forming the segments
	std::vector<StelObjectP> asterism;
	//! In case this describes a single-star asterism (i.e. just one line segment that starts and ends at the same star),
	//! or we have a line segment with such single star somewhere within the asterism,
	//! we will draw a circle with this opening radius.
	double singleStarAsterismRadius;

	SphericalCap boundingCap;

	//! Define whether lines and names must be drawn
	LinearFader lineFader, rayHelperFader, nameFader;

	//! Currently we only need one color for all asterisms, this may change at some point
	static Vec3f lineColor;
	static Vec3f rayHelperColor;
	static Vec3f labelColor;
};

#endif // ASTERISM_HPP
