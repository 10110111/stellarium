/*
 * Copyright (C) 2020 Georg Zotti
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

#ifndef MAYAHAABCALENDAR_HPP
#define MAYAHAABCALENDAR_HPP

#include "Calendar.hpp"

//! The Maya Haab was a 365-day Solar calendar without intercalation.
//! Similar to the Egyptian calendar, after 18 months of 20 days there was a short "month" of 5 extra days.
//! The implementation follows CC.
class MayaHaabCalendar : public Calendar
{
	Q_OBJECT

public:
	MayaHaabCalendar(double jd);

	~MayaHaabCalendar() override {}

public slots:
	void retranslate() override;

	//! Set a calendar date from the Julian day number
	void setJD(double JD) override;

	//! set date from a vector of calendar date elements sorted from the largest to the smallest.
	//! month[1..19]-day[0..19]
	//! We face a problem as the year is not counted. We can only find the date before current JD which matches the parts.
	void setDate(const QVector<int> &parts) override;

	//! get a stringlist of calendar date elements sorted from the largest to the smallest.
	//! monthName-day[0..19]
	QStringList getDateStrings() const override;

	//! get a formatted complete string for a date
	QString getFormattedDateString() const override;

	//! get tzolkin name index of Haab year bearer (name of 0 Pop) from Haab date
	//! This must be one of 2, 7, 12, 17. (TODO: write a test!)
	static int mayanYearBearerFromFixed(int rd);

	//! get RD of a given calendar round date on or before rd. They repeat every 18980 days.
	static int mayanCalendarRoundOnOrBefore(const QVector<int> &haab, const QVector<int> &tzolkin, int rd);

	inline static int mayanHaabOrdinal(const QVector<int> &haab) {return (haab.at(0)-1)*20+haab.at(1);}
	static int mayanHaabOnOrBefore(const QVector<int> &haab, int rd);
	static QVector<int> mayanHaabFromFixed(int rd);
private:
	static QMap<int, QString> monthNames;
	static const int mayanHaabEpoch;
};

#endif
