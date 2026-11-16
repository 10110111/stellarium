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
 
#ifndef SUNSPOTSDIALOG_HPP
#define SUNSPOTSDIALOG_HPP

#include <memory>
#include <QTimer>
#include "Sunspots.hpp"
#include "StelDialog.hpp"

class Ui_SunspotsDialog;

class SunspotsDialog : public StelDialog
{
	Q_OBJECT

public:
	SunspotsDialog(Sunspots* sunspots);
	~SunspotsDialog() override;

protected:
	void createDialogContent() override;

public slots:
	void retranslate() override;

private:
	void setAboutText();
	void updateStatus();
	void updateImageTimes();
	void updateCacheList();
	void filterCacheList();

	Sunspots* sunspots_ = nullptr;
	std::unique_ptr<Ui_SunspotsDialog> ui_;
	QTimer presentationTimer_;
	QTimer timeChangeTimer_;
};

#endif // SUNSPOTSDIALOG_HPP
