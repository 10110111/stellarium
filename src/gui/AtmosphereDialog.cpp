/*
 * Stellarium
 * Copyright (C) 2011 Georg Zotti
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
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Suite 500, Boston, MA  02110-1335, USA.
*/

#include "StelApp.hpp"
#include "StelCore.hpp"
#include "StelSkyDrawer.hpp"
#include "StelTranslator.hpp"
#include "AtmosphereDialog.hpp"
#include "ui_atmosphereDialog.h"

#include <QFileInfo>
#include <QFileDialog>

namespace
{
constexpr char MODEL_PROPERTY[]      = "LandscapeMgr.atmosphereModel";
constexpr char MODEL_PATH_PROPERTY[] = "LandscapeMgr.atmosphereModelPath";
}

AtmosphereDialog::AtmosphereDialog()
	: StelDialog("Atmosphere")
	, refraction(Q_NULLPTR)
	, extinction(Q_NULLPTR)

{
	ui = new Ui_atmosphereDialogForm;
}

AtmosphereDialog::~AtmosphereDialog()
{
	delete ui;
}

void AtmosphereDialog::retranslate()
{
	if (dialog)
		ui->retranslateUi(dialog);
}


void AtmosphereDialog::createDialogContent()
{
	ui->setupUi(dialog);
	
	//Signals and slots
	connect(&StelApp::getInstance(), SIGNAL(languageChanged()), this, SLOT(retranslate()));
	connect(ui->closeStelWindow, SIGNAL(clicked()), this, SLOT(close()));
	connect(ui->TitleBar, SIGNAL(movedTo(QPoint)), this, SLOT(handleMovedTo(QPoint)));

	connectDoubleProperty(ui->pressureDoubleSpinBox,"StelSkyDrawer.atmospherePressure");
	connectDoubleProperty(ui->temperatureDoubleSpinBox,"StelSkyDrawer.atmosphereTemperature");
	connectDoubleProperty(ui->extinctionDoubleSpinBox,"StelSkyDrawer.extinctionCoefficient");

	connect(ui->atmosphereModel, &QComboBox::currentTextChanged, this, &AtmosphereDialog::onModelChoiceChanged);
	connect(ui->showMySky_pathToModelBrowseBtn, &QPushButton::clicked, this, &AtmosphereDialog::browsePathToModel);
	connect(ui->showMySky_pathToModelEdit, &QLineEdit::textChanged, this, &AtmosphereDialog::onPathToModelChanged);
	connect(ui->showMySky_debugOptionsEnabled, &QCheckBox::toggled, this,
			[this](const bool enabled){ ui->showMySky_debugOptionsGroup->setVisible(enabled); });
	connectBoolProperty(ui->showMySky_zeroOrderEnabled, "LandscapeMgr.flagAtmosphereZeroOrderScattering");
	connectBoolProperty(ui->showMySky_singleScatteringEnabled, "LandscapeMgr.flagAtmosphereSingleScattering");
	connectBoolProperty(ui->showMySky_multipleScatteringEnabled, "LandscapeMgr.flagAtmosphereMultipleScattering");

	setCurrentValues();
}

void AtmosphereDialog::onModelChoiceChanged(const QString& model)
{
	const bool isShowMySky = model.toLower()=="showmysky";
	ui->showMySky_optionsGroup->setVisible(isShowMySky);
	if(!isShowMySky || this->hasValidModelPath())
		StelApp::getInstance().getStelPropertyManager()->setStelPropertyValue(MODEL_PROPERTY, model);

}

void AtmosphereDialog::browsePathToModel()
{
	const auto path=QFileDialog::getExistingDirectory(nullptr, q_("Open ShowMySky model"));
	if(path.isNull()) return;
	ui->showMySky_pathToModelEdit->setText(path);
	StelApp::getInstance().getStelPropertyManager()->setStelPropertyValue(MODEL_PATH_PROPERTY, path);
}

bool AtmosphereDialog::hasValidModelPath() const
{
	const auto text=ui->showMySky_pathToModelEdit->text();
	if(text.isEmpty()) return false;
	return QFileInfo(text+"/params.atmo").exists();
}

void AtmosphereDialog::onPathToModelChanged()
{
	if(hasValidModelPath())
		ui->showMySky_pathToModelEdit->setStyleSheet("");
	else
		ui->showMySky_pathToModelEdit->setStyleSheet("color:red;");
}

void AtmosphereDialog::setCurrentValues()
{
	const auto mgr = StelApp::getInstance().getStelPropertyManager();
	const auto currentModel = mgr->getProperty(MODEL_PROPERTY)->getValue().toString();
	for(int i = 0; i < ui->atmosphereModel->count(); ++i)
	{
		if(ui->atmosphereModel->itemText(i).toLower()==currentModel.toLower())
		{
			ui->atmosphereModel->setCurrentIndex(i);
			break;
		}
	}
	onModelChoiceChanged(ui->atmosphereModel->currentText());

	const auto currentModelPath = mgr->getProperty(MODEL_PATH_PROPERTY)->getValue().toString();
	ui->showMySky_pathToModelEdit->setText(currentModelPath);

	ui->showMySky_debugOptionsEnabled->setChecked(false);
}
