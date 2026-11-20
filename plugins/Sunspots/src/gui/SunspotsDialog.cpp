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

#include "SunspotsDialog.hpp"
#include "Sunspots.hpp"
#include "StelApp.hpp"
#include "StelCore.hpp"
#include "StelGui.hpp"
#include "StelModuleMgr.hpp"
#include "StelUtils.hpp"
#include "ui_sunspotsDialog.h"
#include <QHBoxLayout>
#include <QLabel>
#include <QLocale>
#include <QProgressBar>
#include <QPushButton>
#include <QTreeWidgetItem>
#include <algorithm>

namespace
{
class CachedImageItem : public QTreeWidgetItem
{
public:
	using QTreeWidgetItem::QTreeWidgetItem;

	bool operator<(const QTreeWidgetItem& other) const override
	{
		if (treeWidget()->sortColumn() == 0)
		{
			// Sort by the observation instant, independently of the localized date format.
			const auto time = data(0, Qt::UserRole + 1).toDateTime();
			const auto otherTime = other.data(0, Qt::UserRole + 1).toDateTime();
			if (time != otherTime) return time < otherTime;
			return data(0, Qt::UserRole).toString() < other.data(0, Qt::UserRole).toString();
		}
		return QTreeWidgetItem::operator<(other);
	}
};
}

SunspotsDialog::SunspotsDialog(Sunspots* sunspots)
	: StelDialog("Sunspots")
	, sunspots_(sunspots)
	, ui_(new Ui_SunspotsDialog{})
{
	presentationTimer_.setInterval(1000);
	connect(&presentationTimer_, &QTimer::timeout, this, &SunspotsDialog::updateImageTimes);
	timeChangeTimer_.setSingleShot(true);
	connect(&timeChangeTimer_, &QTimer::timeout, this, &SunspotsDialog::updateImageTimes);
	connect(this, &StelDialog::visibleChanged, this, [this](bool visible)
	{
		if (visible)
		{
			sunspots_->refreshCache();
			updateStatus();
			presentationTimer_.start();
		}
		else presentationTimer_.stop();
	});
}

SunspotsDialog::~SunspotsDialog()
{
}

void SunspotsDialog::retranslate()
{
	if(dialog)
	{
		ui_->retranslateUi(dialog);
		setAboutText();
		updateStatus();
	}
}

void SunspotsDialog::createDialogContent()
{
	ui_->setupUi(dialog);
	ui_->tabs->setCurrentIndex(0);
	ui_->cachedImageList->sortItems(0, Qt::DescendingOrder);

	// Kinetic scrolling
	kineticScrollingList << ui_->about;
	StelGui* gui= dynamic_cast<StelGui*>(StelApp::getInstance().getGui());
	if(gui)
	{
		enableKineticScrolling(gui->getFlagUseKineticScrolling());
		connect(gui, &StelGui::flagUseKineticScrollingChanged, this, &SunspotsDialog::enableKineticScrolling);
	}

	connect(&StelApp::getInstance(), &StelApp::languageChanged, this, &SunspotsDialog::retranslate);
	connect(ui_->titleBar, &TitleBar::closeClicked, this, &StelDialog::close);
	connect(ui_->titleBar, &TitleBar::movedTo, this, &SunspotsDialog::handleMovedTo);

	// Images and settings tabs
	ui_->progressBar->hide();
	ui_->cancelButton->hide();
	ui_->errorLabel->hide();
	connectBoolProperty(ui_->hmiPhotosphereEnabled, "Sunspots.hmiPhotosphereEnabled");
	connectBoolProperty(ui_->autoDownloadCB, "Sunspots.autoDownload");
	connectIntProperty(ui_->maxImageAgeSpinBox, "Sunspots.maxImageAgeMinutes");
	connect(ui_->downloadForSimTimeBtn, &QPushButton::clicked, sunspots_, &Sunspots::downloadForSimTime);
	connect(ui_->downloadLatestBtn, &QPushButton::clicked, sunspots_, &Sunspots::downloadLatest);
	connect(ui_->cancelButton, &QPushButton::clicked, sunspots_, &Sunspots::cancelCurrentOperation);
	connect(ui_->cacheListFilter, &QLineEdit::textChanged, this, &SunspotsDialog::filterCacheList);
	connect(ui_->cachedImageList, &QTreeWidget::itemDoubleClicked, this, [this](QTreeWidgetItem* item, int)
	{
		// Copy the identity before the cache refresh can rebuild the list.
		const QString baseName = item->data(0, Qt::UserRole).toString();
		sunspots_->displayCachedImage(baseName);
	});
	connect(sunspots_, &Sunspots::statusChanged, this, &SunspotsDialog::updateStatus);
	connect(sunspots_, &Sunspots::cacheChanged, this, &SunspotsDialog::updateCacheList);
	const auto core = StelApp::getInstance().getCore();
	connect(core, &StelCore::currentTimeZoneChanged, this, &SunspotsDialog::updateImageTimes);
	// A time-sync signal can expose an intermediate clock while planet information is being
	// calculated. Read the restored clock after the calculation returns, once per event-loop turn.
	connect(core, &StelCore::timeSyncOccurred, &timeChangeTimer_, qOverload<>(&QTimer::start));
	connect(core, &StelCore::timeRateChanged, &timeChangeTimer_, qOverload<>(&QTimer::start));
	updateStatus();

	// About tab
	setAboutText();
	if(gui)
	{
		ui_->about->document()->setDefaultStyleSheet(QString(gui->getStelStyle().htmlStyleSheet));
	}
}

void SunspotsDialog::updateStatus()
{
	const bool busy = sunspots_->isBusy();
	const bool fetching = sunspots_->progressKind() == Sunspots::ProgressKind::Fetching;
	ui_->hmiPhotosphereEnabled->setEnabled(!sunspots_->isCanceling());
	ui_->autoDownloadCB->setEnabled(!busy);
	ui_->maxImageAgeSpinBox->setEnabled(!busy);
	ui_->downloadForSimTimeBtn->setEnabled(!busy);
	ui_->downloadLatestBtn->setEnabled(!busy);
	ui_->manageCacheBtn->setEnabled(false);
	ui_->manageCacheBtn->setToolTip(q_("Not implemented yet"));
	ui_->findImagesBtn->setEnabled(false);
	ui_->findImagesBtn->setToolTip(q_("Not implemented yet"));
	ui_->progressBar->setVisible(fetching);
	ui_->cancelButton->setVisible(fetching);
	ui_->cancelButton->setEnabled(!sunspots_->isCanceling());
	ui_->progressBar->setRange(0, 100);
	ui_->progressBar->setFormat(sunspots_->progressText());
	ui_->progressBar->setValue(sunspots_->progressValue());
	ui_->errorLabel->setText(sunspots_->errorText());
	ui_->errorLabel->setVisible(!sunspots_->errorText().isEmpty());
	updateImageTimes();
}

void SunspotsDialog::updateImageTimes()
{
	if (!dialog) return;
	const auto core = StelApp::getInstance().getCore();
	const bool realTime = core->getIsTimeNow() && core->getRealTimeSpeed();
	ui_->downloadForSimTimeBtn->setVisible(!realTime && core->getJD() < StelUtils::getJDFromSystem());
	ui_->displayedTime->setText(sunspots_->displayedTimeText());
	ui_->displayedImageAge->setText(sunspots_->displayedImageAgeText());
	updateCacheList();
}

void SunspotsDialog::updateCacheList()
{
	QList<Sunspots::Record> records = sunspots_->cachedRecords();
	const auto active = sunspots_->activeRecord();
	const QString activeBaseName = active.time.isValid() ? active.cacheBaseName() : QString();
	if (active.time.isValid() && std::none_of(records.cbegin(), records.cend(), [&activeBaseName](const auto& record)
	    { return record.cacheBaseName() == activeBaseName; })) records.append(active);
	const bool firstInsertion = ui_->cachedImageList->topLevelItemCount() == 0 && !records.isEmpty();
	// Match rows by identity rather than position so sorting never replaces their
	// progress controls. Suspend sorting while changing item data and restore the
	// user's selected column and direction once the update is complete.
	const bool sorting = ui_->cachedImageList->isSortingEnabled();
	ui_->cachedImageList->setSortingEnabled(false);
	QMap<QString, QTreeWidgetItem*> items;
	for (int i = 0; i < ui_->cachedImageList->topLevelItemCount(); ++i)
	{
		auto item = ui_->cachedImageList->topLevelItem(i);
		items.insert(item->data(0, Qt::UserRole).toString(), item);
	}
	const QString displayed = sunspots_->displayedCacheBaseName();
	for (const auto& record : records)
	{
		const QString baseName = record.cacheBaseName();
		auto item = items.take(baseName);
		if (!item) item = new CachedImageItem(ui_->cachedImageList);
		item->setData(0, Qt::UserRole, baseName);
		item->setData(0, Qt::UserRole + 1, record.time);
		item->setText(0, sunspots_->observationTimeText(record));
		item->setToolTip(0, baseName);
		if (baseName == activeBaseName)
		{
			auto widget = ui_->cachedImageList->itemWidget(item, 1);
			if (!widget)
			{
				widget = new QWidget(ui_->cachedImageList);
				widget->setAutoFillBackground(true);
				auto layout = new QHBoxLayout(widget);
				layout->setContentsMargins(0, 0, 0, 0);
				auto progress = new QProgressBar(widget);
				progress->setObjectName("imageProgress");
				layout->addWidget(progress, 1);
				auto stage = new QLabel(widget);
				stage->setObjectName("imageStage");
				layout->addWidget(stage, 1);
				auto cancel = new QPushButton(widget);
				cancel->setObjectName("cancelImageBtn");
				layout->addWidget(cancel);
				connect(cancel, &QPushButton::clicked, sunspots_, &Sunspots::cancelCurrentOperation);
				ui_->cachedImageList->setItemWidget(item, 1, widget);
			}
			auto progress = widget->findChild<QProgressBar*>("imageProgress");
			auto stage = widget->findChild<QLabel*>("imageStage");
			auto cancel = widget->findChild<QPushButton*>("cancelImageBtn");
			const bool saving = sunspots_->progressKind() == Sunspots::ProgressKind::Saving;
			const bool loading = sunspots_->progressKind() == Sunspots::ProgressKind::Loading;
			const bool textOnly = saving || loading || sunspots_->isCanceling();
			const QString status = saving ? q_("Saving...") : sunspots_->progressText();
			progress->setVisible(!textOnly);
			progress->setFormat(status);
			progress->setValue(sunspots_->progressValue());
			stage->setVisible(textOnly);
			stage->setText(status);
			cancel->setText(q_("Cancel"));
			cancel->setEnabled(!sunspots_->isCanceling());
			item->setSizeHint(1, widget->sizeHint());
			// Keep searchable text even though a widget occupies this column.
			QString searchable = status;
			searchable.replace("%p%", QString::number(sunspots_->progressValue()) + "%");
			item->setText(1, searchable);
		}
		else
		{
			if (ui_->cachedImageList->itemWidget(item, 1)) ui_->cachedImageList->removeItemWidget(item, 1);
			item->setData(1, Qt::SizeHintRole, QVariant());
			item->setText(1, baseName == displayed ? q_("Displayed") :
			                 record.hasCachedPng ? q_("Ready") : q_("Downloaded"));
		}
	}
	qDeleteAll(items);
	ui_->cachedImageList->setSortingEnabled(sorting);
	const int count = sunspots_->cachedRecords().size();
	const QString size = QLocale().formattedDataSize(sunspots_->cacheSizeBytes(), 1, QLocale::DataSizeIecFormat);
	ui_->cacheSummaryLabel->setText(count == 0 ? q_("No images cached") :
	                               q_("%1 images cached (%2)").arg(count).arg(size));
	filterCacheList();
	if (firstInsertion) ui_->cachedImageList->resizeColumnToContents(0);
}

void SunspotsDialog::filterCacheList()
{
	const QString filter = ui_->cacheListFilter->text();
	for (int i = 0; i < ui_->cachedImageList->topLevelItemCount(); ++i)
	{
		auto item = ui_->cachedImageList->topLevelItem(i);
		item->setHidden(!item->text(0).contains(filter, Qt::CaseInsensitive) &&
		                !item->text(1).contains(filter, Qt::CaseInsensitive));
	}
}

void SunspotsDialog::setAboutText()
{
	QString html = "<html><head></head><body>"
	"<h2>" + q_("Sunspots Plug-in") + "</h2>"
	"<table class='layout' width=\"90%\">"
		"<tr width=\"30%\"><td><strong>" + q_("Version") + ":</strong></td><td>" + SUNSPOTS_PLUGIN_VERSION + "</td></tr>"
		"<tr><td><strong>" + q_("License") + ":</strong></td><td>" + SUNSPOTS_PLUGIN_LICENSE + "</td></tr>"
		"<tr><td><strong>" + q_("Author") + ":</strong></td><td>Ruslan Kabatsayev &lt;b7.10110111+stellarium@gmail.com&gt;</td></tr>"
	"</table>";
	// Overview
	html += "<h3>" + q_("Overview") + "</h3>";

	html += "<p>" + q_("Downloads HMI continuum images for the selected simulation time and displays them on the Sun.") + "</p>";

	html += StelApp::getInstance().getModuleMgr().getStandardSupportLinksInfo("Sunspots plugin");
	html += "</body></html>";

	ui_->about->setHtml(html);
}
