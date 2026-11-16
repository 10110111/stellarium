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
#include "StelMainView.hpp"
#include "StelModuleMgr.hpp"
#include "StelUtils.hpp"
#include "ui_sunspotsDialog.h"
#include <QCheckBox>
#include <QHBoxLayout>
#include <QLabel>
#include <QLocale>
#include <QMessageBox>
#include <QProgressBar>
#include <QPushButton>
#include <QSignalBlocker>
#include <QTreeWidgetItem>
#include <algorithm>

namespace
{
enum CacheColumn { DateColumn, SourceColumn, StatusColumn };

class CachedImageItem : public QTreeWidgetItem
{
public:
	using QTreeWidgetItem::QTreeWidgetItem;

	bool operator<(const QTreeWidgetItem& other) const override
	{
		if (treeWidget()->sortColumn() == DateColumn)
		{
			// Sort by the observation instant, independently of the localized date format.
			const auto time = data(DateColumn, Qt::UserRole + 1).toDouble();
			const auto otherTime = other.data(DateColumn, Qt::UserRole + 1).toDouble();
			if (time != otherTime) return time < otherTime;
			return data(DateColumn, Qt::UserRole).toString() <
			       other.data(DateColumn, Qt::UserRole).toString();
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

bool SunspotsDialog::confirmCacheRemoval(bool keepDisplayedImage, bool keepJobImages)
{
	QString message = keepDisplayedImage ? q_("Delete all cached images except the currently displayed image?") :
	                                       q_("Delete all cached images?");
	if (keepJobImages) message += "\n\n" + q_("Files needed by the running operation will be kept.");
	message += "\n\n" + q_("Deleted FITS, JPEG2000, PNG, and CSV files cannot be recovered.");
	return QMessageBox::question(&StelMainView::getInstance(), q_("Clear cache"), message,
	                             QMessageBox::Yes | QMessageBox::Cancel, QMessageBox::Cancel) == QMessageBox::Yes;
}

bool SunspotsDialog::confirmImageRemoval(bool activeJob)
{
	const QString message = activeJob
		? q_("Cancel downloading or texture preparation and delete all files for this entry?")
		: q_("Delete this entry and all its files?");
	return QMessageBox::question(&StelMainView::getInstance(), q_("Delete cached image"), message,
	                             QMessageBox::Yes | QMessageBox::Cancel, QMessageBox::Cancel) == QMessageBox::Yes;
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
	ui_->cachedImageList->sortItems(DateColumn, Qt::DescendingOrder);

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
	ui_->cacheEmptyError->hide();
	connectBoolProperty(ui_->textureEnabled, "Sunspots.textureEnabled");
	connectBoolProperty(ui_->autoDownloadCB, "Sunspots.autoDownload");
	connectBoolProperty(ui_->continuousUseCacheOnly, "Sunspots.continuousUseCacheOnly");
	connectBoolProperty(ui_->deleteSourceAfterMakingTexture, "Sunspots.deleteSourceAfterMakingTexture");
	connectIntProperty(ui_->maxImageAgeSpinBox, "Sunspots.maxImageAgeMinutes");
	connect(ui_->hmiRadioButton, &QRadioButton::clicked, this,
	        [this] { sunspots_->setInstrument(Sunspots::Instrument::HMI); });
	connect(ui_->phiRadioButton, &QRadioButton::clicked, this,
	        [this] { sunspots_->setInstrument(Sunspots::Instrument::PHI); });
	connect(sunspots_, &Sunspots::instrumentChanged, this, &SunspotsDialog::updateStatus);
	connect(ui_->downloadForSimTimeBtn, &QPushButton::clicked, sunspots_, &Sunspots::downloadForSimTime);
	connect(ui_->downloadLatestBtn, &QPushButton::clicked, sunspots_, &Sunspots::downloadLatest);
	connect(ui_->cancelButton, &QPushButton::clicked, sunspots_, &Sunspots::cancelCurrentOperation);
	connect(ui_->clearCacheBtn, &QPushButton::clicked, this, [this]
	{
		sunspots_->clearCache();
		updateEntryButtons();
	});
	connect(ui_->cacheListFilter, &QLineEdit::textChanged, this, &SunspotsDialog::filterCacheList);
	connect(ui_->onlyCurrentInstrument, &QCheckBox::toggled, this, &SunspotsDialog::filterCacheList);
	connect(ui_->cachedImageList, &QTreeWidget::itemSelectionChanged, this, &SunspotsDialog::updateEntryButtons);
	connect(ui_->useEntryBtn, &QPushButton::clicked, this, &SunspotsDialog::useSelectedEntry);
	connect(ui_->deleteEntryBtn, &QPushButton::clicked, this, [this]
	{
		const auto selected = ui_->cachedImageList->selectedItems();
		if (selected.isEmpty()) return;
		// Pass the identity by value: confirmation can change or remove the selected row.
		sunspots_->deleteCachedImage(selected.first()->data(DateColumn, Qt::UserRole).toString());
		updateEntryButtons();
	});
	connect(ui_->goToEntryTimeBtn, &QPushButton::clicked, this, [this]
	{
		const auto selected = ui_->cachedImageList->selectedItems();
		if (selected.isEmpty()) return;
		StelApp::getInstance().getCore()->setJDE(selected.first()->data(DateColumn, Qt::UserRole + 1).toDouble());
	});
	connect(ui_->cachedImageList, &QTreeWidget::itemDoubleClicked, this, &SunspotsDialog::useSelectedEntry);
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
	ui_->instrumentGroup->setEnabled(!busy);
	const QSignalBlocker blockHmi(ui_->hmiRadioButton), blockPhi(ui_->phiRadioButton);
	ui_->hmiRadioButton->setChecked(sunspots_->getInstrument() == Sunspots::Instrument::HMI);
	ui_->phiRadioButton->setChecked(sunspots_->getInstrument() == Sunspots::Instrument::PHI);
	ui_->textureEnabled->setEnabled(!sunspots_->isCanceling());
	ui_->maxImageAgeSpinBox->setEnabled(!busy);
	for (auto button : {ui_->downloadForSimTimeBtn, ui_->downloadLatestBtn})
	{
		button->setEnabled(!busy);
		button->setToolTip(busy ? q_("An image update is in progress.") : QString());
	}
	ui_->progressBar->setVisible(fetching);
	ui_->cancelButton->setVisible(fetching);
	ui_->cancelButton->setEnabled(!sunspots_->isCanceling());
	ui_->progressBar->setRange(0, 100);
	ui_->progressBar->setFormat(sunspots_->progressText());
	ui_->progressBar->setValue(sunspots_->progressValue());
	ui_->errorLabel->setText(sunspots_->errorText());
	ui_->errorLabel->setVisible(!sunspots_->errorText().isEmpty());
	ui_->cacheEmptyError->setVisible(sunspots_->hasEmptyCacheError());
	updateImageTimes();
}

void SunspotsDialog::updateImageTimes()
{
	if (!dialog) return;
	const auto core = StelApp::getInstance().getCore();
	const bool realTime = core->getIsTimeNow() && core->getRealTimeSpeed();
	ui_->downloadForSimTimeBtn->setVisible(!realTime && core->getJD() < StelUtils::getJDFromSystem());
	ui_->groupBox->setTitle(sunspots_->displayedImageTitle());
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
		items.insert(item->data(DateColumn, Qt::UserRole).toString(), item);
	}
	const QString displayed = sunspots_->displayedCacheBaseName();
	QTreeWidgetItem* insertedActiveItem = nullptr;
	for (const auto& record : records)
	{
		const QString baseName = record.cacheBaseName();
		auto item = items.take(baseName);
		if (!item)
		{
			item = new CachedImageItem(ui_->cachedImageList);
			if (baseName == activeBaseName) insertedActiveItem = item;
		}
		item->setData(DateColumn, Qt::UserRole, baseName);
		item->setData(DateColumn, Qt::UserRole + 1, record.julianDayTT());
		item->setText(DateColumn, sunspots_->observationTimeText(record));
		item->setToolTip(DateColumn, baseName);
		item->setText(SourceColumn, record.sourceName());
		item->setToolTip(SourceColumn, ui_->cachedImageList->headerItem()->toolTip(SourceColumn));
		item->setData(SourceColumn, Qt::UserRole, static_cast<int>(record.instrument));
		if (baseName == activeBaseName)
		{
			auto widget = ui_->cachedImageList->itemWidget(item, StatusColumn);
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
				ui_->cachedImageList->setItemWidget(item, StatusColumn, widget);
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
			item->setSizeHint(StatusColumn, widget->sizeHint());
			// Keep searchable text even though a widget occupies this column.
			QString searchable = status;
			searchable.replace("%p%", QString::number(sunspots_->progressValue()) + "%");
			item->setText(StatusColumn, searchable);
		}
		else
		{
			if (ui_->cachedImageList->itemWidget(item, StatusColumn))
				ui_->cachedImageList->removeItemWidget(item, StatusColumn);
			item->setData(StatusColumn, Qt::SizeHintRole, QVariant());
			item->setText(StatusColumn, baseName == displayed ? q_("Displayed") :
			                           record.hasCachedPng ? q_("Ready") : q_("Downloaded"));
		}
	}
	qDeleteAll(items);
	ui_->cachedImageList->setSortingEnabled(sorting);
	const int count = sunspots_->cachedRecords().size();
	const QString size = QLocale().formattedDataSize(sunspots_->cacheSizeBytes(), 1, QLocale::DataSizeIecFormat);
	ui_->cacheSizeLabel->setText(size);
	ui_->cacheEmptyError->setVisible(sunspots_->hasEmptyCacheError());
	ui_->cacheSummaryLabel->setText(count == 0 ? q_("No images cached") :
	                               q_("%1 images cached (%2)").arg(count).arg(size));
	filterCacheList();
	if (firstInsertion)
	{
		ui_->cachedImageList->resizeColumnToContents(DateColumn);
		ui_->cachedImageList->resizeColumnToContents(SourceColumn);
	}
	if (insertedActiveItem) ui_->cachedImageList->scrollToItem(insertedActiveItem);
}

void SunspotsDialog::filterCacheList()
{
	const QString filter = ui_->cacheListFilter->text();
	const bool onlyCurrentInstrument = ui_->onlyCurrentInstrument->isChecked();
	const int instrument = static_cast<int>(sunspots_->getInstrument());
	for (int i = 0; i < ui_->cachedImageList->topLevelItemCount(); ++i)
	{
		auto item = ui_->cachedImageList->topLevelItem(i);
		const bool matchesInstrument = !onlyCurrentInstrument ||
		                               item->data(SourceColumn, Qt::UserRole).toInt() == instrument;
		const bool matchesText = item->text(DateColumn).contains(filter, Qt::CaseInsensitive) ||
		                         item->text(SourceColumn).contains(filter, Qt::CaseInsensitive) ||
		                         item->text(StatusColumn).contains(filter, Qt::CaseInsensitive);
		item->setHidden(!matchesInstrument || !matchesText);
	}
	updateEntryButtons();
}

void SunspotsDialog::updateEntryButtons()
{
	const auto selected = ui_->cachedImageList->selectedItems();
	const QString noSelection = q_("No image is selected in the cache list.");
	ui_->goToEntryTimeBtn->setEnabled(!selected.isEmpty());
	ui_->goToEntryTimeBtn->setToolTip(selected.isEmpty() ? noSelection : QString());
	QString reason = noSelection;
	const bool canUse = !selected.isEmpty() &&
	    sunspots_->canDisplayCachedImage(selected.first()->data(DateColumn, Qt::UserRole).toString(), &reason);
	ui_->useEntryBtn->setEnabled(canUse);
	ui_->useEntryBtn->setToolTip(canUse ? QString() : reason);
	reason = noSelection;
	const QString baseName = selected.isEmpty() ? QString() :
	                        selected.first()->data(DateColumn, Qt::UserRole).toString();
	const bool canDelete = !selected.isEmpty() && sunspots_->canDeleteCachedImage(baseName, &reason);
	ui_->deleteEntryBtn->setEnabled(canDelete);
	const auto active = sunspots_->activeRecord();
	ui_->deleteEntryBtn->setToolTip(!canDelete ? reason :
	    active.time.isValid() && active.cacheBaseName() == baseName ? q_("Cancel and delete cached image") :
	                                                               q_("Delete cached image"));
}

void SunspotsDialog::useSelectedEntry()
{
	const auto selected = ui_->cachedImageList->selectedItems();
	if (selected.isEmpty()) return;
	// Copy the identity before the cache refresh can rebuild the list.
	const QString baseName = selected.first()->data(DateColumn, Qt::UserRole).toString();
	sunspots_->displayCachedImage(baseName);
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

	html += "<p>" + q_("Downloads HMI or PHI continuum images for the selected simulation time "
	                  "and displays them on the Sun.") + "</p>";

	html += "<h3>" + q_("Data sources and products") + "</h3>";
	html += "<p>" + q_("HMI observes the Sun from NASA's Solar Dynamics Observatory (SDO). "
	                  "The plugin obtains its continuum FITS images from the Joint Science Operations Center (JSOC).")
	        + "</p><ul>";
	html += "<li>" + q_("HMI: science data from the hmi.Ic_noLimbDark_720s series, "
	                   "with limb darkening removed.") + "</li>";
	html += "<li>" + q_("HMI (NRT): preliminary Near Real-Time data from the hmi.Ic_noLimbDark_720s_nrt series, "
	                   "used when the science series does not cover the requested time.") + "</li></ul>";

	html += "<p>" + q_("PHI's Full Disk Telescope observes the Sun from Solar Orbiter, providing a different "
	                  "viewpoint that can include regions hidden from Earth. The plugin obtains PHI continuum "
	                  "FITS images from ESA's Solar Orbiter Archive (SOAR), and JPEG2000 quicklooks from the "
	                  "Royal Observatory of Belgium's Helioviewer service.") + "</p><ul>";
	html += "<li>" + q_("PHI: Level 2 (L2) science FITS data from SOAR.") + "</li>";
	html += "<li>" + q_("PHI (LL): Low Latency (LL02) FITS data from SOAR, used when L2 data does not cover "
	                   "the requested time.") + "</li>";
	html += "<li>" + q_("PHI (JP2): Low Latency (LL03) JPEG2000 quicklook images from Helioviewer, used when "
	                   "neither FITS product covers the requested time. Intensity information lost through "
	                   "clipping and lossy compression cannot be recovered.") + "</li></ul>";
	html += "<p>" + q_("These names appear in the Source column. Download latest compares the available products "
	                  "for the selected instrument and chooses its newest observation.") + "</p>";

	html += StelApp::getInstance().getModuleMgr().getStandardSupportLinksInfo("Sunspots plugin");
	html += "</body></html>";

	ui_->about->setHtml(html);
}
