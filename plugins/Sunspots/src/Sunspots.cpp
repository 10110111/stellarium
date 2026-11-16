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

#include "Sunspots.hpp"
#include "SunspotsDialog.hpp"
#include "StelApp.hpp"
#include "StelGui.hpp"
#include "StelGuiItems.hpp"
#include "StelLocaleMgr.hpp"
#include "StelModuleMgr.hpp"
#include "StelObjectMgr.hpp"
#include "StelProgressController.hpp"
#include "StelTranslator.hpp"
#include "StelCore.hpp"
#include "StelFileMgr.hpp"
#include "StelUtils.hpp"
#include "SolarSystem.hpp"
#include "SunspotsImage.hpp"
#include "VecMath.hpp"
#include <QDebug>
#include <QDir>
#include <QFile>
#include <QFileInfo>
#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonObject>
#include <QNetworkReply>
#include <QNetworkRequest>
#include <QSaveFile>
#include <QTimer>
#include <QUrlQuery>
#include <QRunnable>
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <fitsio.h>

Q_LOGGING_CATEGORY(sunspots, "stel.plugin.Sunspots", QtInfoMsg)

namespace
{
const QString finalSeries = QStringLiteral("hmi.Ic_noLimbDark_720s");
const QString nrtSeries = QStringLiteral("hmi.Ic_noLimbDark_720s_nrt");

QDateTime parseObservationTime(const QString& value)
{
	QDateTime time = QDateTime::fromString(value.left(19), "yyyy.MM.dd_HH:mm:ss");
	// T_OBS contains TAI calendar fields, but QDateTime has no TAI time specification.
	// Mark them as UTC solely to make QDateTime comparisons and StelUtils::qDateTimeToJd()
	// use these fields consistently. This does not convert TAI to UTC: the QDateTime
	// represents a fictitious UTC instant and must not be used as a real UTC timestamp.
	time.setTimeSpec(Qt::UTC);
	return time;
}

QString observationRange(const QList<Sunspots::Record>& records)
{
	QDateTime first, last;
	for (const Sunspots::Record& record : records)
	{
		if (!record.time.isValid()) continue;
		if (!first.isValid() || record.time < first) first = record.time;
		if (!last.isValid() || record.time > last) last = record.time;
	}
	if (!first.isValid()) return "(none)";
	return first.toString("yyyy.MM.dd_HH:mm:ss") + "_TAI .. " +
	       last.toString("yyyy.MM.dd_HH:mm:ss") + "_TAI";
}

QString observationOrNone(const Sunspots::Record& record)
{
	return record.time.isValid() ? record.observation : QStringLiteral("(none)");
}

QString observationTimeString(double jde)
{
	const QDateTime tai = StelUtils::jdToQDateTime(StelUtils::ttToTai(jde), Qt::UTC);
	return tai.toString("yyyy.MM.dd_HH:mm:ss") + "_TAI";
}

QList<Sunspots::Record> parseRecords(const QByteArray& data, const QString& series, QString* error)
{
	QJsonParseError parseError;
	const QJsonDocument document = QJsonDocument::fromJson(data, &parseError);
	if (parseError.error != QJsonParseError::NoError || !document.isObject())
	{
		*error = q_("Invalid JSOC response: %1").arg(parseError.errorString());
		return {};
	}
	const QJsonObject root = document.object();
	if (root.value("status").toInt(-1) != 0)
	{
		*error = q_("JSOC query failed: %1").arg(root.value("error").toString());
		return {};
	}
	QMap<QString, QJsonArray> columns;
	for (const QJsonValue& value : root.value("keywords").toArray())
	{
		const QJsonObject column = value.toObject();
		columns.insert(column.value("name").toString(), column.value("values").toArray());
	}
	for (const QJsonValue& value : root.value("segments").toArray())
	{
		const QJsonObject column = value.toObject();
		columns.insert(column.value("name").toString(), column.value("values").toArray());
	}
	QList<Sunspots::Record> records;
	const int count = columns.value("T_OBS").size();
	for (int i = 0; i < count; ++i)
	{
		Sunspots::Record record;
		record.series = series;
		for (auto it = columns.cbegin(); it != columns.cend(); ++it)
		{
			const QJsonValue value = it.value().at(i);
			record.metadata.insert(it.key(), value.isString() ? value.toString() : value.toVariant().toString());
		}
		record.observation = record.metadata.value("T_OBS");
		record.time = parseObservationTime(record.observation);
		const QString segmentPath = record.metadata.value("continuum");
		if (!record.time.isValid() || segmentPath.isEmpty() || segmentPath == "NoDataFile")
		{
			qCWarning(sunspots).noquote() << "Skipping JSOC record in" << series << "at"
			                              << observationOrNone(record) << ":"
			                              << (!record.time.isValid() ? "invalid T_OBS" : "missing continuum segment");
			continue;
		}
		const QUrl segmentUrl(segmentPath);
		record.segment = segmentUrl.isRelative() ? QUrl("https://jsoc1.stanford.edu/").resolved(segmentUrl) : segmentUrl;
		if (record.segment.scheme() != "https" && record.segment.scheme() != "http")
		{
			qCWarning(sunspots).noquote() << "Skipping JSOC record in" << series << "at" << record.observation
			                              << ": invalid continuum URL" << segmentPath;
			continue;
		}
		record.metadata.remove("continuum");
		record.metadata.insert("SERIES", series);
		records.append(record);
	}
	return records;
}
}

QString Sunspots::Record::cacheBaseName() const
{
	QString timePart = observation.left(19);
	timePart.replace(':', '-');
	return timePart + "_" + series;
}

StelModule* SunspotsStelPluginInterface::getStelModule() const
{
	return new Sunspots();
}

StelPluginInfo SunspotsStelPluginInterface::getPluginInfo() const
{
	// Allow to load the resources when used as a static plugin
	Q_INIT_RESOURCE(Sunspots);

	StelPluginInfo info;
	info.id = "Sunspots";
	info.displayedName = N_("Sunspots");
	info.authors = "Ruslan Kabatsayev";
	info.contact = STELLARIUM_DEV_URL;
	info.description = N_("Replaces the Sun texture with an HMI image of the photosphere");
	info.version = SUNSPOTS_PLUGIN_VERSION;
	info.license = SUNSPOTS_PLUGIN_LICENSE;
	return info;
}

Sunspots::Sunspots()
	: conf_(StelApp::getInstance().getSettings())
{
	setObjectName("Sunspots");
	workers_.setMaxThreadCount(1);

	// We don't depend on cfitsio being reentrant because we only use one parsing job at a time,
	// but if another part of Stellarium ever wants to use cfitsio too, this may be significant.
	qCInfo(sunspots).nospace() << "cfitsio is " << (fits_is_reentrant() ? "" : "NOT ") << "reentrant";

	autoUpdateTimer_.setSingleShot(true);
	connect(&autoUpdateTimer_, &QTimer::timeout, this, &Sunspots::pollTime);
	timeChangeTimer_.setSingleShot(true);
	connect(&timeChangeTimer_, &QTimer::timeout, this, &Sunspots::pollTime);
}

Sunspots::~Sunspots()
{
	cancellationToken_->store(true);
	for (QNetworkReply* reply : queryReplies_)
	{
		reply->disconnect(this);
		reply->abort();
		reply->deleteLater();
	}
	if (downloadReply_)
	{
		downloadReply_->disconnect(this);
		downloadReply_->abort();
		downloadReply_->deleteLater();
	}
	removeDownloadProgressBar();
	workers_.waitForDone();
	if (hmiPhotosphereEnabled)
		GETSTELMODULE(SolarSystem)->resetTextures("Sun");
}

void Sunspots::init()
{
	dialog_.reset(new SunspotsDialog(this));

	loadSettings();
	try
	{
		const QString cacheDir = StelFileMgr::getUserDir() + "/modules/Sunspots/cache";
		StelFileMgr::makeSureDirExistsAndIsWritable(cacheDir);
		cacheDir_ = cacheDir;
	}
	catch (const std::runtime_error& error)
	{
		showError(QString::fromUtf8(error.what()));
	}
	refreshCache();

	addAction("actionShow_Sunspots", N_("Sunspots"), N_("Show the photosphere from HMI"), "hmiPhotosphereEnabled");
	addAction("actionShow_Sunspots_dialog", N_("Sunspots"), N_("Show Sunspots dialog"), "dialogVisible");
	// Make sure that opening/closing the dialog by other means than dialogVisible property is reflected on this property
	connect(dialog_.get(), &StelDialog::visibleChanged, this, [this](const bool visible)
	        {
	            if(visible != dialogVisible)
	            {
	                dialogVisible = visible;
	                emit dialogToggled(visible);
	            }
	        });

	toolbarButton_ = new StelButton(nullptr,
	                                QPixmap(":/Sunspots/bt_Sunspots_On.png"),
	                                QPixmap(":/Sunspots/bt_Sunspots_Off.png"),
	                                QPixmap(":/graphicGui/miscGlow32x32.png"),
	                                "actionShow_Sunspots",
	                                false,
	                                "actionShow_Sunspots_dialog");
	// Cancellation cannot resume; undo a click made while the worker is stopping.
	connect(toolbarButton_, &StelButton::triggered, this, [this]
	{
		if (isCanceling()) toolbarButton_->setChecked(hmiPhotosphereEnabled);
	});
	if(const auto gui = dynamic_cast<StelGui*>(StelApp::getInstance().getGui()))
	{
		gui->getButtonBar()->addButton(toolbarButton_, "065-pluginsGroup");
	}
	// Planet information calculations temporarily change the core time and emit these signals
	// before fully restoring the clock. Coalesce them until control returns to the event loop.
	const auto core = StelApp::getInstance().getCore();
	connect(core, &StelCore::timeSyncOccurred, &timeChangeTimer_, qOverload<>(&QTimer::start));
	connect(core, &StelCore::timeRateChanged, &timeChangeTimer_, qOverload<>(&QTimer::start));
	if (hmiPhotosphereEnabled)
		QTimer::singleShot(0, this, &Sunspots::ensureUsableImage);
}

bool Sunspots::configureGui(const bool show)
{
	if(show)
	{
		dialog_->setVisible(show);
	}
	return true;
}

void Sunspots::loadSettings()
{
	// Read the previous setting name so existing preferences survive the HMI rename.
	hmiPhotosphereEnabled = conf_->value("Sunspots/show_hmi_photosphere",
	                                      conf_->value("Sunspots/show_sdo_photosphere", false)).toBool();
	autoDownload_ = conf_->value("Sunspots/auto_download", false).toBool();
	maxImageAgeMinutes_ = qBound(6, conf_->value("Sunspots/max_image_age_minutes", 180).toInt(), 9999);
}

void Sunspots::setDialogVisible(const bool enable)
{
	if(dialogVisible == enable) return;

	dialogVisible = enable;
	dialog_->setVisible(enable);
	emit dialogToggled(enable);
}

void Sunspots::setHmiPhotosphereEnabled(const bool enable)
{
	if (isCanceling() || hmiPhotosphereEnabled == enable) return;

	hmiPhotosphereEnabled = enable;
	conf_->setValue("Sunspots/show_hmi_photosphere", enable);
	emit hmiPhotosphereToggled(enable);
	if (enable)
	{
		if (isBusy()) prepareTextureForJob_ = true;
		else ensureUsableImage();
	}
	else
	{
		displayedRecord_ = {};
		GETSTELMODULE(SolarSystem)->resetTextures("Sun");
		emit statusChanged();
		if (!isBusy()) showError({});
	}
}

void Sunspots::setAutoDownload(bool enable)
{
	// Disabling automatic updates is safe during a job and lets Cancel stop future polling.
	if ((isBusy() && enable) || autoDownload_ == enable) return;
	autoDownload_ = enable;
	conf_->setValue("Sunspots/auto_download", enable);
	emit autoDownloadChanged(enable);
	if (enable && hmiPhotosphereEnabled) pollTime();
	else if (!enable && tooFastWarning_) showError({});
}

void Sunspots::setMaxImageAgeMinutes(int minutes)
{
	minutes = qBound(6, minutes, 9999);
	if (maxImageAgeMinutes_ == minutes) return;
	maxImageAgeMinutes_ = minutes;
	conf_->setValue("Sunspots/max_image_age_minutes", minutes);
	emit maxImageAgeMinutesChanged(minutes);
	if (autoDownload_ && hmiPhotosphereEnabled) pollTime();
}

QString Sunspots::observationTimeText(const Record& record) const
{
	if (!record.time.isValid()) return q_("N/A");
	const auto core = StelApp::getInstance().getCore();
	const double tt = StelUtils::taiToTt(StelUtils::qDateTimeToJd(record.time));
	// Follow Stellarium's TT/UT convention rather than treating the fictitious UTC
	// QDateTime carrying TAI fields as a genuine UTC instant. Evaluate both DeltaT
	// and the time-zone offset at the observation's date, not at simulation time.
	const double ut = tt - core->computeDeltaT(tt) * StelCore::JD_SECOND;
	const double offset = core->getUTCOffset(ut);
	const auto& locale = StelApp::getInstance().getLocaleMgr();
	return locale.getPrintableDateLocal(ut, offset) + " " + locale.getPrintableTimeLocal(ut, offset);
}

QString Sunspots::displayedTimeText() const
{
	return observationTimeText(displayedRecord_);
}

QString Sunspots::displayedCacheBaseName() const
{
	return displayedRecord_.time.isValid() ? displayedRecord_.cacheBaseName() : QString();
}

QString Sunspots::displayedImageAgeText() const
{
	if (!displayedRecord_.time.isValid()) return q_("N/A");
	const auto core = StelApp::getInstance().getCore();
	const bool realTime = core->getIsTimeNow() && core->getRealTimeSpeed();
	double referenceTT = core->getJDE();
	if (realTime)
	{
		const double now = StelUtils::getJDFromSystem();
		referenceTT = now + core->computeDeltaT(now) * StelCore::JD_SECOND;
	}
	const double age = (StelUtils::ttToTai(referenceTT) - StelUtils::qDateTimeToJd(displayedRecord_.time)) * 1440.;
	if (std::abs(age) >= 1440.)
	{
		const qint64 hours = std::llround(std::abs(age) / 60.);
		if (realTime)
		{
			const QString sign = age < 0. ? QStringLiteral("-") : QString();
			return q_("%1 days %2 hours").arg(sign + QString::number(hours / 24)).arg(hours % 24);
		}
		return (age < 0. ? q_("%1 days %2 hours after sim time") : q_("%1 days %2 hours before sim time"))
		       .arg(hours / 24).arg(hours % 24);
	}
	const qint64 minutes = std::llround(std::abs(age));
	if (realTime)
	{
		const QString sign = age < 0. ? QStringLiteral("-") : QString();
		return minutes < 60 ? q_("%1 minutes").arg(sign + QString::number(minutes))
		                    : q_("%1h%2m").arg(sign + QString::number(minutes / 60)).arg(minutes % 60);
	}
	if (minutes < 60)
		return (age < 0. ? q_("%1 minutes after sim time") : q_("%1 minutes before sim time")).arg(minutes);
	return (age < 0. ? q_("%1h%2m after sim time") : q_("%1h%2m before sim time"))
	       .arg(minutes / 60).arg(minutes % 60);
}

QString Sunspots::errorText() const
{
	return tooFastWarning_ ? q_("(too fast simulation, continuous updating suppressed)") : error_;
}

QString Sunspots::progressText() const
{
	switch (progressKind_)
	{
	case ProgressKind::Fetching: return q_("Fetching data lists %p%");
	case ProgressKind::Downloading: return q_("Downloading %p%");
	case ProgressKind::PreparingTexture: return q_("Preparing texture %p%");
	case ProgressKind::Saving: return q_("Saving %p%");
	case ProgressKind::Loading: return q_("Loading...");
	case ProgressKind::Canceling: return q_("Canceling...");
	case ProgressKind::None: return {};
	}
	return {};
}

bool Sunspots::isFast() const
{
	return std::abs(StelApp::getInstance().getCore()->getTimeRate()) > StelCore::JD_SECOND;
}

bool Sunspots::imageNeedsUpdate(const Record& record) const
{
	if (!record.time.isValid()) return true;
	const double simulatedTai = StelUtils::ttToTai(StelApp::getInstance().getCore()->getJDE());
	const double imageTai = StelUtils::qDateTimeToJd(record.time);
	return std::abs(simulatedTai - imageTai) * 1440. >= maxImageAgeMinutes_;
}

Sunspots::Record Sunspots::usableCachedRecord() const
{
	if (cacheDir_.isEmpty()) return {};
	if (!imageNeedsUpdate(lastReadyRecord_) && QFileInfo::exists(dataPath(lastReadyRecord_, ".png")))
		return lastReadyRecord_;

	Record best;
	const double simulatedTai = StelUtils::ttToTai(StelApp::getInstance().getCore()->getJDE());
	const QDir dir(cacheDir_);
	for (const QFileInfo& file : dir.entryInfoList({QStringLiteral("*.png"), QStringLiteral("*.fits")}, QDir::Files))
	{
		if (file.size() == 0) continue;
		const QString cacheBaseName = file.completeBaseName();
		if (cacheBaseName.size() <= 20 || cacheBaseName.at(19) != '_') continue;
		Record record;
		record.series = cacheBaseName.mid(20);
		if (record.series != finalSeries && record.series != nrtSeries) continue;
		record.observation = cacheBaseName.left(19);
		record.observation.replace('-', ':');
		record.observation += "_TAI";
		record.time = parseObservationTime(record.observation);
		if (!record.time.isValid() || record.cacheBaseName() != cacheBaseName || imageNeedsUpdate(record)) continue;
		const double distance = std::abs(StelUtils::qDateTimeToJd(record.time) - simulatedTai);
		const double bestDistance = best.time.isValid()
		                            ? std::abs(StelUtils::qDateTimeToJd(best.time) - simulatedTai) : 0.;
		if (!best.time.isValid() || distance < bestDistance ||
		    (distance == bestDistance && record.series == finalSeries && best.series != finalSeries))
			best = record;
	}
	return best;
}

void Sunspots::refreshCache()
{
	cachedRecords_.clear();
	cacheSizeBytes_ = 0;
	if (!cacheDir_.isEmpty())
	{
		QMap<QString, Record> records;
		const QDir dir(cacheDir_);
		const auto files = dir.entryInfoList({QStringLiteral("*.fits"), QStringLiteral("*.png"),
		                                      QStringLiteral("*.csv")}, QDir::Files);
		for (const QFileInfo& file : files)
		{
			cacheSizeBytes_ += file.size();
			if (file.size() == 0 || file.suffix() == "csv") continue;
			const QString baseName = file.completeBaseName();
			if (baseName.size() <= 20 || baseName.at(19) != '_') continue;
			const bool png = file.suffix() == "png";
			auto existing = records.find(baseName);
			if (existing != records.end())
			{
				existing->hasCachedPng |= png;
				continue;
			}
			Record record;
			record.series = baseName.mid(20);
			if (record.series != finalSeries && record.series != nrtSeries) continue;
			record.observation = baseName.left(19);
			record.observation.replace('-', ':');
			record.observation += "_TAI";
			record.time = parseObservationTime(record.observation);
			if (!record.time.isValid() || record.cacheBaseName() != baseName) continue;
			record.hasCachedPng = png;
			records.insert(baseName, record);
		}
		cachedRecords_ = records.values();
		std::sort(cachedRecords_.begin(), cachedRecords_.end(), [](const Record& a, const Record& b)
		{
			if (a.time != b.time) return a.time > b.time;
			return a.series < b.series;
		});
	}
	emit cacheChanged();
}

void Sunspots::pollTime()
{
	if (!hmiPhotosphereEnabled || !autoDownload_) return;
	if (isFast())
	{
		showTooFastWarning();
		return;
	}
	if (tooFastWarning_) showError({});
	if (isBusy() || autoUpdateTimer_.isActive()) return;
	if (!imageNeedsUpdate(displayedRecord_))
	{
		autoUpdateTimer_.start(60000);
		return;
	}
	ensureUsableImage();
}

void Sunspots::showError(const QString& error)
{
	if (error_ == error && !tooFastWarning_) return;
	tooFastWarning_ = false;
	error_ = error;
	emit statusChanged();
}

void Sunspots::showTooFastWarning()
{
	if (tooFastWarning_) return;
	tooFastWarning_ = true;
	error_.clear();
	emit statusChanged();
}

void Sunspots::beginJob(ProgressKind kind)
{
	Q_ASSERT(!isBusy());
	cancellationToken_ = std::make_shared<std::atomic_bool>(false);
	prepareTextureForJob_ = hmiPhotosphereEnabled;
	jobRecord_ = {};
	queryError_.clear();
	showError({});
	showProgress(kind, 0);
}

void Sunspots::showProgress(ProgressKind kind, int value)
{
	if (kind != ProgressKind::Downloading) removeDownloadProgressBar();
	else if (!downloadProgressBar_)
	{
		downloadProgressBar_ = StelApp::getInstance().addProgressBar();
		downloadProgressBar_->setRange(0, 100);
		downloadProgressBar_->setFormat(QString("%1 %p%").arg(q_("Loading HMI image")));
	}
	progressKind_ = kind;
	progressValue_ = qBound(0, value, 100);
	if (downloadProgressBar_) downloadProgressBar_->setValue(progressValue_);
	emit statusChanged();
}

void Sunspots::removeDownloadProgressBar()
{
	if (!downloadProgressBar_) return;
	StelApp::getInstance().removeProgressBar(downloadProgressBar_);
	downloadProgressBar_ = nullptr;
}

void Sunspots::finishJob()
{
	Q_ASSERT(!workerRunning_ && !downloadReply_ && queryReplies_.isEmpty());
	removeDownloadProgressBar();
	progressKind_ = ProgressKind::None;
	forceImageUpdateForJob_ = false;
	latestForJob_ = false;
	refreshCache();
	// Restart the interval at completion, including when a manual job ran during the interval.
	autoUpdateTimer_.start(60000);
	if (autoDownload_ && hmiPhotosphereEnabled && isFast()) showTooFastWarning();
	emit statusChanged();
}

QString Sunspots::dataPath(const Record& record, const QString& suffix) const
{
	return cacheDir_ + "/" + record.cacheBaseName() + suffix;
}

void Sunspots::ensureUsableImage()
{
	if (isBusy()) return;
	canceledCacheBaseName_.clear();
	if (hmiPhotosphereEnabled && !imageNeedsUpdate(displayedRecord_))
	{
		autoUpdateTimer_.start(60000);
		return;
	}
	const Record cached = usableCachedRecord();
	if (cached.time.isValid())
	{
		if (hmiPhotosphereEnabled && cached.cacheBaseName() != displayedRecord_.cacheBaseName())
		{
			beginJob(ProgressKind::Loading);
			useRecord(cached);
		}
		else autoUpdateTimer_.start(60000);
		return;
	}
	queryAvailability();
}

void Sunspots::downloadForSimTime()
{
	if (isBusy()) return;
	canceledCacheBaseName_.clear();
	queryAvailability(true, true);
}

void Sunspots::downloadLatest()
{
	if (isBusy()) return;
	setAutoDownload(false);
	canceledCacheBaseName_.clear();
	queryAvailability(true, true, true);
}

void Sunspots::displayCachedImage(const QString& cacheBaseName)
{
	if (isBusy()) return;
	refreshCache();
	for (const Record& cached : cachedRecords_)
	{
		if (cached.cacheBaseName() != cacheBaseName) continue;
		const Record record = cached;
		setAutoDownload(false);
		canceledCacheBaseName_.clear();
		beginJob(ProgressKind::Loading);
		if (!hmiPhotosphereEnabled) setHmiPhotosphereEnabled(true);
		useRecord(record);
		return;
	}
}

void Sunspots::cancelCurrentOperation()
{
	if (!isBusy() || isCanceling()) return;
	setAutoDownload(false);
	canceledCacheBaseName_ = jobRecord_.time.isValid() ? jobRecord_.cacheBaseName() : QString();
	qCInfo(sunspots).noquote() << "Canceled operation for" << canceledCacheBaseName_;
	pendingQueries_ = 0;
	for (QNetworkReply* reply : queryReplies_)
	{
		reply->disconnect(this);
		reply->abort();
		reply->deleteLater();
	}
	queryReplies_.clear();
	cancellationToken_->store(true);
	if (downloadReply_)
	{
		downloadReply_->disconnect(this);
		downloadReply_->abort();
		downloadReply_->deleteLater();
		downloadReply_.clear();
	}
	if (downloadFile_) downloadFile_->cancelWriting();
	downloadFile_.reset();
	showError({});
	// Keep the job active until the worker acknowledges cancellation, including during PNG I/O.
	if (workerRunning_) showProgress(ProgressKind::Canceling, 50);
	else finishJob();
}

void Sunspots::queryAvailability(bool userRequested, bool forceImageUpdate, bool latest)
{
	if (isBusy() || cacheDir_.isEmpty()) return;
	if (!userRequested && autoDownload_ && hmiPhotosphereEnabled && isFast())
	{
		showTooFastWarning();
		return;
	}
	beginJob(ProgressKind::Fetching);
	forceImageUpdateForJob_ = forceImageUpdate;
	latestForJob_ = latest;
	if (forceImageUpdate && !hmiPhotosphereEnabled) setHmiPhotosphereEnabled(true);
	totalQueries_ = pendingQueries_ = latest ? 2 : 4;
	for (auto& series : queryResults_)
		for (auto& side : series) side.clear();
	const QString time = observationTimeString(StelApp::getInstance().getCore()->getJDE());
	requestedTime_ = parseObservationTime(time);
	requestRecord(finalSeries, true, time);
	requestRecord(nrtSeries, true, time);
	if (!latest)
	{
		requestRecord(finalSeries, false, time);
		requestRecord(nrtSeries, false, time);
	}
}

void Sunspots::requestRecord(const QString& series, bool before, const QString& time)
{
	const QString first = "2010.01.01_00:00:00_TAI";
	const QString last = "2100.01.01_00:00:00_TAI";
	const QDateTime boundary = parseObservationTime(time);
	const QString overlap = boundary.addSecs(before ? 3600 : -3600).toString("yyyy.MM.dd_HH:mm:ss") + "_TAI";
	const QString recordSet = latestForJob_ ? series :
	                          series + "[" + (before ? first + "-" + overlap : overlap + "-" + last) + "]";
	QUrl url("http://jsoc.stanford.edu/cgi-bin/ajax/jsoc_info");
	QUrlQuery query;
	query.addQueryItem("op", "rs_list");
	query.addQueryItem("ds", recordSet);
	query.addQueryItem("n", before ? "-10" : "10");
	query.addQueryItem("key", "T_REC,T_OBS,QUALITY,CRLT_OBS,CRLN_OBS,CROTA2,DSUN_OBS,RSUN_OBS,RSUN_REF,"
	                            "CRPIX1,CRPIX2,CRVAL1,CRVAL2,CDELT1,CDELT2,CTYPE1,CTYPE2,CUNIT1,CUNIT2");
	query.addQueryItem("seg", "continuum");
	url.setQuery(query);
	qCDebug(sunspots).noquote() << "JSOC query URL:" << url.toString(QUrl::FullyEncoded);
	QNetworkRequest request(url);
	request.setRawHeader("User-Agent", StelUtils::getUserAgentString().toUtf8());
	QNetworkReply* reply = StelApp::getInstance().getNetworkAccessManager()->get(request);
	queryReplies_.insert(reply);
	QTimer::singleShot(30000, reply, [reply]{ if (reply->isRunning()) reply->abort(); });
	connect(reply, &QNetworkReply::finished, this, [this, reply, series, before, recordSet]
	{
		const QByteArray data = reply->readAll();
		const auto networkError = reply->error();
		const QString networkMessage = reply->errorString();
		const int httpStatus = reply->attribute(QNetworkRequest::HttpStatusCodeAttribute).toInt();
		queryReplies_.remove(reply);
		reply->deleteLater();
		QString error;
		if (networkError != QNetworkReply::NoError) error = q_("JSOC query failed: %1").arg(networkMessage);
		const QList<Record> records = error.isEmpty() ? parseRecords(data, series, &error) : QList<Record>();
		if (!error.isEmpty())
		{
			qCWarning(sunspots).noquote() << "JSOC query error for" << recordSet << "(HTTP"
			                              << httpStatus << "):" << error << "Response:"
			                              << QString::fromUtf8(data.left(256)).simplified();
		}
		qCInfo(sunspots).noquote() << "JSOC" << (latestForJob_ ? "latest" : before ? "before" : "after")
		                           << "query for" << series
		                           << "returned" << records.size() << "usable observations; T_OBS range:"
		                           << observationRange(records);
		queryResults_[series == finalSeries ? 0 : 1][before ? 0 : 1] = records;
		if (!error.isEmpty() && queryError_.isEmpty()) queryError_ = error;
		--pendingQueries_;
		showProgress(ProgressKind::Fetching, 100 * (totalQueries_ - pendingQueries_) / totalQueries_);
		if (pendingQueries_ == 0) acceptAvailability();
	});
}

void Sunspots::acceptAvailability()
{
	const QDateTime wanted = requestedTime_;
	auto closest = [&](int series, bool before) -> Record
	{
		Record best;
		for (const auto& side : queryResults_[series])
			for (const Record& record : side)
			{
				if (before ? record.time > wanted : record.time < wanted) continue;
				if (!best.time.isValid() || (before ? record.time > best.time : record.time < best.time)) best = record;
			}
		return best;
	};
	const Record finalBefore = closest(0, true), finalAfter = closest(0, false);
	const Record nrtBefore = closest(1, true), nrtAfter = closest(1, false);
	qCInfo(sunspots).noquote() << "Available T_OBS endpoints for" << wanted.toString("yyyy.MM.dd_HH:mm:ss") + "_TAI: final"
	                           << observationOrNone(finalBefore) << ".." << observationOrNone(finalAfter)
	                           << "; NRT" << observationOrNone(nrtBefore) << ".." << observationOrNone(nrtAfter);
	const int series = finalBefore.time.isValid() && finalAfter.time.isValid() ? 0 : 1;
	const Record earlier = series == 0 ? finalBefore : nrtBefore;
	const Record later = series == 0 ? finalAfter : nrtAfter;
	Record selected;
	if (latestForJob_)
	{
		for (const auto& seriesResults : queryResults_)
			for (const auto& side : seriesResults)
				for (const Record& record : side)
					if (!selected.time.isValid() || record.time > selected.time ||
					    (record.time == selected.time && record.series == finalSeries)) selected = record;
	}
	else if (earlier.time.isValid() && later.time.isValid())
		selected = earlier.time.msecsTo(wanted) <= wanted.msecsTo(later.time) ? earlier : later;
	else selected = earlier.time.isValid() ? earlier : later;
	if (!selected.time.isValid())
	{
		emit statusChanged();
		qCWarning(sunspots).noquote() << "No usable HMI observation selected for"
		                              << wanted.toString("yyyy.MM.dd_HH:mm:ss") + "_TAI";
		showError(queryError_.isEmpty() ? q_("No usable HMI continuum observation is available for this time.")
		                                : queryError_);
		finishJob();
		return;
	}
	qCInfo(sunspots).noquote() << "Selected" << selected.series << selected.observation;
	emit statusChanged();
	if (!tooFastWarning_) showError(queryError_);
	if (!forceImageUpdateForJob_ && autoDownload_ && prepareTextureForJob_ && isFast())
	{
		showTooFastWarning();
		finishJob();
	}
	else if (prepareTextureForJob_ && (forceImageUpdateForJob_ || imageNeedsUpdate(displayedRecord_))) useRecord(selected);
	else finishJob();
}

void Sunspots::useRecord(const Record& record)
{
	if (record.cacheBaseName() == canceledCacheBaseName_ ||
	    record.cacheBaseName() == displayedRecord_.cacheBaseName())
	{
		finishJob();
		return;
	}
	canceledCacheBaseName_.clear();
	jobRecord_ = record;
	if (QFileInfo::exists(dataPath(record, ".png")))
	{
		startImageLoad(record);
		return;
	}
	if (QFileInfo::exists(dataPath(record, ".fits")))
	{
		startTexturePreparation(record);
		return;
	}
	startDownload(record);
}

void Sunspots::startDownload(const Record& record)
{
	Q_ASSERT(isBusy() && !workerRunning_ && !downloadReply_ && queryReplies_.isEmpty());
	if (!record.segment.isValid())
	{
		showError(q_("HMI observation has no downloadable continuum segment."));
		finishJob();
		return;
	}
	downloadFile_.reset(new QSaveFile(dataPath(record, ".fits")));
	if (!downloadFile_->open(QIODevice::WriteOnly))
	{
		showError(q_("Cannot create FITS cache file: %1").arg(downloadFile_->errorString()));
		downloadFile_.reset();
		finishJob();
		return;
	}
	downloadedBytes_ = 0;
	downloadReceivedBytes_ = 0;
	downloadTotalBytes_ = 0;
	downloadWriteFailed_ = false;
	downloadTimedOut_ = false;
	showProgress(ProgressKind::Downloading, 0);
	QNetworkRequest request(record.segment);
	request.setRawHeader("User-Agent", StelUtils::getUserAgentString().toUtf8());
	qCDebug(sunspots).noquote() << "FITS download URL:" << request.url().toString(QUrl::FullyEncoded);
	QNetworkReply* reply = StelApp::getInstance().getNetworkAccessManager()->get(request);
	downloadReply_ = reply;
	QTimer* idleTimer = new QTimer(reply);
	idleTimer->setSingleShot(true);
	idleTimer->setInterval(120000);
	connect(idleTimer, &QTimer::timeout, this, [this, reply]
	{
		if (reply->isRunning())
		{
			downloadTimedOut_ = true;
			reply->abort();
		}
	});
	connect(reply, &QNetworkReply::readyRead, this, [this, reply, idleTimer]
	{
		const QByteArray chunk = reply->readAll();
		if (downloadFile_->write(chunk) != chunk.size())
		{
			downloadWriteFailed_ = true;
			showError(q_("Cannot write FITS cache file: %1").arg(downloadFile_->errorString()));
			reply->abort();
			return;
		}
		downloadedBytes_ += chunk.size();
		if (!chunk.isEmpty()) idleTimer->start();
	});
	connect(reply, &QNetworkReply::downloadProgress, this, [this, reply, idleTimer](qint64 received, qint64 total)
	{
		if (received > downloadReceivedBytes_) idleTimer->start();
		downloadReceivedBytes_ = qMax(downloadReceivedBytes_, received);
		if (total > 0) downloadTotalBytes_ = total;
		if (downloadTotalBytes_ > 0)
			showProgress(ProgressKind::Downloading, 100 * downloadReceivedBytes_ / downloadTotalBytes_);
	});
	connect(reply, &QNetworkReply::finished, this, [this, reply, record]
	{
		if (downloadFile_)
		{
			const QByteArray rest = reply->readAll();
			if (!rest.isEmpty())
			{
				if (downloadFile_->write(rest) != rest.size())
				{
					downloadWriteFailed_ = true;
					showError(q_("Cannot write FITS cache file: %1").arg(downloadFile_->errorString()));
				}
				downloadedBytes_ += rest.size();
			}
		}
		const bool received = reply->error() == QNetworkReply::NoError && !downloadWriteFailed_ &&
		                      downloadedBytes_ > 0 && downloadFile_;
		const bool good = received && downloadFile_->commit();
		if (received && !good)
			showError(q_("Cannot save FITS cache file: %1").arg(downloadFile_->errorString()));
		if (!good && downloadFile_) downloadFile_->cancelWriting();
		if (!good && error_.isEmpty())
			showError(downloadTimedOut_ ? q_("FITS download timed out.") :
			          q_("FITS download failed: %1").arg(reply->errorString()));
		reply->deleteLater();
		downloadReply_.clear();
		downloadFile_.reset();
		if (good)
		{
			if (saveMetadata(record, record.metadata))
			{
				refreshCache();
				startTexturePreparation(record);
			}
			else finishJob();
		}
		else finishJob();
	});
	idleTimer->start();
}

void Sunspots::startTexturePreparation(const Record& record)
{
	Q_ASSERT(isBusy() && !workerRunning_ && !downloadReply_ && queryReplies_.isEmpty());
	showProgress(ProgressKind::PreparingTexture, 0);
	const auto sun = GETSTELMODULE(SolarSystem)->getSun();
	const double jde = StelUtils::taiToTt(StelUtils::qDateTimeToJd(record.time));
	const double rotation = sun->getSiderealTime(jde, jde) + 90.;
	const Mat4d modelToEcliptic = sun->getRotEquatorialToVsop87() * Mat4d::zrotation(rotation * M_PI / 180.);
	const Vec3d observer = StelApp::getInstance().getCore()->getObserverHeliocentricEclipticPos();
	const QString fitsPath = dataPath(record, ".fits");
	const QString pngPath = dataPath(record, ".png");
	QMap<QString, QString> metadata = readMetadata(record);
	for (auto it = record.metadata.cbegin(); it != record.metadata.cend(); ++it)
		if (!metadata.contains(it.key())) metadata.insert(it.key(), it.value());
	const auto token = cancellationToken_;
	workerRunning_ = true;
	workers_.start(QRunnable::create([this, record, fitsPath, pngPath, metadata, modelToEcliptic, observer, token]
	{
		const auto result = SunspotsImage::prepareTexture(fitsPath, pngPath, metadata, modelToEcliptic, observer,
			[this, token](SunspotsImage::PreparationStage stage, int value)
			{
				QMetaObject::invokeMethod(this, [this, stage, value, token]
				{
					if (token->load()) return;
					showProgress(stage == SunspotsImage::PreparationStage::Saving
					             ? ProgressKind::Saving : ProgressKind::PreparingTexture, value);
				}, Qt::QueuedConnection);
			}, [token] { return token->load(); });
		QMetaObject::invokeMethod(this, [this, result, record, token]
		{
			workerRunning_ = false;
			if (result.cancelled || token->load())
			{
				finishJob();
				return;
			}
			if (!result.error.isEmpty())
			{
				showError(result.error);
				finishJob();
			}
			else if (saveMetadata(record, result.metadata))
			{
				refreshCache();
				startImageLoad(record);
			}
			else finishJob();
		}, Qt::QueuedConnection);
	}));
}

void Sunspots::startImageLoad(const Record& record)
{
	Q_ASSERT(isBusy() && !workerRunning_ && !downloadReply_ && queryReplies_.isEmpty());
	const QString path = dataPath(record, ".png");
	const auto token = cancellationToken_;
	showProgress(ProgressKind::Loading, 0);
	workerRunning_ = true;
	workers_.start(QRunnable::create([this, record, path, token]
	{
		QString error;
		const QImage image = SunspotsImage::loadTexture(path, &error, [token] { return token->load(); });
		QMetaObject::invokeMethod(this, [this, record, path, error, image, token]
		{
			workerRunning_ = false;
			if (token->load())
			{
				finishJob();
				return;
			}
			if (!error.isEmpty())
			{
				showError(error);
				QFile::remove(path);
				finishJob();
				return;
			}
			lastReadyRecord_ = record;
			if (hmiPhotosphereEnabled)
			{
				GETSTELMODULE(SolarSystem)->setTextureForPlanet("Sun", image);
				displayedRecord_ = record;
			}
			showError(queryError_);
			finishJob();
		}, Qt::QueuedConnection);
	}));
}

bool Sunspots::saveMetadata(const Record& record, const QMap<QString, QString>& metadata)
{
	QSaveFile csv(dataPath(record, ".csv"));
	if (!csv.open(QIODevice::WriteOnly))
	{
		showError(q_("Cannot write metadata sidecar: %1").arg(csv.errorString()));
		return false;
	}
	for (auto it = metadata.cbegin(); it != metadata.cend(); ++it)
	{
		QString value = it.value();
		value.replace('"', "\"\"");
		if (value.contains(',') || value.contains('"') || value.contains('\n')) value = '"' + value + '"';
		const QByteArray row = (it.key() + ',' + value + '\n').toUtf8();
		if (csv.write(row) != row.size())
		{
			showError(q_("Cannot write metadata sidecar: %1").arg(csv.errorString()));
			return false;
		}
	}
	if (!csv.commit())
	{
		showError(q_("Cannot save metadata sidecar: %1").arg(csv.errorString()));
		return false;
	}
	return true;
}

QMap<QString, QString> Sunspots::readMetadata(const Record& record) const
{
	QFile csv(dataPath(record, ".csv"));
	if (!csv.open(QIODevice::ReadOnly)) return {};
	QMap<QString, QString> metadata;
	while (!csv.atEnd())
	{
		const QString line = QString::fromUtf8(csv.readLine()).trimmed();
		const int comma = line.indexOf(',');
		if (comma <= 0) continue;
		QString value = line.mid(comma + 1);
		if (value.startsWith('"') && value.endsWith('"'))
		{
			value = value.mid(1, value.size() - 2);
			value.replace("\"\"", "\"");
		}
		metadata.insert(line.left(comma), value);
	}
	return metadata;
}
