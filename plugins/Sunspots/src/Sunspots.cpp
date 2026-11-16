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
#include <QDebug>
#include <QDir>
#include <QFile>
#include <QFileInfo>
#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonObject>
#include <QNetworkReply>
#include <QNetworkRequest>
#include <QRegularExpression>
#include <QSaveFile>
#include <QScopedValueRollback>
#include <QTimer>
#include <QTimeZone>
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
const QString phiSeries = QStringLiteral("solo_L2_phi-fdt-icnt");
const QString phiLowLatencySeries = QStringLiteral("solo_LL02_phi-fdt-icnt");
const QString phiJpeg2000Series = QStringLiteral("solo_LL03_phi-fdt-icnt");

int seriesPriority(const QString& series)
{
	// Prefer science FITS, then preliminary FITS, then JPEG2000 quicklooks on equal times/distances.
	if (series == finalSeries || series == phiSeries) return 0;
	return series == phiJpeg2000Series ? 2 : 1;
}

QTimer* startQueryTimeout(QNetworkReply* reply)
{
	// Archive queries can take longer than thirty seconds even after receiving HTTP headers.
	// Keep a separate timer per reply so concurrent queries report their own timeouts.
	auto timer = new QTimer(reply);
	timer->setSingleShot(true);
	timer->setInterval(120000);
	QObject::connect(timer, &QTimer::timeout, reply, &QNetworkReply::abort);
	timer->start();
	return timer;
}

QDateTime parseUtc(const QString& value)
{
	QString utc = value;
	utc.replace(' ', 'T');
	if (!utc.endsWith('Z')) utc += 'Z';
	return QDateTime::fromString(utc, Qt::ISODateWithMs).toUTC();
}

QDateTime parseObservationTime(const QString& value)
{
	QDateTime time = QDateTime::fromString(value.left(19), "yyyy.MM.dd_HH:mm:ss");
	// T_OBS contains TAI calendar fields, but QDateTime has no TAI time specification.
	// Mark them as UTC solely to make QDateTime comparisons and StelUtils::qDateTimeToJd()
	// use these fields consistently. This does not convert TAI to UTC: the QDateTime
	// represents a fictitious UTC instant and must not be used as a real UTC timestamp.
	time.setTimeZone(QTimeZone("UTC"));
	return time;
}

QString intermediateBaseName(const Sunspots::Record& record)
{
	if (record.instrument == Sunspots::Instrument::PHI)
	{
		const QString id = QUrlQuery(record.segment).queryItemValue("data_item_id");
		return id.isEmpty() ? record.series + "_" + record.startTime.toString("yyyyMMdd'T'HHmmss") : id;
	}
	return record.series + "_" + record.time.toString("yyyyMMdd'T'HHmmss'_TAI'");
}

Sunspots::Record recordFromIntermediateFileName(const QString& fileName)
{
	static const QRegularExpression name(
		R"(^(hmi\.Ic_noLimbDark_720s(?:_nrt)?|solo_(?:L2|LL02|LL03)_phi-fdt-icnt)_(\d{8}T\d{6})(_TAI)?\.(fits|jp2)$)");
	const auto match = name.match(fileName);
	if (!match.hasMatch()) return {};
	Sunspots::Record record;
	record.series = match.captured(1);
	record.instrument = record.series.startsWith("solo_") ? Sunspots::Instrument::PHI : Sunspots::Instrument::HMI;
	const bool hmi = record.instrument == Sunspots::Instrument::HMI;
	if (hmi != !match.captured(3).isEmpty()) return {};
	record.time = QDateTime::fromString(match.captured(2), "yyyyMMdd'T'HHmmss");
	record.time.setTimeZone(QTimeZone("UTC"));
	if (!record.time.isValid()) return {};
	if (hmi) record.observation = record.time.toString("yyyy.MM.dd_HH:mm:ss'_TAI'");
	else record.startTime = record.time;
	return intermediateBaseName(record) + record.sourceSuffix() == fileName ? record : Sunspots::Record{};
}

Sunspots::Record recordFromCacheFileName(const QString& fileName)
{
	const QFileInfo file(fileName);
	if (file.suffix() != "fits" && file.suffix() != "png" && file.suffix() != "csv" && file.suffix() != "jp2") return {};
	const QString baseName = file.completeBaseName();
	static const QRegularExpression phiName(
		R"(^(\d{4}\.\d{2}\.\d{2}_\d{2}-\d{2}-\d{2}\.\d{3}Z)_)"
		R"((solo_(?:L2|LL02|LL03)_phi-fdt-icnt)([+-]\d+\.\d{3})s$)");
	const auto match = phiName.match(baseName);
	if (match.hasMatch() && file.suffix() != "csv")
	{
		Sunspots::Record record;
		record.instrument = Sunspots::Instrument::PHI;
		record.series = match.captured(2);
		record.startTime = QDateTime::fromString(match.captured(1), "yyyy.MM.dd_HH-mm-ss.zzz'Z'");
		record.startTime.setTimeZone(QTimeZone("UTC"));
		bool ok = false;
		QString milliseconds = match.captured(3);
		milliseconds.remove('.');
		const qint64 offset = milliseconds.toLongLong(&ok);
		if (!record.startTime.isValid() || !ok) return {};
		record.time = record.startTime.addMSecs(offset);
		if (!record.time.isValid()) return {};
		record.exactTime = true;
		record.observation = record.time.toString(Qt::ISODateWithMs);
		return record.cacheBaseName() == baseName &&
		       (file.suffix() == "png" || "." + file.suffix() == record.sourceSuffix())
		       ? record : Sunspots::Record{};
	}
	if (baseName.size() <= 20 || baseName.at(19) != '_') return {};
	Sunspots::Record record;
	record.series = baseName.mid(20);
	if ((record.series != finalSeries && record.series != nrtSeries) || file.suffix() == "jp2") return {};
	record.observation = baseName.left(19);
	record.observation.replace('-', ':');
	record.observation += "_TAI";
	record.time = parseObservationTime(record.observation);
	if (!record.time.isValid() || record.time.toString("yyyy.MM.dd_HH-mm-ss") + "_" + record.series != baseName)
		return {};
	return record;
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
	if (instrument == Instrument::PHI)
	{
		const qint64 milliseconds = startTime.msecsTo(time);
		const QString offset = (milliseconds < 0 ? "-" : "+") +
		                       QString::number(std::abs(milliseconds) / 1000) + "." +
		                       QString::number(std::abs(milliseconds) % 1000).rightJustified(3, '0');
		return startTime.toString("yyyy.MM.dd_HH-mm-ss.zzz'Z'") + "_" + series +
		       (exactTime ? offset + "s" : "_pending");
	}
	QString timePart = observation.left(19);
	timePart.replace(':', '-');
	return timePart + "_" + series;
}

double Sunspots::Record::julianDayTT() const
{
	const double jd = StelUtils::qDateTimeToJd(time);
	return instrument == Instrument::HMI ? StelUtils::taiToTt(jd) :
	       jd + StelApp::getInstance().getCore()->computeDeltaT(jd) * StelCore::JD_SECOND;
}

QString Sunspots::Record::instrumentName() const
{
	return instrument == Instrument::PHI ? QStringLiteral("PHI") : QStringLiteral("HMI");
}

QString Sunspots::Record::sourceName() const
{
	if (series == nrtSeries) return q_("HMI (NRT)");
	if (series == phiLowLatencySeries) return q_("PHI (LL)");
	if (series == phiJpeg2000Series) return q_("PHI (JP2)");
	return instrumentName();
}

QString Sunspots::Record::sourceSuffix() const
{
	return series == phiJpeg2000Series ? QStringLiteral(".jp2") : QStringLiteral(".fits");
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
	info.description = N_("Replaces the Sun texture with an HMI or PHI image of the photosphere");
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
	if (textureEnabled)
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

	addAction("actionShow_Sunspots", N_("Sunspots"), N_("Show the observed photosphere"), "textureEnabled");
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
		if (isCanceling()) toolbarButton_->setChecked(textureEnabled);
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
	if (textureEnabled)
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
	textureEnabled = conf_->value("Sunspots/show_texture", false).toBool();
	autoDownload_ = conf_->value("Sunspots/auto_download", false).toBool();
	instrument_ = conf_->value("Sunspots/instrument", "HMI").toString() == "PHI" ? Instrument::PHI : Instrument::HMI;
	continuousUseCacheOnly_ = conf_->value("Sunspots/continuous_use_cache_only", false).toBool();
	deleteSourceAfterMakingTexture_ = conf_->value("Sunspots/delete_source_after_making_texture", false).toBool();
	maxImageAgeMinutes_ = qBound(6, conf_->value("Sunspots/max_image_age_minutes", 180).toInt(), 9999);
}

void Sunspots::setDialogVisible(const bool enable)
{
	if(dialogVisible == enable) return;

	dialogVisible = enable;
	dialog_->setVisible(enable);
	emit dialogToggled(enable);
}

void Sunspots::setTextureEnabled(const bool enable)
{
	if (isCanceling() || textureEnabled == enable) return;

	textureEnabled = enable;
	conf_->setValue("Sunspots/show_texture", enable);
	emit textureToggled(enable);
	if (enable)
	{
		if (isBusy()) prepareTextureForJob_ = true;
		else ensureUsableImage();
	}
	else
	{
		cacheEmptyError_ = false;
		displayedRecord_ = {};
		displayedImageDeleted_ = false;
		GETSTELMODULE(SolarSystem)->resetTextures("Sun");
		emit statusChanged();
		if (!isBusy()) showError({});
	}
}

void Sunspots::setAutoDownload(bool enable)
{
	// Changes affect subsequent jobs; Cancel can still disable future polling immediately.
	if (autoDownload_ == enable) return;
	autoDownload_ = enable;
	conf_->setValue("Sunspots/auto_download", enable);
	emit autoDownloadChanged(enable);
	if (!isBusy())
	{
		if (enable && textureEnabled) pollTime();
		else if (!enable && tooFastWarning_) showError({});
	}
}

void Sunspots::setInstrument(Instrument instrument)
{
	if (isBusy() || cacheClearPending_ || instrument_ == instrument ||
	    (instrument != Instrument::HMI && instrument != Instrument::PHI)) return;
	instrument_ = instrument;
	conf_->setValue("Sunspots/instrument", instrument == Instrument::PHI ? "PHI" : "HMI");
	cacheEmptyError_ = false;
	canceledCacheBaseName_.clear();
	emit instrumentChanged(instrument);
	if (textureEnabled) ensureUsableImage();
	emit statusChanged();
}

void Sunspots::setContinuousUseCacheOnly(bool enable)
{
	if (continuousUseCacheOnly_ == enable) return;
	continuousUseCacheOnly_ = enable;
	conf_->setValue("Sunspots/continuous_use_cache_only", enable);
	if (!enable) cacheEmptyError_ = false;
	emit continuousUseCacheOnlyChanged(enable);
	emit statusChanged();
}

void Sunspots::setDeleteSourceAfterMakingTexture(bool enable)
{
	if (deleteSourceAfterMakingTexture_ == enable) return;
	deleteSourceAfterMakingTexture_ = enable;
	conf_->setValue("Sunspots/delete_source_after_making_texture", enable);
	emit deleteSourceAfterMakingTextureChanged(enable);
}

void Sunspots::clearCache()
{
	if (cacheClearPending_ || cacheDir_.isEmpty()) return;
	QScopedValueRollback<bool> clearing(cacheClearPending_, true);
	cacheClearProtectedBaseNames_.clear();
	if (isBusy() && jobRecord_.time.isValid()) protectJobFiles(jobRecord_);
	// The modal confirmation runs an event loop. Existing jobs may select an image or finish
	// there; useRecord() adds their identities to the protection set until clearing completes.
	if (dialog_->confirmCacheRemoval(textureEnabled && displayedRecord_.time.isValid(), isBusy()))
	{
		if (textureEnabled && displayedRecord_.time.isValid())
			cacheClearProtectedBaseNames_.insert(displayedRecord_.cacheBaseName());
		removeCachedFiles(cacheClearProtectedBaseNames_);
	}
	cacheClearProtectedBaseNames_.clear();
	autoUpdateTimer_.start(60000);
	// A job can finish while the confirmation's nested event loop is running.
	QTimer::singleShot(0, this, &Sunspots::resumeSelectionAfterDeletion);
}

bool Sunspots::canDeleteCachedImage(const QString& cacheBaseName, QString* disabledReason) const
{
	QString reason;
	const bool active = isBusy() && jobRecord_.time.isValid() && jobRecord_.cacheBaseName() == cacheBaseName;
	if (cacheClearPending_) reason = q_("A cache deletion operation is in progress.");
	else if (active && isCanceling()) reason = q_("Cancellation of this image's operation is in progress.");
	else if (cacheDir_.isEmpty() || (!active && std::none_of(cachedRecords_.cbegin(), cachedRecords_.cend(),
	    [&cacheBaseName](const Record& record) { return record.cacheBaseName() == cacheBaseName; })))
		reason = q_("This image is no longer available in the cache.");
	if (disabledReason) *disabledReason = reason;
	return reason.isEmpty();
}

void Sunspots::deleteCachedImage(QString cacheBaseName)
{
	if (!canDeleteCachedImage(cacheBaseName)) return;
	const bool active = isBusy() && jobRecord_.time.isValid() && jobRecord_.cacheBaseName() == cacheBaseName;
	const Record record = active ? jobRecord_ : recordFromCacheFileName(cacheBaseName + ".png");
	if (!record.time.isValid() || record.cacheBaseName() != cacheBaseName) return;
	bool confirmed;
	{
		// Keep automatic selection from starting another job during the modal event loop.
		// The existing job may finish or rename its PHI file while confirmation is open.
		QScopedValueRollback<bool> confirming(cacheClearPending_, true);
		confirmed = dialog_->confirmImageRemoval(active);
	}
	cacheClearProtectedBaseNames_.clear();
	autoUpdateTimer_.start(60000);
	if (!confirmed)
	{
		resumeSelectionAfterDeletion();
		return;
	}
	const QString currentBaseName = active ? jobRecord_.cacheBaseName() : cacheBaseName;
	if (isBusy() && jobRecord_.time.isValid() && jobRecord_.cacheBaseName() == currentBaseName)
	{
		jobCacheBaseNameToDelete_ = cacheBaseName;
		cancelCurrentOperation();
		return;
	}
	QStringList errors;
	for (const QString& baseName : QSet<QString>{cacheBaseName, currentBaseName})
		errors.append(removeCachedImageFiles(baseName));
	refreshCache();
	if (textureEnabled && continuousUseCacheOnly_ && !hasCachedImagesForInstrument()) cacheEmptyError_ = true;
	showError(errors.join('\n'));
	emit statusChanged();
	resumeSelectionAfterDeletion();
}

QStringList Sunspots::removeCachedImageFiles(const QString& cacheBaseName)
{
	if (displayedRecord_.time.isValid() && displayedRecord_.cacheBaseName() == cacheBaseName)
	{
		displayedImageDeleted_ = true;
		selectionAfterDeletionPending_ = true;
	}
	QStringList errors;
	QStringList fileNames;
	for (const QString& suffix : {QStringLiteral(".png"), QStringLiteral(".fits"),
	                              QStringLiteral(".jp2"), QStringLiteral(".csv")})
		fileNames.append(cacheBaseName + suffix);
	const Record record = recordFromCacheFileName(cacheBaseName + ".png");
	if (record.time.isValid()) fileNames.append(intermediateBaseName(record) + record.sourceSuffix());
	// The active PHI row can still have its provisional identity while its header is being read.
	if (jobRecord_.time.isValid() && jobRecord_.cacheBaseName() == cacheBaseName)
		fileNames.append(intermediateBaseName(jobRecord_) + jobRecord_.sourceSuffix());
	fileNames.removeDuplicates();
	for (const QString& fileName : fileNames)
	{
		QFile file(cacheDir_ + "/" + fileName);
		if (file.exists() && !file.remove())
			errors.append(q_("Cannot remove cached file %1: %2")
			              .arg(QFileInfo(file).fileName(), file.errorString()));
	}
	if (lastReadyRecord_.cacheBaseName() == cacheBaseName &&
	    !QFileInfo::exists(cacheDir_ + "/" + cacheBaseName + ".png") &&
	    !QFileInfo::exists(cacheDir_ + "/" + cacheBaseName + record.sourceSuffix()))
		lastReadyRecord_ = {};
	return errors;
}

void Sunspots::resumeSelectionAfterDeletion()
{
	if (!selectionAfterDeletionPending_ || isBusy() || cacheClearPending_) return;
	selectionAfterDeletionPending_ = false;
	if (textureEnabled) ensureUsableImage();
}

void Sunspots::removeCachedFiles(const QSet<QString>& protectedBaseNames)
{
	QStringList errors;
	const QDir dir(cacheDir_);
	for (const QFileInfo& file : dir.entryInfoList(QDir::Files | QDir::Hidden | QDir::System))
	{
		const Record record = recordFromCacheFileName(file.fileName());
		if ((!record.time.isValid() && !recordFromIntermediateFileName(file.fileName()).time.isValid()) ||
		    protectedBaseNames.contains(file.completeBaseName())) continue;
		QFile cachedFile(file.absoluteFilePath());
		if (!cachedFile.remove())
		{
			const QString error = q_("Cannot remove cached file %1: %2")
			                      .arg(file.fileName(), cachedFile.errorString());
			errors.append(error);
		}
	}
	if (lastReadyRecord_.time.isValid() && !QFileInfo::exists(dataPath(lastReadyRecord_, ".png")) &&
	    !QFileInfo::exists(dataPath(lastReadyRecord_, lastReadyRecord_.sourceSuffix()))) lastReadyRecord_ = {};
	refreshCache();
	if (textureEnabled && continuousUseCacheOnly_ && !hasCachedImagesForInstrument()) cacheEmptyError_ = true;
	showError(errors.join('\n'));
	emit statusChanged();
}

void Sunspots::removeImageSource(const Record& record)
{
	// Read the setting here, rather than capturing it at job start: it can be changed during preparation.
	if (!deleteSourceAfterMakingTexture_) return;
	QFile source(dataPath(record, record.sourceSuffix()));
	if (source.exists() && !source.remove())
	{
		textureCleanupError_ = q_("Cannot remove source image %1: %2")
		                       .arg(source.fileName(), source.errorString());
		showError(textureCleanupError_);
	}
}

void Sunspots::setMaxImageAgeMinutes(int minutes)
{
	minutes = qBound(6, minutes, 9999);
	if (maxImageAgeMinutes_ == minutes) return;
	maxImageAgeMinutes_ = minutes;
	conf_->setValue("Sunspots/max_image_age_minutes", minutes);
	emit maxImageAgeMinutesChanged(minutes);
	if (autoDownload_ && textureEnabled) pollTime();
}

QString Sunspots::observationTimeText(const Record& record) const
{
	if (!record.time.isValid()) return q_("N/A");
	const auto core = StelApp::getInstance().getCore();
	const double tt = record.julianDayTT();
	// PHI already provides a UTC instant. For HMI, follow Stellarium's TT/UT convention
	// rather than interpreting the QDateTime carrying TAI fields as genuine UTC.
	// Evaluate DeltaT and the zone offset at the observation's date, not simulation time.
	const double ut = record.instrument == Instrument::PHI ? StelUtils::qDateTimeToJd(record.time) :
	                  tt - core->computeDeltaT(tt) * StelCore::JD_SECOND;
	const double offset = core->getUTCOffset(ut);
	const auto& locale = StelApp::getInstance().getLocaleMgr();
	return locale.getPrintableDateLocal(ut, offset) + " " + locale.getPrintableTimeLocal(ut, offset);
}

QString Sunspots::displayedImageTitle() const
{
	const Instrument instrument = displayedRecord_.time.isValid() ? displayedRecord_.instrument : instrument_;
	return instrument == Instrument::PHI ? q_("Currently displayed PHI image") :
	                                      q_("Currently displayed HMI image");
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
	const double age = (referenceTT - displayedRecord_.julianDayTT()) * 1440.;
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
	if (!record.time.isValid() || record.instrument != instrument_) return true;
	return std::abs(StelApp::getInstance().getCore()->getJDE() - record.julianDayTT()) * 1440. >= maxImageAgeMinutes_;
}

bool Sunspots::hasCachedImagesForInstrument() const
{
	return std::any_of(cachedRecords_.cbegin(), cachedRecords_.cend(), [this](const Record& record)
	                  { return record.instrument == instrument_; });
}

Sunspots::Record Sunspots::usableCachedRecord(bool allowStale) const
{
	if (cacheDir_.isEmpty()) return {};
	if (!allowStale && !imageNeedsUpdate(lastReadyRecord_) && QFileInfo::exists(dataPath(lastReadyRecord_, ".png")))
		return lastReadyRecord_;

	Record best;
	const double simulatedTT = StelApp::getInstance().getCore()->getJDE();
	const QDir dir(cacheDir_);
	for (const QFileInfo& file : dir.entryInfoList(
	             {QStringLiteral("*.png"), QStringLiteral("*.fits"), QStringLiteral("*.jp2")}, QDir::Files))
	{
		if (file.size() == 0) continue;
		const Record record = recordFromCacheFileName(file.fileName());
		if (!record.time.isValid() || record.instrument != instrument_ ||
		    (!allowStale && imageNeedsUpdate(record))) continue;
		const double distance = std::abs(record.julianDayTT() - simulatedTT);
		const double bestDistance = best.time.isValid()
		                            ? std::abs(best.julianDayTT() - simulatedTT) : 0.;
		if (!best.time.isValid() || distance < bestDistance ||
		    (distance == bestDistance && ((seriesPriority(record.series) < seriesPriority(best.series)) ||
		     (record.series == best.series && record.time < best.time))))
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
		                                      QStringLiteral("*.csv"), QStringLiteral("*.jp2")}, QDir::Files);
		for (const QFileInfo& file : files)
		{
			cacheSizeBytes_ += file.size();
			if (file.size() == 0 || file.suffix() == "csv") continue;
			const QString baseName = file.completeBaseName();
			Record record = recordFromCacheFileName(file.fileName());
			if (!record.time.isValid()) continue;
			const bool png = file.suffix() == "png";
			auto existing = records.find(baseName);
			if (existing != records.end())
			{
				existing->hasCachedPng |= png;
				continue;
			}
			record.hasCachedPng = png;
			records.insert(baseName, record);
		}
		cachedRecords_ = records.values();
		std::sort(cachedRecords_.begin(), cachedRecords_.end(), [](const Record& a, const Record& b)
		{
			if (a.julianDayTT() != b.julianDayTT()) return a.julianDayTT() > b.julianDayTT();
			return a.series < b.series;
		});
	}
	if (hasCachedImagesForInstrument()) cacheEmptyError_ = false;
	emit cacheChanged();
}

void Sunspots::pollTime()
{
	if (!textureEnabled || !autoDownload_ || isBusy() || cacheClearPending_) return;
	if (isFast())
	{
		showTooFastWarning();
		return;
	}
	if (tooFastWarning_) showError({});
	if (autoUpdateTimer_.isActive()) return;
	if (!displayedImageDeleted_ && !imageNeedsUpdate(displayedRecord_))
	{
		autoUpdateTimer_.start(60000);
		return;
	}
	ensureUsableImage();
}

void Sunspots::showError(const QString& error)
{
	if (!error.isEmpty()) qCWarning(sunspots).noquote() << error;
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
	qCWarning(sunspots).noquote() << errorText();
	emit statusChanged();
}

void Sunspots::beginJob(ProgressKind kind)
{
	Q_ASSERT(!isBusy());
	cancellationToken_ = std::make_shared<std::atomic_bool>(false);
	prepareTextureForJob_ = textureEnabled;
	autoDownloadForJob_ = autoDownload_;
	cacheOnlyForJob_ = continuousUseCacheOnly_;
	jobRecord_ = {};
	queryError_.clear();
	textureCleanupError_.clear();
	cacheEmptyError_ = false;
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
		downloadProgressBar_->setFormat(QString("%1 %p%").arg(jobRecord_.instrument == Instrument::PHI ?
		                               q_("Loading PHI image") : q_("Loading HMI image")));
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
	QStringList deletionErrors;
	if (!jobCacheBaseNameToDelete_.isEmpty())
	{
		// PHI header finalization may have renamed the completed source while cancellation settled.
		const QSet<QString> baseNames{jobCacheBaseNameToDelete_, jobRecord_.cacheBaseName()};
		jobCacheBaseNameToDelete_.clear();
		for (const QString& baseName : baseNames) deletionErrors.append(removeCachedImageFiles(baseName));
	}
	refreshCache();
	if (!deletionErrors.isEmpty()) showError(deletionErrors.join('\n'));
	// Restart the interval at completion, including when a manual job ran during the interval.
	autoUpdateTimer_.start(60000);
	if (autoDownload_ && textureEnabled && isFast()) showTooFastWarning();
	emit statusChanged();
	resumeSelectionAfterDeletion();
}

QString Sunspots::dataPath(const Record& record, const QString& suffix) const
{
	return cacheDir_ + "/" + record.cacheBaseName() + suffix;
}

QString Sunspots::intermediatePath(const Record& record) const
{
	return cacheDir_ + "/" + intermediateBaseName(record) + record.sourceSuffix();
}

void Sunspots::protectJobFiles(const Record& record)
{
	if (!cacheClearPending_) return;
	cacheClearProtectedBaseNames_.insert(record.cacheBaseName());
	cacheClearProtectedBaseNames_.insert(intermediateBaseName(record));
}

void Sunspots::ensureUsableImage()
{
	if (isBusy() || cacheClearPending_ || cacheDir_.isEmpty()) return;
	if (autoDownload_ && textureEnabled && isFast())
	{
		showTooFastWarning();
		return;
	}
	canceledCacheBaseName_.clear();
	if (textureEnabled && !displayedImageDeleted_ && !imageNeedsUpdate(displayedRecord_))
	{
		autoUpdateTimer_.start(60000);
		return;
	}
	refreshCache();
	const Record cached = usableCachedRecord(continuousUseCacheOnly_);
	if (cached.time.isValid())
	{
		if (textureEnabled && cached.cacheBaseName() != displayedRecord_.cacheBaseName())
		{
			beginJob(ProgressKind::Loading);
			useRecord(cached);
		}
		else autoUpdateTimer_.start(60000);
		return;
	}
	if (continuousUseCacheOnly_)
	{
		cacheEmptyError_ = true;
		showError({});
		emit statusChanged();
		autoUpdateTimer_.start(60000);
		return;
	}
	queryAvailability();
}

void Sunspots::downloadForSimTime()
{
	if (isBusy() || cacheClearPending_) return;
	canceledCacheBaseName_.clear();
	queryAvailability(true, true);
}

void Sunspots::downloadLatest()
{
	if (isBusy() || cacheClearPending_) return;
	setAutoDownload(false);
	canceledCacheBaseName_.clear();
	queryAvailability(true, true, true);
}

bool Sunspots::canDisplayCachedImage(const QString& cacheBaseName, QString* disabledReason) const
{
	QString reason;
	if (isBusy()) reason = q_("An image update is in progress.");
	else if (cacheClearPending_) reason = q_("A cache deletion operation is in progress.");
	else if (!displayedImageDeleted_ && cacheBaseName == displayedRecord_.cacheBaseName())
		reason = q_("This image is already displayed.");
	else if (std::none_of(cachedRecords_.cbegin(), cachedRecords_.cend(), [&cacheBaseName](const Record& record)
	{
		return record.cacheBaseName() == cacheBaseName;
	})) reason = q_("This image is no longer available in the cache.");
	if (disabledReason) *disabledReason = reason;
	return reason.isEmpty();
}

void Sunspots::displayCachedImage(const QString& cacheBaseName)
{
	if (isBusy() || cacheClearPending_) return;
	refreshCache();
	if (!canDisplayCachedImage(cacheBaseName)) return;
	for (const Record& cached : cachedRecords_)
	{
		if (cached.cacheBaseName() != cacheBaseName) continue;
		const Record record = cached;
		setAutoDownload(false);
		// Change source without starting automatic selection before this explicit job.
		instrument_ = record.instrument;
		conf_->setValue("Sunspots/instrument", record.instrumentName());
		emit instrumentChanged(instrument_);
		canceledCacheBaseName_.clear();
		beginJob(ProgressKind::Loading);
		cacheOnlyForJob_ = true;
		if (!textureEnabled) setTextureEnabled(true);
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
	if (isBusy() || cacheClearPending_ || cacheDir_.isEmpty()) return;
	if (!userRequested && autoDownload_ && textureEnabled && isFast())
	{
		showTooFastWarning();
		return;
	}
	refreshCache();
	beginJob(ProgressKind::Fetching);
	// An explicit Download command bypasses the automatic cache-only policy.
	if (userRequested) cacheOnlyForJob_ = false;
	forceImageUpdateForJob_ = forceImageUpdate;
	latestForJob_ = latest;
	requestedTimeTT_ = StelApp::getInstance().getCore()->getJDE();
	if (forceImageUpdate && !textureEnabled) setTextureEnabled(true);
	totalQueries_ = pendingQueries_ = latest ? 2 : 4;
	requestedPhiJpeg2000_ = false;
	for (auto& series : queryResults_)
		for (auto& side : series) side.clear();
	if (instrument_ == Instrument::PHI)
	{
		totalQueries_ = pendingQueries_ = 2;
		requestPhiRecords(false);
		requestPhiRecords(true);
		return;
	}
	const QString time = observationTimeString(requestedTimeTT_);
	requestRecord(finalSeries, true, time);
	requestRecord(nrtSeries, true, time);
	if (!latest)
	{
		requestRecord(finalSeries, false, time);
		requestRecord(nrtSeries, false, time);
	}
}

void Sunspots::requestPhiRecords(bool lowLatency)
{
	const QString series = lowLatency ? phiLowLatencySeries : phiSeries;
	const QString table = lowLatency ? "v_ll_data_item" : "v_sc_data_item";
	const QString level = lowLatency ? "LL02" : "L2";
	QUrl url("https://soar.esac.esa.int/soar-sl-tap/tap/sync");
	QUrlQuery query;
	query.addQueryItem("REQUEST", "doQuery");
	query.addQueryItem("LANG", "ADQL");
	query.addQueryItem("FORMAT", "json");
	query.addQueryItem("MAXREC", latestForJob_ ? "10" : "100000");
	const QString select = latestForJob_ ? "SELECT TOP 10 " : "SELECT ";
	query.addQueryItem("QUERY", select + "begin_time,end_time,filename,data_item_id FROM " + table +
	                          " WHERE descriptor='phi-fdt-icnt' AND level='" + level + "'" +
	                          (latestForJob_ ? " ORDER BY begin_time DESC,filename DESC" : ""));
	url.setQuery(query);
	qCInfo(sunspots).noquote() << "SOAR query URL:" << url.toString(QUrl::FullyEncoded);
	QNetworkRequest request(url);
	request.setRawHeader("User-Agent", StelUtils::getUserAgentString().toUtf8());
	QNetworkReply* reply = StelApp::getInstance().getNetworkAccessManager()->get(request);
	queryReplies_.insert(reply);
	const auto timeout = startQueryTimeout(reply);
	connect(reply, &QNetworkReply::finished, this, [this, reply, series, lowLatency, timeout]
	{
		const bool timedOut = !timeout->isActive();
		timeout->stop();
		const QByteArray data = reply->isOpen() ? reply->readAll() : QByteArray();
		queryReplies_.remove(reply);
		reply->deleteLater();
		QList<Record> records;
		QString error;
		try
		{
			if (timedOut)
				throw std::runtime_error(q_("SOAR query timed out.").toStdString());
			if (reply->error() != QNetworkReply::NoError)
				throw std::runtime_error(q_("SOAR query failed: %1").arg(reply->errorString()).toStdString());
			QJsonParseError parseError;
			const auto document = QJsonDocument::fromJson(data, &parseError);
			if (parseError.error != QJsonParseError::NoError || !document.isObject())
				throw std::runtime_error(q_("Invalid SOAR response: %1").arg(parseError.errorString()).toStdString());
			const auto root = document.object();
			QMap<QString, int> columns;
			const auto metadata = root.value("metadata").toArray();
			for (int i = 0; i < metadata.size(); ++i)
				columns.insert(metadata[i].toObject().value("name").toString(), i);
			for (const char* name : {"begin_time", "end_time", "filename", "data_item_id"})
				if (!columns.contains(name))
					throw std::runtime_error(q_("SOAR response is missing column %1.").arg(name).toStdString());
			if (!root.value("data").isArray())
				throw std::runtime_error(q_("SOAR response is missing its data table.").toStdString());
			QMap<QDateTime, Record> observations;
			QMap<QDateTime, QString> fileNames;
			for (const auto& value : root.value("data").toArray())
			{
				const auto row = value.toArray();
				auto field = [&](const QString& key)
				{
					const int column = columns.value(key);
					return column < row.size() ? row.at(column).toString() : QString();
				};
				Record record;
				record.instrument = Instrument::PHI;
				record.series = series;
				record.startTime = parseUtc(field("begin_time"));
				const auto end = parseUtc(field("end_time"));
				const QString fileName = field("filename"), id = field("data_item_id");
				if (!record.startTime.isValid() || !end.isValid() || end < record.startTime || id.isEmpty() ||
				    !fileName.startsWith(series + "_") || !fileName.endsWith(".fits"))
				{
					qCWarning(sunspots) << "Skipping malformed SOAR row:" << value;
					continue;
				}
				record.time = record.startTime.addMSecs(record.startTime.msecsTo(end) / 2);
				record.observation = record.time.toString(Qt::ISODateWithMs);
				record.segment = QUrl("https://soar.esac.esa.int/soar-sl-tap/data");
				QUrlQuery download;
				download.addQueryItem("retrieval_type", "LAST_PRODUCT");
				download.addQueryItem("product_type", lowLatency ? "LOW_LATENCY" : "SCIENCE");
				download.addQueryItem("data_item_id", id);
				record.segment.setQuery(download);
				// The listing identifies cached files by product and start time, independently of its
				// estimated midpoint. Cached DATE-AVG replaces that estimate when available.
				for (const Record& cached : cachedRecords_)
					if (cached.series == series && cached.startTime == record.startTime)
					{
						record.time = cached.time;
						record.exactTime = true;
						break;
					}
				if (fileName > fileNames.value(record.startTime))
				{
					observations.insert(record.startTime, record);
					fileNames.insert(record.startTime, fileName);
				}
			}
			records = observations.values();
		}
		catch (const std::runtime_error& e)
		{
			error = QString::fromUtf8(e.what());
			qCWarning(sunspots).noquote() << "SOAR query error for" << series << "(HTTP"
			                              << reply->attribute(QNetworkRequest::HttpStatusCodeAttribute).toInt()
			                              << "):" << error << "Response:" << QString::fromUtf8(data.left(256));
		}
		const QString first = records.isEmpty() ? "(none)" : records.first().startTime.toString(Qt::ISODateWithMs);
		const QString last = records.isEmpty() ? "(none)" : records.last().startTime.toString(Qt::ISODateWithMs);
		qCInfo(sunspots).noquote() << "SOAR" << series << "returned" << records.size()
		                           << "usable observations; UTC start-time range:"
		                           << first << ".." << last;
		queryResults_[lowLatency ? 1 : 0][0] = records;
		if (!error.isEmpty() && queryError_.isEmpty()) queryError_ = error;
		--pendingQueries_;
		showProgress(ProgressKind::Fetching, 100 * (totalQueries_ - pendingQueries_) / totalQueries_);
		if (pendingQueries_ == 0) acceptAvailability();
	});
}

void Sunspots::requestPhiJpeg2000(const QString& date, bool findLatest)
{
	QUrl url(QStringLiteral("https://api.swhv.oma.be/hv_docpage/v2/") +
	         (findLatest ? "getDataSources/" : "getJP2Image/"));
	if (!findLatest)
	{
		QUrlQuery query;
		query.addQueryItem("date", date);
		query.addQueryItem("sourceId", "500");
		query.addQueryItem("jpip", "true");
		query.addQueryItem("json", "true");
		url.setQuery(query);
	}
	qCInfo(sunspots).noquote() << "PHI JPEG2000 query URL:" << url.toString(QUrl::FullyEncoded);
	QNetworkRequest request(url);
	request.setRawHeader("User-Agent", StelUtils::getUserAgentString().toUtf8());
	QNetworkReply* reply = StelApp::getInstance().getNetworkAccessManager()->get(request);
	queryReplies_.insert(reply);
	const auto timeout = startQueryTimeout(reply);
	connect(reply, &QNetworkReply::finished, this, [this, reply, findLatest, timeout]
	{
		const bool timedOut = !timeout->isActive();
		timeout->stop();
		const QByteArray data = reply->isOpen() ? reply->readAll() : QByteArray();
		queryReplies_.remove(reply);
		reply->deleteLater();
		QString error;
		bool nextQueryStarted = false;
		try
		{
			if (timedOut)
				throw std::runtime_error(q_("PHI JPEG2000 query timed out.").toStdString());
			if (reply->error() != QNetworkReply::NoError)
				throw std::runtime_error(q_("PHI JPEG2000 query failed: %1")
				                         .arg(reply->errorString()).toStdString());
			QJsonParseError parseError;
			const auto document = QJsonDocument::fromJson(data, &parseError);
			if (parseError.error != QJsonParseError::NoError || !document.isObject())
				throw std::runtime_error(q_("Invalid PHI JPEG2000 response: %1")
				                         .arg(parseError.errorString()).toStdString());
			if (findLatest)
			{
				QJsonObject source;
				std::function<void(const QJsonObject&)> findSource = [&](const QJsonObject& object)
				{
					if (object.value("sourceId").toInt() == 500) { source = object; return; }
					for (auto it = object.begin(); it != object.end(); ++it)
						if (it.value().isObject()) findSource(it.value().toObject());
				};
				findSource(document.object());
				const auto first = parseUtc(source.value("start").toString());
				const auto last = parseUtc(source.value("end").toString());
				if (!first.isValid() || !last.isValid() || last < first)
					throw std::runtime_error(q_("PHI JPEG2000 source has no usable time range.").toStdString());
				qCInfo(sunspots).noquote() << "PHI JPEG2000 UTC time range:"
				                           << first.toString(Qt::ISODateWithMs) << ".."
				                           << last.toString(Qt::ISODateWithMs);
				requestPhiJpeg2000(last.toString(Qt::ISODateWithMs), false);
				nextQueryStarted = true;
			}
			else
			{
				const QUrl uri(document.object().value("uri").toString());
				const QString fileName = QFileInfo(uri.path()).fileName();
				static const QRegularExpression name(R"(^solo_LL03_phi-fdt-icnt_(\d{8}T\d{6})_V[^/]+\.jp2$)");
				const auto match = name.match(fileName);
				if (!uri.isValid() || uri.scheme() != "jpip" || !match.hasMatch())
					throw std::runtime_error(q_("PHI JPEG2000 response has an invalid image filename.").toStdString());
				Record record;
				record.instrument = Instrument::PHI;
				record.series = phiJpeg2000Series;
				record.startTime = QDateTime::fromString(match.captured(1), "yyyyMMdd'T'HHmmss");
				record.startTime.setTimeZone(QTimeZone("UTC"));
				if (!record.startTime.isValid())
					throw std::runtime_error(q_("PHI JPEG2000 response has an invalid observation time.").toStdString());
				record.time = record.startTime;
				for (const Record& cached : cachedRecords_)
					if (cached.series == record.series &&
					    cached.startTime.toSecsSinceEpoch() == record.startTime.toSecsSinceEpoch())
					{
						record.startTime = cached.startTime;
						record.time = cached.time;
						record.exactTime = true;
						break;
					}
				record.observation = record.time.toString(Qt::ISODateWithMs);
				record.segment = QUrl("https://api.swhv.oma.be/hv_docpage/v2/getJP2Image/");
				QUrlQuery download;
				// Freeze the selected observation rather than querying again at the original requested time.
				download.addQueryItem("date", record.startTime.toString("yyyy-MM-dd'T'HH:mm:ss'Z'"));
				download.addQueryItem("sourceId", "500");
				record.segment.setQuery(download);
				queryResults_[2][0].append(record);
				qCInfo(sunspots).noquote() << "PHI JPEG2000 selected archive filename:" << fileName
				                           << "; UTC observation:" << record.observation;
			}
		}
		catch (const std::runtime_error& e)
		{
			error = QString::fromUtf8(e.what());
			qCWarning(sunspots).noquote() << "PHI JPEG2000 query error (HTTP"
			                              << reply->attribute(QNetworkRequest::HttpStatusCodeAttribute).toInt()
			                              << "):" << error << "Response:" << QString::fromUtf8(data.left(256));
		}
		if (!error.isEmpty() && queryError_.isEmpty()) queryError_ = error;
		// A failed range lookup also finishes the image lookup that depended on it.
		pendingQueries_ -= findLatest && !nextQueryStarted ? 2 : 1;
		showProgress(ProgressKind::Fetching, 100 * (totalQueries_ - pendingQueries_) / totalQueries_);
		if (pendingQueries_ == 0) acceptAvailability();
	});
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
	qCInfo(sunspots).noquote() << "JSOC query URL:" << url.toString(QUrl::FullyEncoded);
	QNetworkRequest request(url);
	request.setRawHeader("User-Agent", StelUtils::getUserAgentString().toUtf8());
	QNetworkReply* reply = StelApp::getInstance().getNetworkAccessManager()->get(request);
	queryReplies_.insert(reply);
	const auto timeout = startQueryTimeout(reply);
	connect(reply, &QNetworkReply::finished, this, [this, reply, series, before, recordSet, timeout]
	{
		const bool timedOut = !timeout->isActive();
		timeout->stop();
		const QByteArray data = reply->isOpen() ? reply->readAll() : QByteArray();
		const auto networkError = reply->error();
		const QString networkMessage = reply->errorString();
		const int httpStatus = reply->attribute(QNetworkRequest::HttpStatusCodeAttribute).toInt();
		queryReplies_.remove(reply);
		reply->deleteLater();
		QString error;
		if (timedOut) error = q_("JSOC query timed out.");
		else if (networkError != QNetworkReply::NoError) error = q_("JSOC query failed: %1").arg(networkMessage);
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
	const double wanted = requestedTimeTT_;
	const bool phi = instrument_ == Instrument::PHI;
	const double wantedUT = wanted - StelApp::getInstance().getCore()->computeDeltaT(wanted) * StelCore::JD_SECOND;
	const QString wantedText = phi ? StelUtils::jdToQDateTime(wantedUT, Qt::UTC).toString(Qt::ISODateWithMs) :
	                                 observationTimeString(wanted);
	auto closest = [&](int series, bool before) -> Record
	{
		Record best;
		for (const auto& side : queryResults_[series])
			for (const Record& record : side)
			{
				if (before ? record.julianDayTT() > wanted : record.julianDayTT() < wanted) continue;
				if (!best.time.isValid() || (before ? record.time > best.time : record.time < best.time)) best = record;
			}
		return best;
	};
	const Record finalBefore = closest(0, true), finalAfter = closest(0, false);
	const Record nrtBefore = closest(1, true), nrtAfter = closest(1, false);
	qCInfo(sunspots).noquote() << "Available" << (phi ? "PHI" : "HMI") << "observation endpoints for"
	                           << wantedText << (phi ? ": L2" : ": final")
	                           << observationOrNone(finalBefore) << "<|>" << observationOrNone(finalAfter)
	                           << (phi ? "; LL02" : "; NRT") << observationOrNone(nrtBefore)
	                           << "<|>" << observationOrNone(nrtAfter);
	const bool finalCoversTime = finalBefore.time.isValid() && finalAfter.time.isValid();
	const bool lowLatencyCoversTime = nrtBefore.time.isValid() && nrtAfter.time.isValid();
	if (phi && !requestedPhiJpeg2000_ && (latestForJob_ || (!finalCoversTime && !lowLatencyCoversTime)))
	{
		requestedPhiJpeg2000_ = true;
		pendingQueries_ = latestForJob_ ? 2 : 1;
		totalQueries_ += pendingQueries_;
		showProgress(ProgressKind::Fetching, 100 * (totalQueries_ - pendingQueries_) / totalQueries_);
		requestPhiJpeg2000(wantedText, latestForJob_);
		return;
	}
	const int series = finalCoversTime ? 0 : (phi && !lowLatencyCoversTime ? 2 : 1);
	const Record earlier = closest(series, true);
	const Record later = closest(series, false);
	if (phi && requestedPhiJpeg2000_)
		qCInfo(sunspots).noquote() << "Available PHI JPEG2000 observation endpoints for" << wantedText
		                           << ":" << observationOrNone(closest(2, true)) << "<|>"
		                           << observationOrNone(closest(2, false));
	Record selected;
	if (latestForJob_)
	{
		for (const auto& seriesResults : queryResults_)
			for (const auto& side : seriesResults)
				for (const Record& record : side)
				{
					const QDateTime time = phi ? record.startTime : record.time;
					const QDateTime selectedTime = phi ? selected.startTime : selected.time;
					if (!selected.time.isValid() || time > selectedTime ||
					    (time == selectedTime && seriesPriority(record.series) < seriesPriority(selected.series)))
						selected = record;
				}
	}
	else if (earlier.time.isValid() && later.time.isValid())
		selected = wanted - earlier.julianDayTT() <= later.julianDayTT() - wanted ? earlier : later;
	else selected = earlier.time.isValid() ? earlier : later;
	if (!selected.time.isValid())
	{
		emit statusChanged();
		qCWarning(sunspots).noquote() << "No usable" << (phi ? "PHI" : "HMI") << "observation selected for"
		                              << wantedText;
		const QString noObservation = phi
			? q_("No usable PHI continuum observation is available for this time.")
			: q_("No usable HMI continuum observation is available for this time.");
		showError(queryError_.isEmpty() ? noObservation : queryError_);
		finishJob();
		return;
	}
	qCInfo(sunspots).noquote() << "Selected" << selected.series << selected.observation;
	emit statusChanged();
	if (!tooFastWarning_) showError(queryError_);
	if (!forceImageUpdateForJob_ && autoDownloadForJob_ && prepareTextureForJob_ && isFast())
	{
		showTooFastWarning();
		finishJob();
	}
	else if (prepareTextureForJob_ &&
	         (forceImageUpdateForJob_ || displayedImageDeleted_ || imageNeedsUpdate(displayedRecord_)))
		useRecord(selected);
	else finishJob();
}

void Sunspots::useRecord(const Record& record)
{
	if (record.cacheBaseName() == canceledCacheBaseName_ ||
	    (!displayedImageDeleted_ && record.cacheBaseName() == displayedRecord_.cacheBaseName()))
	{
		finishJob();
		return;
	}
	canceledCacheBaseName_.clear();
	jobRecord_ = record;
	protectJobFiles(record);
	if (QFileInfo::exists(dataPath(record, ".png")))
	{
		startImageLoad(record);
		return;
	}
	if (QFileInfo::exists(dataPath(record, record.sourceSuffix())))
	{
		startTexturePreparation(record);
		return;
	}
	if (cacheOnlyForJob_)
	{
		refreshCache();
		cacheEmptyError_ = !hasCachedImagesForInstrument();
		showError(q_("Cached image is no longer available: %1").arg(record.cacheBaseName()));
		finishJob();
	}
	else acquireSource(record);
}

void Sunspots::acquireSource(const Record& record)
{
	Q_ASSERT(isBusy() && !workerRunning_ && !downloadReply_ && queryReplies_.isEmpty());
	if (!record.segment.isValid())
	{
		showError(q_("Observation has no downloadable continuum image."));
		finishJob();
		return;
	}
	const QString source = intermediatePath(record);
	qCInfo(sunspots).noquote() << "Source image download URL:" << record.segment.toString(QUrl::FullyEncoded);
	qCInfo(sunspots).noquote() << "Source image intermediate file:" << QFileInfo(source).absoluteFilePath();
	// Only accept the supported plain filenames, including when the archive supplied the PHI identifier.
	if (!recordFromIntermediateFileName(intermediateBaseName(record) + record.sourceSuffix()).time.isValid())
	{
		showError(q_("Invalid intermediate source image filename: %1").arg(source));
		finishJob();
		return;
	}
	// Keep the geometry even if acquisition fails, so a manually supplied FITS needs no sidecar from the user.
	if (record.instrument == Instrument::HMI && !saveMetadata(record, record.metadata))
	{
		finishJob();
		return;
	}
	if (QFileInfo::exists(source))
	{
		qCInfo(sunspots).noquote() << "Using intermediate source image file:" << source;
		finalizeSource(record);
	}
	else startDownload(record);
}

void Sunspots::startDownload(const Record& record)
{
	Q_ASSERT(isBusy() && !workerRunning_ && !downloadReply_ && queryReplies_.isEmpty());
	downloadFile_.reset(new QSaveFile(intermediatePath(record)));
	if (!downloadFile_->open(QIODevice::WriteOnly))
	{
		showError(q_("Cannot create source image cache file: %1").arg(downloadFile_->errorString()));
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
			showError(q_("Cannot write source image cache file: %1").arg(downloadFile_->errorString()));
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
			const QByteArray rest = reply->isOpen() ? reply->readAll() : QByteArray();
			if (!rest.isEmpty())
			{
				if (downloadFile_->write(rest) != rest.size())
				{
					downloadWriteFailed_ = true;
					showError(q_("Cannot write source image cache file: %1").arg(downloadFile_->errorString()));
				}
				downloadedBytes_ += rest.size();
			}
		}
		const bool received = reply->error() == QNetworkReply::NoError && !downloadWriteFailed_ &&
		                      downloadedBytes_ > 0 && downloadFile_;
		const bool good = received && downloadFile_->commit();
		if (received && !good)
			showError(q_("Cannot save source image cache file: %1").arg(downloadFile_->errorString()));
		if (!good && downloadFile_) downloadFile_->cancelWriting();
		if (!good && !downloadWriteFailed_ && !received)
			showError(downloadTimedOut_ ? q_("Source image download timed out.") :
			          reply->error() != QNetworkReply::NoError
			          ? q_("Source image download failed: %1").arg(reply->errorString())
			          : q_("Source image download is empty."));
		reply->deleteLater();
		downloadReply_.clear();
		downloadFile_.reset();
		if (good) finalizeSource(record);
		else finishJob();
	});
	idleTimer->start();
}

void Sunspots::finalizeSource(const Record& record)
{
	Q_ASSERT(!workerRunning_ && !downloadReply_);
	showProgress(ProgressKind::PreparingTexture, 0);
	const QString source = intermediatePath(record);
	workerRunning_ = true;
	workers_.start(QRunnable::create([this, record, source]
	{
		Record finished = record;
		QString error;
		try
		{
			const auto metadata = SunspotsImage::readMetadata(source);
			for (auto it = metadata.cbegin(); it != metadata.cend(); ++it)
				finished.metadata.insert(it.key(), it.value());
			if (record.instrument == Instrument::PHI)
			{
				if (finished.metadata.value("TIMESYS") != "UTC")
					throw std::runtime_error(
						q_("PHI source image must use UTC observation times.").toStdString());
				finished.startTime = parseUtc(finished.metadata.value("DATE-BEG"));
				finished.time = parseUtc(finished.metadata.value("DATE-AVG"));
				if (!finished.startTime.isValid() || !finished.time.isValid() ||
				    (record.series == phiJpeg2000Series
				     ? finished.startTime.toSecsSinceEpoch() != record.startTime.toSecsSinceEpoch()
				     : finished.startTime != record.startTime) ||
				    finished.metadata.value("INSTRUME") != "PHI")
					throw std::runtime_error(
						q_("PHI source image has missing or inconsistent observation times.").toStdString());
				finished.exactTime = true;
				finished.observation = finished.time.toString(Qt::ISODateWithMs);
			}
		}
		catch (const std::runtime_error& e) { error = QString::fromUtf8(e.what()); }
		QMetaObject::invokeMethod(this, [this, source, finished, error]
		{
			workerRunning_ = false;
			if (!error.isEmpty()) { showError(error); finishJob(); return; }
			const QString destination = dataPath(finished, finished.sourceSuffix());
			QFile file(source);
			if (!file.rename(destination))
			{
				showError(q_("Cannot finalize source image cache file %1: %2")
				          .arg(destination, file.errorString()));
				finishJob();
				return;
			}
			jobRecord_ = finished;
			protectJobFiles(finished);
			refreshCache();
			// Preserve the completed source even when cancellation arrives during header reading.
			if (cancellationToken_->load()) finishJob();
			else startTexturePreparation(finished);
		}, Qt::QueuedConnection);
	}));
}

void Sunspots::startTexturePreparation(const Record& record)
{
	Q_ASSERT(isBusy() && !workerRunning_ && !downloadReply_ && queryReplies_.isEmpty());
	showProgress(ProgressKind::PreparingTexture, 0);
	const QString sourcePath = dataPath(record, record.sourceSuffix());
	const QString pngPath = dataPath(record, ".png");
	QMap<QString, QString> metadata;
	if (record.instrument == Instrument::HMI) metadata = readMetadata(record);
	for (auto it = record.metadata.cbegin(); it != record.metadata.cend(); ++it)
		if (!metadata.contains(it.key())) metadata.insert(it.key(), it.value());
	const auto token = cancellationToken_;
	workerRunning_ = true;
	workers_.start(QRunnable::create([this, record, sourcePath, pngPath, metadata, token]
	{
		const auto result = SunspotsImage::prepareTexture(sourcePath, pngPath, metadata,
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
			else if (record.instrument == Instrument::PHI || saveMetadata(record, result.metadata))
			{
				removeImageSource(record);
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
			if (textureEnabled)
			{
				GETSTELMODULE(SolarSystem)->setTextureForPlanet("Sun", image);
				displayedRecord_ = record;
				displayedImageDeleted_ = false;
			}
			QStringList errors;
			if (!queryError_.isEmpty()) errors.append(queryError_);
			if (!textureCleanupError_.isEmpty()) errors.append(textureCleanupError_);
			showError(errors.join('\n'));
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
