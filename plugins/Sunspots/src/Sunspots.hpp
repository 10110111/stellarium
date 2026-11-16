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

#ifndef SUNSPOTS_HPP
#define SUNSPOTS_HPP

#include <atomic>
#include <memory>
#include "StelModule.hpp"
#include <QDateTime>
#include <QLoggingCategory>
#include <QMap>
#include <QPointer>
#include <QSet>
#include <QStringList>
#include <QThreadPool>
#include <QTimer>
#include <QUrl>

Q_DECLARE_LOGGING_CATEGORY(sunspots)

class QSettings;
class StelButton;
class StelProgressController;
class SunspotsDialog;
class QNetworkReply;
class QSaveFile;
class QImage;
class Sunspots : public StelModule
{
	Q_OBJECT
	Q_PROPERTY(bool dialogVisible READ isDialogVisible WRITE setDialogVisible NOTIFY dialogToggled)
	Q_PROPERTY(bool textureEnabled READ isTextureEnabled WRITE setTextureEnabled NOTIFY textureToggled)
	Q_PROPERTY(Instrument instrument READ getInstrument WRITE setInstrument NOTIFY instrumentChanged)
	Q_PROPERTY(bool autoDownload READ isAutoDownload WRITE setAutoDownload NOTIFY autoDownloadChanged)
	Q_PROPERTY(bool continuousUseCacheOnly READ isContinuousUseCacheOnly WRITE setContinuousUseCacheOnly
	           NOTIFY continuousUseCacheOnlyChanged)
	Q_PROPERTY(bool deleteSourceAfterMakingTexture READ isDeleteSourceAfterMakingTexture
	           WRITE setDeleteSourceAfterMakingTexture
	           NOTIFY deleteSourceAfterMakingTextureChanged)
	Q_PROPERTY(int maxImageAgeMinutes READ getMaxImageAgeMinutes WRITE setMaxImageAgeMinutes
	           NOTIFY maxImageAgeMinutesChanged)
public:
	//! Source of photospheric continuum observations.
	enum class Instrument { HMI, PHI };
	Q_ENUM(Instrument)

	struct Record
	{
		Instrument instrument = Instrument::HMI;
		QString series;
		QString observation;
		//! HMI stores TAI calendar fields; PHI stores genuine UTC. Use julianDayTT() for comparisons.
		QDateTime time;
		//! Archive identity for PHI, known before downloading its exact DATE-AVG.
		QDateTime startTime;
		//! Whether PHI's DATE-AVG has been read from the source or recovered from the cache filename.
		bool exactTime = false;
		QUrl segment;
		QMap<QString, QString> metadata;
		//! Whether a completed PNG was found when enumerating cached observations.
		bool hasCachedPng = false;
		//! Shared cache filename without an extension; only HMI requires a CSV sidecar.
		QString cacheBaseName() const;
		//! Observation instant in TT, following Stellarium's configured UTC/TT conversion.
		double julianDayTT() const;
		//! Instrument identifier used for settings and logging.
		QString instrumentName() const;
		//! Source column label, including the NRT, LL, or JPEG2000 qualifier for preliminary products.
		QString sourceName() const;
		//! Extension of the downloaded image, including the leading dot.
		QString sourceSuffix() const;
	};
	//! Stages of the single active availability/image job.
	enum class ProgressKind { None, Fetching, Downloading, PreparingTexture, Saving, Loading, Canceling };

	Sunspots();
	~Sunspots() override;

	void init() override;
	bool configureGui(bool show = true) override;
	bool isDialogVisible() const { return dialogVisible; }
	bool isTextureEnabled() const { return textureEnabled; }
	//! Instrument used for availability lookups and automatic cache selection.
	Instrument getInstrument() const { return instrument_; }
	bool isAutoDownload() const { return autoDownload_; }
	//! Whether automatic selection is restricted to cached images.
	bool isContinuousUseCacheOnly() const { return continuousUseCacheOnly_; }
	//! Whether to remove the downloaded source after successfully caching its texture and metadata.
	bool isDeleteSourceAfterMakingTexture() const { return deleteSourceAfterMakingTexture_; }
	//! Whether a cache-only selection failed because the cache was empty.
	bool hasEmptyCacheError() const
	{
		return cacheEmptyError_ && textureEnabled && continuousUseCacheOnly_;
	}
	int getMaxImageAgeMinutes() const { return maxImageAgeMinutes_; }
	//! Translated title identifying the displayed image's instrument, or the selected instrument when none is shown.
	QString displayedImageTitle() const;
	//! Date and time of the image currently applied to the Sun, or N/A when none is shown.
	QString displayedTimeText() const;
	//! Signed age in real time, or distance before/after simulation time otherwise.
	QString displayedImageAgeText() const;
	//! Format an observation in the selected simulation time zone, using its date's offset.
	QString observationTimeText(const Record& record) const;
	//! Completed cache entries, with source and PNG files for an observation grouped together.
	const QList<Record>& cachedRecords() const { return cachedRecords_; }
	//! Total size in bytes of all FITS, JP2, PNG, and CSV files at the last cache scan.
	qint64 cacheSizeBytes() const { return cacheSizeBytes_; }
	//! Observation being prepared, or an invalid record while idle or fetching availability.
	Record activeRecord() const { return isBusy() ? jobRecord_ : Record{}; }
	//! Identity of the observation currently displayed, or an empty string.
	QString displayedCacheBaseName() const;
	//! Whether this completed cache entry can replace the displayed image without interrupting a job.
	//! If supplied, disabledReason receives a translated explanation, or an empty string when usable.
	bool canDisplayCachedImage(const QString& cacheBaseName, QString* disabledReason = nullptr) const;
	//! Whether a cached or active entry can be deleted; cancellation already in progress cannot be repeated.
	//! If supplied, disabledReason receives a translated explanation, or an empty string when deletable.
	bool canDeleteCachedImage(const QString& cacheBaseName, QString* disabledReason = nullptr) const;
	QString errorText() const;
	QString progressText() const;
	int progressValue() const { return progressValue_; }
	//! Current stage, also used by the cache row's progress controls.
	ProgressKind progressKind() const { return progressKind_; }
	//! Whether an update or texture preparation job is running.
	bool isBusy() const { return progressKind_ != ProgressKind::None; }
	//! Whether cancellation is waiting for the image worker to stop.
	bool isCanceling() const { return progressKind_ == ProgressKind::Canceling; }

public slots:
	void setDialogVisible(bool enable);
	void setTextureEnabled(bool enable);
	//! Select the source for subsequent jobs. Ignored while a job or cache clearing is active.
	void setInstrument(Instrument instrument);
	void setAutoDownload(bool enable);
	//! Change the policy for subsequent automatic selections, including while a job is running.
	void setContinuousUseCacheOnly(bool enable);
	//! Change the cleanup policy, which is read when a running job reaches source cleanup.
	void setDeleteSourceAfterMakingTexture(bool enable);
	//! Confirm cache clearing, preserving the displayed image and files needed by ongoing work.
	void clearCache();
	//! Confirm removal of an observation's files, canceling and waiting for its active job first.
	//! If it is displayed and Show is enabled, select a replacement using the current selection policy.
	void deleteCachedImage(QString cacheBaseName);
	//! Set the maximum absolute distance between simulation time and the displayed observation, in minutes.
	void setMaxImageAgeMinutes(int minutes);
	//! Fetch availability and apply the nearest observation regardless of the current image's age.
	void downloadForSimTime();
	//! Disable automatic updates and apply the latest observation from the selected instrument.
	void downloadLatest();
	//! Disable automatic updates and display this cached observation without an availability query.
	void displayCachedImage(const QString& cacheBaseName);
	//! Rescan completed source/PNG cache entries, for example when opening the dialog.
	void refreshCache();
	//! Cancel the current job and automatic updates without changing the displayed texture.
	void cancelCurrentOperation();

signals:
	void dialogToggled(bool enabled);
	void textureToggled(bool enabled);
	void instrumentChanged(Instrument instrument);
	void autoDownloadChanged(bool enabled);
	void continuousUseCacheOnlyChanged(bool enabled);
	void deleteSourceAfterMakingTextureChanged(bool enabled);
	void maxImageAgeMinutesChanged(int minutes);
	void statusChanged();
	//! Completed cache entries may have been added or removed.
	void cacheChanged();

private:
	void loadSettings();
	void pollTime();
	void ensureUsableImage();
	void queryAvailability(bool userRequested = false, bool forceImageUpdate = false, bool latest = false);
	void requestRecord(const QString& series, bool before, const QString& time);
	void requestPhiRecords(bool lowLatency);
	void requestPhiJpeg2000(const QString& date, bool findLatest);
	void finalizeSource(const Record& record);
	void acceptAvailability();
	void useRecord(const Record& record);
	void acquireSource(const Record& record);
	void startDownload(const Record& record);
	void startTexturePreparation(const Record& record);
	void startImageLoad(const Record& record);
	void showError(const QString& error);
	void showTooFastWarning();
	void beginJob(ProgressKind kind);
	void showProgress(ProgressKind kind, int value);
	void removeDownloadProgressBar();
	void finishJob();
	void removeCachedFiles(const QSet<QString>& protectedBaseNames);
	QStringList removeCachedImageFiles(const QString& cacheBaseName);
	void resumeSelectionAfterDeletion();
	void removeImageSource(const Record& record);
	QString dataPath(const Record& record, const QString& suffix) const;
	QString intermediatePath(const Record& record) const;
	void protectJobFiles(const Record& record);
	bool saveMetadata(const Record& record, const QMap<QString, QString>& metadata);
	QMap<QString, QString> readMetadata(const Record& record) const;
	bool isFast() const;
	bool imageNeedsUpdate(const Record& record) const;
	bool hasCachedImagesForInstrument() const;
	Record usableCachedRecord(bool allowStale = false) const;

	QSettings* conf_;
	bool dialogVisible = false;
	bool textureEnabled = false;
	Instrument instrument_ = Instrument::HMI;
	bool autoDownload_ = false;
	bool continuousUseCacheOnly_ = false;
	bool deleteSourceAfterMakingTexture_ = false;
	bool cacheEmptyError_ = false;
	bool cacheClearPending_ = false;
	//! Jobs that finish or select their observation during the confirmation must remain protected.
	QSet<QString> cacheClearProtectedBaseNames_;
	int maxImageAgeMinutes_ = 180;
	//! Keep preparing the texture if Show is switched off after the job has begun.
	bool prepareTextureForJob_ = false;
	bool forceImageUpdateForJob_ = false;
	bool latestForJob_ = false;
	bool autoDownloadForJob_ = false;
	bool cacheOnlyForJob_ = false;
	QString cacheDir_;
	QString error_;
	QString queryError_;
	QString textureCleanupError_;
	bool tooFastWarning_ = false;
	ProgressKind progressKind_ = ProgressKind::None;
	int progressValue_ = 0;
	int totalQueries_ = 0;
	int pendingQueries_ = 0;
	bool workerRunning_ = false;
	qint64 downloadedBytes_ = 0;
	qint64 downloadReceivedBytes_ = 0;
	qint64 downloadTotalBytes_ = 0;
	bool downloadWriteFailed_ = false;
	bool downloadTimedOut_ = false;
	QString canceledCacheBaseName_;
	//! Original identity of an active entry to remove once its worker has stopped (PHI may rename it).
	QString jobCacheBaseNameToDelete_;
	bool selectionAfterDeletionPending_ = false;
	//! Keep the old texture/date visible while forcing selection to replace a deleted displayed entry.
	bool displayedImageDeleted_ = false;
	QList<Record> cachedRecords_;
	qint64 cacheSizeBytes_ = 0;
	Record jobRecord_;
	//! Snapshot of the simulation time for all queries belonging to one availability job.
	double requestedTimeTT_ = 0.;
	Record displayedRecord_;
	Record lastReadyRecord_;
	//! Final/preliminary/quicklook products; PHI catalogs use only the first side of each product.
	QList<Record> queryResults_[3][2];
	bool requestedPhiJpeg2000_ = false;
	QPointer<QNetworkReply> downloadReply_;
	QSet<QNetworkReply*> queryReplies_;
	std::unique_ptr<QSaveFile> downloadFile_;
	QThreadPool workers_;
	QTimer autoUpdateTimer_;
	QTimer timeChangeTimer_;
	std::shared_ptr<std::atomic_bool> cancellationToken_ = std::make_shared<std::atomic_bool>(false);
	StelButton* toolbarButton_ = nullptr;
	StelProgressController* downloadProgressBar_ = nullptr;
	std::unique_ptr<SunspotsDialog> dialog_;
};


#include <QObject>
#include "StelPluginInterface.hpp"

class SunspotsStelPluginInterface : public QObject, public StelPluginInterface
{
	Q_OBJECT
	Q_PLUGIN_METADATA(IID StelPluginInterface_iid)
	Q_INTERFACES(StelPluginInterface)
public:
	StelModule* getStelModule() const override;
	StelPluginInfo getPluginInfo() const override;
};

#endif
