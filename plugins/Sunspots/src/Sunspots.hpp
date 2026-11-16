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
	Q_PROPERTY(bool hmiPhotosphereEnabled READ isHmiPhotosphereEnabled WRITE setHmiPhotosphereEnabled NOTIFY hmiPhotosphereToggled)
	Q_PROPERTY(bool autoDownload READ isAutoDownload WRITE setAutoDownload NOTIFY autoDownloadChanged)
	Q_PROPERTY(int maxImageAgeMinutes READ getMaxImageAgeMinutes WRITE setMaxImageAgeMinutes
	           NOTIFY maxImageAgeMinutesChanged)
public:
	struct Record
	{
		QString series;
		QString observation;
		QDateTime time;
		QUrl segment;
		QMap<QString, QString> metadata;
		//! Whether a completed PNG was found when enumerating cached observations.
		bool hasCachedPng = false;
		//! Filename shared by the FITS, CSV, and PNG cache entries, without an extension.
		QString cacheBaseName() const;
	};
	//! Stages of the single active availability/image job.
	enum class ProgressKind { None, Fetching, Downloading, PreparingTexture, Saving, Loading, Canceling };

	Sunspots();
	~Sunspots() override;

	void init() override;
	bool configureGui(bool show = true) override;
	bool isDialogVisible() const { return dialogVisible; }
	bool isHmiPhotosphereEnabled() const { return hmiPhotosphereEnabled; }
	bool isAutoDownload() const { return autoDownload_; }
	int getMaxImageAgeMinutes() const { return maxImageAgeMinutes_; }
	//! Date and time of the image currently applied to the Sun, or N/A when none is shown.
	QString displayedTimeText() const;
	//! Signed age in real time, or distance before/after simulation time otherwise.
	QString displayedImageAgeText() const;
	//! Format an observation in the selected simulation time zone, using its date's offset.
	QString observationTimeText(const Record& record) const;
	//! Completed cache entries, with FITS and PNG files for an observation grouped together.
	const QList<Record>& cachedRecords() const { return cachedRecords_; }
	//! Total size in bytes of all FITS, PNG, and CSV files at the last cache scan.
	qint64 cacheSizeBytes() const { return cacheSizeBytes_; }
	//! Observation being prepared, or an invalid record while idle or fetching availability.
	Record activeRecord() const { return isBusy() ? jobRecord_ : Record{}; }
	//! Identity of the observation currently displayed, or an empty string.
	QString displayedCacheBaseName() const;
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
	void setHmiPhotosphereEnabled(bool enable);
	void setAutoDownload(bool enable);
	//! Set the maximum absolute distance between simulation time and the displayed image's T_OBS, in minutes.
	void setMaxImageAgeMinutes(int minutes);
	//! Fetch availability and apply the nearest observation regardless of the current image's age.
	void downloadForSimTime();
	//! Disable automatic updates and apply the latest observation across both HMI series.
	void downloadLatest();
	//! Disable automatic updates and display this cached observation without an availability query.
	void displayCachedImage(const QString& cacheBaseName);
	//! Rescan completed FITS/PNG cache entries, for example when opening the dialog.
	void refreshCache();
	//! Cancel the current job and automatic updates without changing the displayed texture.
	void cancelCurrentOperation();

signals:
	void dialogToggled(bool enabled);
	void hmiPhotosphereToggled(bool enabled);
	void autoDownloadChanged(bool enabled);
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
	void acceptAvailability();
	void useRecord(const Record& record);
	void startDownload(const Record& record);
	void startTexturePreparation(const Record& record);
	void startImageLoad(const Record& record);
	void showError(const QString& error);
	void showTooFastWarning();
	void beginJob(ProgressKind kind);
	void showProgress(ProgressKind kind, int value);
	void removeDownloadProgressBar();
	void finishJob();
	QString dataPath(const Record& record, const QString& suffix) const;
	bool saveMetadata(const Record& record, const QMap<QString, QString>& metadata);
	QMap<QString, QString> readMetadata(const Record& record) const;
	bool isFast() const;
	bool imageNeedsUpdate(const Record& record) const;
	Record usableCachedRecord() const;

	QSettings* conf_;
	bool dialogVisible = false;
	bool hmiPhotosphereEnabled = false;
	bool autoDownload_ = false;
	int maxImageAgeMinutes_ = 180;
	//! Keep preparing the texture if Show is switched off after the job has begun.
	bool prepareTextureForJob_ = false;
	bool forceImageUpdateForJob_ = false;
	bool latestForJob_ = false;
	QString cacheDir_;
	QString error_;
	QString queryError_;
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
	QList<Record> cachedRecords_;
	qint64 cacheSizeBytes_ = 0;
	Record jobRecord_;
	QDateTime requestedTime_;
	Record displayedRecord_;
	Record lastReadyRecord_;
	QList<Record> queryResults_[2][2];
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
