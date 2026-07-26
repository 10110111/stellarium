#pragma once

#include <QTimer>
#include <QCache>
#include <QWidget>
#include <StelHips.hpp>

class HipsMap : public QWidget
{
	enum class Status
	{
		Unrequested,
		Loading,
		Ready,
	};
	struct Tile
	{
		Status status;
		QWeakPointer<StelTexture> texture;
	};
	QTimer timer;
	std::vector<std::vector<Tile>> tiles; // tiles[order][npix]
	bool everShown = false;
public:
	HipsMap();
	void tileCreated(unsigned order, unsigned npix, const StelTextureSP& texture);
protected:
	void paintEvent(QPaintEvent* event) override;
};
