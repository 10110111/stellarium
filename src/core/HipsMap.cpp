#include "HipsMap.hpp"
#include <QPainter>
#include "StelHealpix.hpp"

HipsMap::HipsMap()
{
	resize(500, 500);
	setWindowFlags(Qt::WindowStaysOnTopHint);
	connect(&timer, &QTimer::timeout, this, qOverload<>(&HipsMap::update));
	timer.start(100);
}

void HipsMap::tileCreated(unsigned order, unsigned npix, const StelTextureSP& texture)
{
	if(tiles.size() < order + 1)
		tiles.resize(order + 1);
	if(tiles[order].size() < npix + 1)
		tiles[order].resize(npix + 1);
	tiles[order][npix] = {Status::Loading, texture.toWeakRef()};
	if(!everShown)
	{
		show();
		everShown = true;
	}
	update();
}

void HipsMap::paintEvent(QPaintEvent* event)
{
	/*
	 * Order-0 layout:
	 *                         Npix8  Npix4
	 *                  Npix9  Npix5  Npix0
	 *          Npix10  Npix6  Npix1
	 *  Npix11  Npix7   Npix2
	 *          Npix3
	 */
	static constexpr std::pair<double, double> order0layout[12] = {
		{4, 1},
		{3, 2},
		{2, 3},
		{1, 4},

		{4, 0},
		{3, 1},
		{2, 2},
		{1, 3},

		{3, 0},
		{2, 1},
		{1, 2},
		{0, 3},
	};

	QPainter p(this);
	auto f = font();
	f.setPixelSize(7);
	p.setFont(f);
	p.setBrush(Qt::black);
	p.drawRect(rect());
	const auto scaleX = width() / 5.;
	const auto scaleY = height() / 5.;
	for(unsigned order = 0; order < tiles.size(); ++order)
	{
		const unsigned npixMax = (1 << (2*order)) * 12;
		for(unsigned npix = 0; npix < npixMax; ++npix)
		{
			if(npix >= tiles[order].size()) break;
			const int nside = 1 << order;
			int ix, iy, face;
			healpix_nest2xyf(nside, npix, &ix, &iy, &face);
			const double x = order0layout[face].first  + iy / double(nside);
			const double y = order0layout[face].second + ix / double(nside);
			QColor brush;
			auto& tile = tiles[order][npix];
			StelTextureSP tex = tile.texture.toStrongRef();
			if(!tex)
				tile.status = Status::Unrequested;
			else if(tex->canBind() || tex->dataLoaded())
				tile.status = Status::Ready;
			else
				tile.status = Status::Loading;
			const QColor readyColor(0, std::min(255u, 32 * (order+1)), 0);
			const QColor loadingColor(std::min(255u, 32 * (order+1)),
			                          std::min(128u, 16 * (order+1)),
			                          0);
			switch(tile.status)
			{
			case Status::Unrequested:
				if(order > 0) continue;
				brush = Qt::white;
				break;
			case Status::Loading:
				brush = loadingColor;
				break;
			case Status::Ready:
				brush = readyColor;
				break;
			default:
				brush = Qt::magenta;
			}
			p.setBrush(brush);
			p.setPen(Qt::transparent);
			const QRectF rect(QPointF(x * scaleX, y * scaleY), QSizeF(scaleX / nside, scaleY / nside));
			p.drawRect(rect);
			if(tile.status == Status::Loading && tex)
			{
				p.setBrush(readyColor);
				const auto [bytesDone, bytesToDo] = tex->getDownloadProgress();
				const auto progress = double(bytesDone) / bytesToDo;
				const QRectF progressBar(rect.topLeft(), QSizeF(rect.width() * progress, rect.height()));
				p.drawRect(progressBar);
			}
			p.setPen(Qt::white);
			p.drawText(rect, QString::number(npix), Qt::AlignCenter|Qt::AlignVCenter);
		}
	}
}
