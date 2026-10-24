/*
 * Stellarium Telescope Control Plug-in
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

#include "Connection.hpp"
#include "SerialPort.hpp"
#include "Server.hpp"
#include "TelescopeLogging.hpp"

#include <QDir>
#include <QEventLoop>
#include <QTemporaryDir>
#include <QTimer>
#include <QtTest>

class LoggingServer : public Server
{
public:
	explicit LoggingServer(int slot) : Server(0, slot) {}

private:
	void gotoReceived(unsigned int, int) override {}
	void syncReceived(unsigned int, int) override {}
};

class LoggingConnection : public Connection
{
public:
	explicit LoggingConnection(Server &server) : Connection(server, INVALID_SOCKET) {}

	void failReading()
	{
		performReading();
	}
};

class TestTelescopeLogging : public QObject
{
	Q_OBJECT

private slots:
	void init()
	{
		QLoggingCategory::setFilterRules("stel.Telescopes.debug=false\nstel.Telescopes.info=true\n"
			"stel.Telescopes.warning=true");
	}

	void cleanup()
	{
		QLoggingCategory::setFilterRules(QString());
	}

	void categoryFiltering()
	{
		int evaluations = 0;
		qCDebug(Telescopes) << TelescopeLogPrefix(23) << ++evaluations;
		QCOMPARE(evaluations, 0);

		QLoggingCategory::setFilterRules("stel.Telescopes.debug=true");
		QTest::ignoreMessage(QtDebugMsg, "[23] 1");
		qCDebug(Telescopes) << TelescopeLogPrefix(23) << ++evaluations;
		QCOMPARE(evaluations, 1);
	}

	void prefixPreservesFormatting()
	{
		QTest::ignoreMessage(QtInfoMsg, "[23] message 42");
		qCInfo(Telescopes) << TelescopeLogPrefix(23) << "message" << 42;
		QTest::ignoreMessage(QtInfoMsg, "[23] message42");
		qCInfo(Telescopes).nospace() << TelescopeLogPrefix(23) << "message" << 42;
		QTest::ignoreMessage(QtInfoMsg, "message 42");
		qCInfo(Telescopes) << TelescopeLogPrefix(0) << "message" << 42;
	}

	void interleavedTransportDiagnostics()
	{
		LoggingServer first(23);
		LoggingServer second(7);
		LoggingConnection firstConnection(first);
		LoggingConnection secondConnection(second);

		QTest::ignoreMessage(QtWarningMsg, QRegularExpression(
			"^\\[7\\] Connection::performReading: readNonblocking failed: "));
		QTest::ignoreMessage(QtWarningMsg, QRegularExpression(
			"^\\[23\\] Connection::performReading: readNonblocking failed: "));
		QEventLoop loop;
		QTimer::singleShot(0, &loop, [&]()
		{
			firstConnection.failReading();
			loop.quit();
		});
		secondConnection.failReading();
		loop.exec();
	}

	void constructorDiagnostics()
	{
		QTemporaryDir directory;
		QVERIFY(directory.isValid());
		LoggingServer server(23);
		QTest::ignoreMessage(QtWarningMsg, QRegularExpression(
			"^\\[23\\] SerialPort::SerialPort\\(.*\\): (open|CreateFile)\\(\\) failed: "));
		SerialPort port(server, qPrintable(directory.filePath("missing-device")));
		QVERIFY(port.isClosed());
		QVERIFY(QDir(directory.path()).entryList(QDir::Files).isEmpty());
	}
};

QTEST_GUILESS_MAIN(TestTelescopeLogging)
#include "testTelescopeLogging.moc"
