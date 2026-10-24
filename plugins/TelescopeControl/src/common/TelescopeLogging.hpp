/*
Stellarium Telescope Control Plug-in

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this library; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Suite 500, Boston, MA  02110-1335, USA.
*/

#ifndef TELESCOPELOGGING_HPP
#define TELESCOPELOGGING_HPP

#include <QDebug>
#include <QLoggingCategory>

Q_DECLARE_LOGGING_CATEGORY(Telescopes)

//! Slot identification for telescope diagnostics. A nonpositive slot omits the prefix.
struct TelescopeLogPrefix
{
	explicit TelescopeLogPrefix(int slot) : slot(slot) {}
	int slot;
};

//! Writes the slot prefix, with its format defined in one place, preserving stream spacing.
QDebug operator<<(QDebug debug, TelescopeLogPrefix prefix);

#endif // TELESCOPELOGGING_HPP
