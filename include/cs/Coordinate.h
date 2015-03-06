/**
 * Copyright (C) 2009-2015  Przemysław Dobrowolski
 *
 * This file is part of the Configuration Space Library (libcs), a library
 * for creating configuration spaces of various motion planning problems.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef LIBCS_COORDINATE_H
#define LIBCS_COORDINATE_H

#include <ostream>
#include <vector>

namespace CS
{
// Cell index is a bit vector
//
// Tested:
// vector<bool> - 64 sec
// vector<int>  - 38 sec
// vector<char> - 36 sec
typedef std::vector<char> Coordinate;

std::ostream &operator <<(std::ostream &stream, const Coordinate &coordinate);
} // namespace CS

#endif // LIBCS_COORDINATE_H
