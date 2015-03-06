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
#ifndef LIBCS_LINEAR_SYSTEM_H
#define LIBCS_LINEAR_SYSTEM_H

namespace CS
{
// A x = b
template<class FT_>
bool linear_solve(FT_ a00, FT_ a01,
                  FT_ a10, FT_ a11,
                  FT_ b0, FT_ b1,
                  FT_ &x0, FT_ &x1);

// A x = b
template<class FT_>
bool linear_solve(FT_ a00, FT_ a01, FT_ a02,
                  FT_ a10, FT_ a11, FT_ a12,
                  FT_ a20, FT_ a21, FT_ a22,
                  FT_ b0, FT_ b1, FT_ b2,
                  FT_ &x0, FT_ &x1, FT_ &x2);

} // namespace CS

#include "Linear_system.ipp"

#endif // LIBCS_LINEAR_SYSTEM_H
