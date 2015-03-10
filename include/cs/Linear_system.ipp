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
#include "Linear_system.h"

namespace CS
{
template<class FT_>
bool linear_solve(FT_ a00, FT_ a01,
                  FT_ a10, FT_ a11,
                  FT_ b0, FT_ b1,
                  FT_ &x0, FT_ &x1)
{
    FT_ det = a00 * a11 - a01 * a10;

    if (det == FT_(0))
        return false;

    x0 = ( a11 * b0 - a01 * b1 ) / det;
    x1 = ( a00 * b1 - a10 * b0 ) / det;

    return true;
}

template<class FT_>
bool linear_solve(FT_ a00, FT_ a01, FT_ a02,
                  FT_ a10, FT_ a11, FT_ a12,
                  FT_ a20, FT_ a21, FT_ a22,
                  FT_ b0, FT_ b1, FT_ b2,
                  FT_ &x0, FT_ &x1, FT_ &x2)
{
    FT_ det = a02 * a11 * a20 - a01 * a12 * a20 - a02 * a10 * a21 + a00 * a12 * a21 + a01 * a10 * a22 - a00 * a11 * a22;

    if (det == FT_(0))
        return false;

    x0 = ( (a12 * a21 - a11 * a22) * b0 + (a01 * a22 - a02 * a21 ) * b1 + (a02 * a11 - a01 * a12) * b2) / det;
    x1 = ( (a10 * a22 - a12 * a20) * b0 + (a02 * a20 - a00 * a22 ) * b1 + (a00 * a12 - a02 * a10) * b2) / det;
    x2 = ( (a11 * a20 - a10 * a21) * b0 + (a00 * a21 - a01 * a20 ) * b1 + (a01 * a10 - a00 * a11) * b2) / det;

    return true;
}
} // namespace CS
