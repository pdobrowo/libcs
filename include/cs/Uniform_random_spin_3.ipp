/**
 * Copyright (C) 2009-2013  Przemysław Dobrowolski
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
#include "Uniform_random_spin_3.h"

namespace CS
{
template<typename FT_>
FT_ uniform_rand()
{
    // generate uniform random on [0; 1]
    return FT_(rand()) / FT_(RAND_MAX);
}

template<class FT_>
void uniform_random_spin_3(Spin_3<FT_> &out)
{
    using Math::pi;

    // based on:
    //
    // K. Shoemake.
    // Uniform random rotations.
    // In D. Kirk, editor, Graphics Gems III, pages 124-132. Academic, New York, 1992.
    FT_ u1, u2, u3;

    u1 = uniform_rand<FT_>();
    u2 = uniform_rand<FT_>();
    u3 = uniform_rand<FT_>();

    const FT_ ONE    = 1.0;
    const FT_ TWO    = 2.0;
    const FT_ TWO_PI = TWO * pi<FT_>();

    // take math functions from std
    using std::sqrt;
    using std::sin;
    using std::cos;

    // also take math functions from private implementation for CGAL types
    using Math::sqrt;
    using Math::sin;
    using Math::cos;

    // calculate
    out = Spin_3<FT_>(sqrt(ONE - u1) * sin(TWO_PI * u2),
                     sqrt(ONE - u1) * cos(TWO_PI * u2),
                     sqrt(u1) * sin(TWO_PI * u3),
                     sqrt(u1) * cos(TWO_PI * u3));
}
} // namespace CS
