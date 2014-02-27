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
#include <cs/Spin_3.h>

namespace CS
{
template<>
void check_spin_3_norm<float>(const float &s12, const float &s23, const float &s31, const float &s0)
{
    // double is an inexact type - check norm with an error bound
    const float MAX_ERROR = 10e-6;
    (void)MAX_ERROR;

    const float length = std::sqrt(s12 * s12 + s23 * s23 + s31 * s31 + s0 * s0);
    (void)length;

    assert(std::fabs(1.0 - length) <= MAX_ERROR);
}

template<>
void check_spin_3_norm<double>(const double &s12, const double &s23, const double &s31, const double &s0)
{
    // double is an inexact type - check norm with an error bound
    const double MAX_ERROR = 10e-10;
    (void)MAX_ERROR;

    const double length = std::sqrt(s12 * s12 + s23 * s23 + s31 * s31 + s0 * s0);
    (void)length;

    assert(std::fabs(1.0 - length) <= MAX_ERROR);
}

template<>
void check_spin_3_norm<CGAL::Gmpfr>(const CGAL::Gmpfr &s12, const CGAL::Gmpfr &s23, const CGAL::Gmpfr &s31, const CGAL::Gmpfr &s0)
{
    // CGAL::Gmpfr is an inexact type - check norm with an error bound
    const CGAL::Gmpfr MAX_ERROR = 10e-64;

    const CGAL::Gmpfr length = Math::sqrt(s12 * s12 + s23 * s23 + s31 * s31 + s0 * s0);
    //Math::print(length);

    assert(Math::fabs(CGAL::Gmpfr(1.0) - length) <= MAX_ERROR);
}
} // namespace CS
