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
#ifndef LIBCS_RANDOM_SPIN_CIRCLE_3_H
#define LIBCS_RANDOM_SPIN_CIRCLE_3_H

#include <CGAL/Kernel/global_functions.h>
#include <CGAL/Gmpfr.h>
#include "Spin_3.h"
#include "Quartic.h"
#include "Cubic.h"
#include "Quadratic.h"
#include <cassert>

namespace CS
{
// Random_spin_circle_3:
//
// Note: this class is designed to work with inexact types only
//
// g - random spin rotation
// g := g12 * e12 + g23 * e23 + g31 * e31 + g0
//
// Gamma_g(alpha) := (g12 * cos(alpha) + g0  * sin(alpha)) * e12 +
//                   (g23 * cos(alpha) - g31 * sin(alpha)) * e23 +
//                   (g31 * cos(alpha) + g23 * sin(alpha)) * e31 +
//                   (g0  * cos(alpha) - g12 * sin(alpha))
//
// intersection with quadric is:
//
// P * Cos[a]^2 + 2 * Q * Cos[a] * Sin[a] + R * Sin[a]^2 = 0
//
// where:
// P = a44 * g0^2 + 2 a14 * g0 * g12 + a11 * g12^2 + 2 a24 * g0 * g23 + 2 a12 * g12 * g23 +
//     a22 * g23^2 + 2 a34 * g0 * g31 + 2 a13 * g12 * g31 + 2 a23 * g23 * g31 + a33 * g31^2
//
// Q = a14 * g0^2 + a11 * g0 * g1- a44 * g0 * g1- a14 * g12^2 + a1g0 * g23 + a34 * g0 * g23 + a13 * g1g23 - a24 * g1g23 +
//     a23 * g23^2 + a13 * g0 * g31 - a24 * g0 * g31 - a1g1g31 - a34 * g1g31 - a2g23 * g31 + a33 * g23 * g31 - a23 * g31^2
//
// R = a11 * g0^2 - 2 a14 * g0 * g12 + a44 * g12^2 + 2 a13 * g0 * g23 - 2 a34 * g12 * g23 +
//     a33 * g23^2 - 2 a12 * g0 * g31 + 2 a24 * g12 * g31 - 2 a23 * g23 * g31 + a22 * g31^2
//
// quadric coefficients are defined as in Spin_quadric_3
//
// intersection equation is reduced to:
//
// P * t^4 - 4 * Q * t^3 + (4 * R - 2 * P) * t^2 + 4 * Q * t + P = 0
//
// alpha = 2 * atan(t)
//
template<class FT_>
int solve_general_quartic(const FT_ &a, const FT_ &b, const FT_ &c, const FT_ &d, const FT_ &e,
                          FT_ *x0, FT_ *x1, FT_ *x2, FT_ *x3);

template<class Kernel_>
class Random_spin_circle_3
{
    typedef typename Kernel_::FT                FT;

    typedef typename Kernel_::Spin_3            Spin_3;
    typedef typename Kernel_::Spin_quadric_3    Spin_quadric_3;

public:
    Random_spin_circle_3(const Spin_3 &rotation);

    // intersect with a spin quadric and return number of real intersections
    size_t intersect_quadric(const Spin_quadric_3 &quadric, FT out_solutions[4]) const;

    // circle evaluation at given value
    Spin_3 evaluate(const FT &x) const;

private:
    Spin_3      m_g;
};
} // namespace CS

#include "Random_spin_circle_3.ipp"

#endif // LIBCS_RANDOM_SPIN_CIRCLE_3_H
