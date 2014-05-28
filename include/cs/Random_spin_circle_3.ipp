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
#include "Random_spin_circle_3.h"

namespace CS
{
template<class FT_>
int solve_general_quartic(const FT_ &a, const FT_ &b, const FT_ &c, const FT_ &d, const FT_ &e,
                          FT_ *x0, FT_ *x1, FT_ *x2, FT_ *x3)
{
    if (a != 0)
        return Math::solve_quartic(b / a, c / a, d / a, e / a, x0, x1, x2, x3);
    else if (b != 0)
        return Math::solve_cubic(c / b, d / b, e / b, x0, x1, x2);
    else
        return Math::solve_quadratic(c, d, e, x0, x1);
}

template<class Kernel_>
Random_spin_circle_3<Kernel_>::Random_spin_circle_3(const Spin_3 &rotation)
    : m_g(rotation)
{
}

template<class Kernel_>
size_t Random_spin_circle_3<Kernel_>::intersect_quadric(const Spin_quadric_3 &quadric, FT out_solutions[4]) const
{
    // take quadric coefficients
    const FT &a11 = quadric.a11();
    const FT &a22 = quadric.a22();
    const FT &a33 = quadric.a33();
    const FT &a44 = quadric.a44();
    const FT &a12 = quadric.a12();
    const FT &a13 = quadric.a13();
    const FT &a14 = quadric.a14();
    const FT &a23 = quadric.a23();
    const FT &a24 = quadric.a24();
    const FT &a34 = quadric.a34();

    // take random spin coefficients
    const FT &g12 = m_g.s12();
    const FT &g23 = m_g.s23();
    const FT &g31 = m_g.s31();
    const FT &g0 = m_g.s0();

    // calculate P, Q, R coefficients
    const FT TWO = FT(2.0);
    const FT FOUR = FT(4.0);

    const FT P =       a44 * g0  * g0 +
                 TWO * a14 * g0  * g12 +
                       a11 * g12 * g12 +
                 TWO * a24 * g0  * g23 +
                 TWO * a12 * g12 * g23 +
                       a22 * g23 * g23 +
                 TWO * a34 * g0  * g31 +
                 TWO * a13 * g12 * g31 +
                 TWO * a23 * g23 * g31 +
                       a33 * g31 * g31;

    const FT Q = a14 * g0  * g0 +
                 a11 * g0  * g12 -
                 a44 * g0  * g12 -
                 a14 * g12 * g12 +
                 a12 * g0  * g23 +
                 a34 * g0  * g23 +
                 a13 * g12 * g23 -
                 a24 * g12 * g23 +
                 a23 * g23 * g23 +
                 a13 * g0  * g31 -
                 a24 * g0  * g31 -
                 a12 * g12 * g31 -
                 a34 * g12 * g31 -
                 a22 * g23 * g31 +
                 a33 * g23 * g31 -
                 a23 * g31 * g31;

    const FT R =       a11 * g0  * g0 -
                 TWO * a14 * g0  * g12 +
                       a44 * g12 * g12 +
                 TWO * a13 * g0  * g23 -
                 TWO * a34 * g12 * g23 +
                       a33 * g23 * g23 -
                 TWO * a12 * g0  * g31 +
                 TWO * a24 * g12 * g31 -
                 TWO * a23 * g23 * g31 +
                       a22 * g31 * g31;

    // solve base equation
    FT solutions[4];

    int number_of_solutions = solve_general_quartic(P, -FOUR * Q, FOUR * R - TWO * P, FOUR * Q, P,
                                                    solutions + 0, solutions + 1, solutions + 2, solutions + 3);

    // take inexact atan functions
    using std::atan;
    using Math::atan;

    // convert solutions to alpha-space
    for (int i = 0; i < number_of_solutions; ++i)
        out_solutions[i] = TWO * atan(solutions[i]);

    return static_cast<size_t>(number_of_solutions);
}

template<class Kernel_>
typename Random_spin_circle_3<Kernel_>::Spin_3 Random_spin_circle_3<Kernel_>::evaluate(const FT &x) const
{
    // take random spin coefficients
    const FT &g12 = m_g.s12();
    const FT &g23 = m_g.s23();
    const FT &g31 = m_g.s31();
    const FT &g0 = m_g.s0();

    // take inexact sin/cos functions
    using std::sin;
    using std::cos;
    using Math::sin;
    using Math::cos;

    // sin / cos argument
    const FT cx = cos(x);
    const FT sx = sin(x);

    return Spin_3(g12 * cx + g0  * sx,
                  g23 * cx - g31 * sx,
                  g31 * cx + g23 * sx,
                  g0  * cx - g12 * sx);
}
} // namespace CS
