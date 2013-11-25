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
#include "Cubic.h"

namespace CS
{
namespace Math
{
template<class FT>
int solve_cubic(const FT &a, const FT &b, const FT &c,
                       FT *x0, FT *x1, FT *x2)
{
    const FT ZERO  = 0.0;
    const FT ONE   = 1.0;
    const FT TWO   = 2.0;
    const FT THREE = 3.0;
    const FT NINE  = 9.0;

    FT q = (a * a - THREE * b);
    FT r = (TWO * a * a * a - NINE * a * b + FT(27) * c);

    FT Q = q / NINE;
    FT R = r / FT(54);

    FT Q3 = Q * Q * Q;
    FT R2 = R * R;

    FT CR2 = FT(729) * r * r;
    FT CQ3 = FT(2916) * q * q * q;

    if (R == ZERO && Q == ZERO)
    {
        *x0 = - a / THREE;
        *x1 = - a / THREE;
        *x2 = - a / THREE;
        return 3;
    }
    else if (CR2 == CQ3)
    {
        /* this test is actually R2 == Q3, written in a form suitable
         for exact computation with integers */

        /* Due to finite precision some double roots may be missed, and
         considered to be a pair of complex roots z = x +/- epsilon i
         close to the real axis. */

        FT sqrtQ = sqrt (Q);

        if (R > ZERO)
        {
            *x0 = -TWO * sqrtQ  - a / THREE;
            *x1 = sqrtQ - a / THREE;
            *x2 = sqrtQ - a / THREE;
        }
        else
        {
            *x0 = - sqrtQ  - a / THREE;
            *x1 = - sqrtQ - a / THREE;
            *x2 = TWO * sqrtQ - a / THREE;
        }
        return 3;
    }
    else if (R2 < Q3)
    {
        FT sgnR = (R >= ZERO ? ONE : -ONE);
        FT ratio = sgnR * sqrt (R2 / Q3);
        FT theta = acos (ratio);
        FT norm = -TWO * sqrt (Q);
        *x0 = norm * cos (theta / THREE) - a / THREE;
        *x1 = norm * cos ((theta + TWO * pi<FT>()) / THREE) - a / THREE;
        *x2 = norm * cos ((theta - TWO * pi<FT>()) / THREE) - a / THREE;

        /* Sort *x0, *x1, *x2 into increasing order */

        if (*x0 > *x1)
            std::swap(*x0, *x1);

        if (*x1 > *x2)
        {
            std::swap(*x1, *x2);

            if (*x0 > *x1)
                std::swap(*x0, *x1);
        }

        return 3;
    }
    else
    {
        FT sgnR = (R >= ZERO ? ONE : -ONE);
        FT A = -sgnR * pow (fabs (R) + sqrt (R2 - Q3), ONE/THREE);
        FT B = Q / A;
        *x0 = A + B - a / THREE;
        return 1;
    }
}
} // namespace Math
} // namespace CS
