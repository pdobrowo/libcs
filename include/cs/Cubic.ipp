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
#include "Cubic.h"

namespace CS
{
namespace Math
{
template<class FT_>
int solve_cubic(const FT_ &a, const FT_ &b, const FT_ &c,
                FT_ *x0, FT_ *x1, FT_ *x2)
{
    const FT_ ZERO  = 0.0;
    const FT_ ONE   = 1.0;
    const FT_ TWO   = 2.0;
    const FT_ THREE = 3.0;
    const FT_ NINE  = 9.0;

    FT_ q = (a * a - THREE * b);
    FT_ r = (TWO * a * a * a - NINE * a * b + FT_(27) * c);

    FT_ Q = q / NINE;
    FT_ R = r / FT_(54);

    FT_ Q3 = Q * Q * Q;
    FT_ R2 = R * R;

    FT_ CR2 = FT_(729) * r * r;
    FT_ CQ3 = FT_(2916) * q * q * q;

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

        FT_ sqrtQ = sqrt (Q);

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
        FT_ sgnR = (R >= ZERO ? ONE : -ONE);
        FT_ ratio = sgnR * sqrt (R2 / Q3);
        FT_ theta = acos (ratio);
        FT_ norm = -TWO * sqrt (Q);
        *x0 = norm * cos (theta / THREE) - a / THREE;
        *x1 = norm * cos ((theta + TWO * pi<FT_>()) / THREE) - a / THREE;
        *x2 = norm * cos ((theta - TWO * pi<FT_>()) / THREE) - a / THREE;

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
        FT_ sgnR = (R >= ZERO ? ONE : -ONE);
        FT_ A = -sgnR * pow (fabs (R) + sqrt (R2 - Q3), ONE/THREE);
        FT_ B = Q / A;
        *x0 = A + B - a / THREE;
        return 1;
    }
}
} // namespace Math
} // namespace CS
