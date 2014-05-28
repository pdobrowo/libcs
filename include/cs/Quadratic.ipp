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
#include "Quadratic.h"

namespace CS
{
namespace Math
{
template<class FT_>
int solve_quadratic(const FT_ &a, const FT_ &b, const FT_ &c,
                    FT_ *x0, FT_ *x1)
{
    const FT_ ZERO = 0.0;
    const FT_ ONE  = 1.0;
    const FT_ FOUR = 4.0;
    const FT_ HALF = 0.5;

    FT_ disc = b * b - FOUR * a * c;

    if (a == ZERO) /* Handle linear case */
    {
        if (b == ZERO)
        {
            return 0;
        }
        else
        {
            *x0 = -c / b;
            return 1;
        }
    }

    if (disc > ZERO)
    {
        if (b == ZERO)
        {
            FT_ r = fabs (HALF * sqrt (disc) / a);
            *x0 = -r;
            *x1 =  r;
        }
        else
        {
            FT_ sgnb = (b > ZERO ? ONE : -ONE);
            FT_ temp = -HALF * (b + sgnb * sqrt (disc));
            FT_ r1 = temp / a;
            FT_ r2 = c / temp;

            if (r1 < r2)
            {
                *x0 = r1;
                *x1 = r2;
            }
            else
            {
                *x0 = r2;
                *x1 = r1;
            }
        }
        return 2;
    }
    else if (disc == ZERO)
    {
        *x0 = -HALF * b / a;
        *x1 = -HALF * b / a;
        return 2;
    }
    else
    {
        return 0;
    }
}
} // namespace Math
} // namespace CS
