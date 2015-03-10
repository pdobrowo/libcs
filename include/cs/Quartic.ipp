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
#include "Quartic.h"

namespace CS
{
namespace Math
{
template<class FT_>
int solve_quartic(const FT_ &a, const FT_ &b, const FT_ &c, const FT_ &d,
                  FT_ *x0, FT_ *x1, FT_ *x2, FT_ *x3)
{
    const FT_ ZERO  = 0.0;
    const FT_ HALF  = 0.5;
    const FT_ ONE   = 1.0;
    const FT_ TWO   = 2.0;
    const FT_ THREE = 3.0;
    const FT_ FOUR  = 4.0;
    const FT_ EIGHT = 8.0;
    const FT_ NINE  = 9.0;

    /*
     * This code is based on a simplification of
     * the algorithm from zsolve_quartic.c for real roots
     */
    FT_ u[3];
    FT_ aa, pp, qq, rr, rc, sc, tc;
    FT_ w1r, w1i, w2r, w2i, w3r;
    FT_ v[3], v1, v2, arg, theta;
    FT_ disc, h;
    int k1 = 0, k2 = 0, mt;
    FT_ zarr[4];

    /* Deal easily with the cases where the quartic is degenerate. The
     * ordering of solutions is done explicitly. */
    if (ZERO == b && ZERO == c)
    {
        if (ZERO == d)
        {
            if (a > ZERO)
            {
                *x0 = -a;
                *x1 = ZERO;
                *x2 = ZERO;
                *x3 = ZERO;
            }
            else
            {
                *x0 = ZERO;
                *x1 = ZERO;
                *x2 = ZERO;
                *x3 = -a;
            }
            return 4;
        }
        else if (ZERO == a)
        {
            if (d > ZERO)
            {
                return 0;
            }
            else
            {
                *x1 = sqrt (sqrt (-d));
                *x0 = -(*x1);
                return 2;
            }
        }
    }

    if (ZERO == c && ZERO == d)
    {
        *x0 = ZERO;
        *x1 = ZERO;

        if (solve_quadratic(ONE, a, b, x2, x3) == 0)
        {
            mt = 3;
        }
        else
        {
            mt = 1;
        }
    }
    else
    {
        /* For non-degenerate solutions, proceed by constructing and
         * solving the resolvent cubic */
        aa = a * a;
        pp = b - (THREE/EIGHT) * aa;
        qq = c - (ONE/TWO) * a * (b - (ONE/FOUR) * aa);
        rr = d - (ONE/FOUR) * (a * c - (ONE/FOUR) * aa * (b - (THREE/FT_(16)) * aa));
        rc = (ONE/TWO) * pp;
        sc = (ONE/FOUR) * ((ONE/FOUR) * pp * pp - rr);
        tc = -((ONE/EIGHT) * qq * (ONE/EIGHT) * qq);

        /* This code solves the resolvent cubic in a convenient fashion
         * for this implementation of the quartic. If there are three real
         * roots, then they are placed directly into u[].  If two are
         * complex, then the real root is put into u[0] and the real
         * and imaginary part of the complex roots are placed into
         * u[1] and u[2], respectively. Additionally, this
         * calculates the discriminant of the cubic and puts it into the
         * variable disc. */
        {
            FT_ qcub = (rc * rc - THREE * sc);
            FT_ rcub = (TWO * rc * rc * rc - NINE * rc * sc + FT_(27.0) * tc);

            FT_ Q = qcub / NINE;
            FT_ R = rcub / FT_(54.0);

            FT_ Q3 = Q * Q * Q;
            FT_ R2 = R * R;

            FT_ CR2 = FT_(729.0) * rcub * rcub;
            FT_ CQ3 = FT_(2916.0) * qcub * qcub * qcub;

            disc = (CR2 - CQ3) / FT_(2125764.0);

            if (ZERO == R && ZERO == Q)
            {
                u[0] = -rc / THREE;
                u[1] = -rc / THREE;
                u[2] = -rc / THREE;
            }
            else if (CR2 == CQ3)
            {
                FT_ sqrtQ = sqrt (Q);
                if (R > ZERO)
                {
                    u[0] = -TWO * sqrtQ - rc / THREE;
                    u[1] = sqrtQ - rc / THREE;
                    u[2] = sqrtQ - rc / THREE;
                }
                else
                {
                    u[0] = -sqrtQ - rc / THREE;
                    u[1] = -sqrtQ - rc / THREE;
                    u[2] = TWO * sqrtQ - rc / THREE;
                }
            }
            else if (CR2 < CQ3)
            {
                FT_ sqrtQ = sqrt (Q);
                FT_ sqrtQ3 = sqrtQ * sqrtQ * sqrtQ;
                FT_ theta = acos (R / sqrtQ3);
                if (R / sqrtQ3 >= ONE) theta = ZERO;
                {
                    FT_ norm = -TWO * sqrtQ;

                    u[0] = norm * cos (theta / THREE) - rc / THREE;
                    u[1] = norm * cos ((theta + TWO * pi<FT_>()) / THREE) - rc / THREE;
                    u[2] = norm * cos ((theta - TWO * pi<FT_>()) / THREE) - rc / THREE;
                }
            }
            else
            {
                FT_ sgnR = (R >= ZERO ? ONE : -ONE);
                FT_ modR = fabs (R);
                FT_ sqrt_disc = sqrt (R2 - Q3);
                FT_ A = -sgnR * pow (modR + sqrt_disc, ONE / THREE);
                FT_ B = Q / A;
                FT_ mod_diffAB = fabs (A - B);

                u[0] = A + B - rc / THREE;
                u[1] = -HALF * (A + B) - rc / THREE;
                u[2] = -(sqrt (THREE) / TWO) * mod_diffAB;
            }
        }
        /* End of solution to resolvent cubic */

        /* Combine the square roots of the roots of the cubic
         * resolvent appropriately. Also, calculate 'mt' which
         * designates the nature of the roots:
         * mt=1 : 4 real roots (disc == 0)
         * mt=2 : 0 real roots (disc < 0)
         * mt=3 : 2 real roots (disc > 0)
         */

        if (ZERO == disc)
            u[2] = u[1];

        if (ZERO >= disc)
        {
            mt = 2;

            /* One would think that we could return 0 here and exit,
             * since mt=2. However, this assignment is temporary and
             * changes to mt=1 under certain conditions below.
             */

            v[0] = fabs (u[0]);
            v[1] = fabs (u[1]);
            v[2] = fabs (u[2]);

            v1 = std::max(std::max(v[0], v[1]), v[2]);
            /* Work out which two roots have the largest moduli */
            k1 = 0; k2 = 0;
            if (v1 == v[0])
            {
                k1 = 0;
                v2 = std::max (v[1], v[2]);
            }
            else if (v1 == v[1])
            {
                k1 = 1;
                v2 = std::max (v[0], v[2]);
            }
            else
            {
                k1 = 2;
                v2 = std::max (v[0], v[1]);
            }

            if (v2 == v[0])
            {
                k2 = 0;
            }
            else if (v2 == v[1])
            {
                k2 = 1;
            }
            else
            {
                k2 = 2;
            }

            if (ZERO <= u[k1])
            {
                w1r=sqrt(u[k1]);
                w1i=ZERO;
            }
            else
            {
                w1r=ZERO;
                w1i=sqrt(-u[k1]);
            }
            if (ZERO <= u[k2])
            {
                w2r=sqrt(u[k2]);
                w2i=ZERO;
            }
            else
            {
                w2r=ZERO;
                w2i=sqrt(-u[k2]);
            }
        }
        else
        {
            mt = 3;

            if (ZERO == u[1] && ZERO == u[2])
            {
                arg = ZERO;
            }
            else
            {
                arg = sqrt(sqrt(u[1] * u[1] + u[2] * u[2]));
            }
            theta = atan2(u[2], u[1]);

            w1r = arg * cos(theta / TWO);
            w1i = arg * sin(theta / TWO);
            w2r = w1r;
            w2i = -w1i;
        }

        /* Solve the quadratic to obtain the roots to the quartic */
        w3r = qq / EIGHT * (w1i * w2i - w1r * w2r) /
                (w1i * w1i + w1r * w1r) / (w2i * w2i + w2r * w2r);
        h = a / FOUR;

        zarr[0] = w1r + w2r + w3r - h;
        zarr[1] = -w1r - w2r + w3r - h;
        zarr[2] = -w1r + w2r - w3r - h;
        zarr[3] = w1r - w2r - w3r - h;

        /* Arrange the roots into the variables z0, z1, z2, z3 */
        if (2 == mt)
        {
            if (u[k1] >= 0 && u[k2] >= 0)
            {
                mt = 1;
                *x0 = zarr[0];
                *x1 = zarr[1];
                *x2 = zarr[2];
                *x3 = zarr[3];
            }
            else
            {
                return 0;
            }
        }
        else
        {
            *x0 = zarr[0];
            *x1 = zarr[1];
        }
    }

    /* Sort the roots as usual */
    if (1 == mt)
    {
        /* Roots are all real, sort them by the real part */
        if (*x0 > *x1)
            std::swap (*x0, *x1);
        if (*x0 > *x2)
            std::swap (*x0, *x2);
        if (*x0 > *x3)
            std::swap (*x0, *x3);

        if (*x1 > *x2)
            std::swap (*x1, *x2);
        if (*x2 > *x3)
        {
            std::swap (*x2, *x3);
            if (*x1 > *x2)
                std::swap (*x1, *x2);
        }
        return 4;
    }
    else
    {
        /* 2 real roots */
        if (*x0 > *x1)
            std::swap (*x0, *x1);
    }

    return 2;
}
} // namespace Math
} // namespace CS
