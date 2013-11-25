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
#ifndef LIBCS_POLYNOMIAL_INTERVAL_H
#define LIBCS_POLYNOMIAL_INTERVAL_H

#include "Benchmark.h"

#if 0 // BAD CODE

// return sign of hom polynomial at value at given by given root
template<class Polynomial_1, class Algebraic_real_1>
CGAL::Gmpfi polynomial_interval_value(
        const Polynomial_1 &polynomial,
        const Algebraic_real_1 &real)
{
    CS_BENCHMARK_POINT();

    CGAL::Gmpfi low(CGAL::Gmpfi(real.low().numerator()) /
                    CGAL::Gmpfi(real.low().denominator()));

    CGAL::Gmpfi high(CGAL::Gmpfi(real.high().numerator()) /
                     CGAL::Gmpfi(real.high().denominator()));

    CGAL::Gmpfi interval(low.left_mpfr(), high.right_mpfr());

    CGAL::Gmpfi value = polynomial.evaluate(interval);

    return value;

//  This is slow:
//  return polynomial.evaluate(CGAL::convert_to_bfi(real));
}

#endif

#include "Polynomial_interval.ipp"

#endif // LIBCS_POLYNOMIAL_INTERVAL_H
