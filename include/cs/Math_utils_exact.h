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
#ifndef LIBCS_MATH_UTILS_EXACT_H
#define LIBCS_MATH_UTILS_EXACT_H

#if !defined(LIBCS_MATH_UTILS_EXACT_INTERNAL_FILE)
#if !defined(LIBCS_SPIN_KERNEL_INCLUDE_FILE)
#error Do not include this file directly, include Spin_kernel_3.h instead
#endif // !defined(LIBCS_SPIN_KERNEL_INCLUDE_FILE)
#endif // !defined(LIBCS_MATH_UTILS_EXACT_INTERNAL_FILE)

#include "Math_utils.h"

#include <CGAL/Gmpfr.h>
#include <mpfr.h>
#include <cassert>

namespace CS
{
// implement missing math operators for CGAL::Gmpfr
namespace Math
{
mpfr_rnd_t gmp_rounding_mode(std::float_round_style r);
CGAL::Gmpfr::Precision_type gmp_result_precision(const CGAL::Gmpfr &x);

CGAL::Gmpfr sin(const CGAL::Gmpfr &x);
CGAL::Gmpfr cos(const CGAL::Gmpfr &x);
CGAL::Gmpfr acos(const CGAL::Gmpfr &x);
CGAL::Gmpfr atan(const CGAL::Gmpfr &x);
CGAL::Gmpfr pow(const CGAL::Gmpfr &x, const CGAL::Gmpfr &y);
CGAL::Gmpfr atan2(const CGAL::Gmpfr &y, const CGAL::Gmpfr &x);
CGAL::Gmpfr sqrt(const CGAL::Gmpfr &x);
CGAL::Gmpfr fabs(const CGAL::Gmpfr &x);

template<>
CGAL::Gmpfr pi<CGAL::Gmpfr>();

void print(const CGAL::Gmpfr &x);
} // namespace Math
} // namespace CS

#include "Math_utils_exact.ipp"

#endif // LIBCS_MATH_UTILS_EXACT_H
