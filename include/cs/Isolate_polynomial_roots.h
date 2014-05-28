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
#ifndef LIBCS_ISOLATE_POLYNOMIAL_ROOTS_H
#define LIBCS_ISOLATE_POLYNOMIAL_ROOTS_H

#include "Benchmark.h"
#include "Convert_polynomial.h"
#include "Uspensky.h"
#include "Polynomial_interval.h"
#include "Contfrac.h"
#include <boost/scoped_array.hpp>

// isolation algorithm
//#define isolate_polynomial_roots  isolate_polynomial_roots_cgal
//#define isolate_polynomial_roots  isolate_polynomial_roots_hanrot
#define isolate_polynomial_roots  isolate_polynomial_roots_contfrac

namespace CS
{
template<class Kernel_, typename OutputIterator_>
void isolate_polynomial_roots_cgal(
        const typename Kernel_::Algebraic_kernel_with_sqrt::Polynomial_1 &polynomial,
        OutputIterator_ output);

template<class Kernel_>
typename Kernel_::Algebraic_kernel_with_sqrt::Algebraic_real_1 interval_to_algebraic_real(
        const typename Kernel_::Algebraic_kernel_with_sqrt::Polynomial_1 &polynomial,
        const interval &z);

template<class Kernel_, typename OutputIterator_>
void isolate_polynomial_roots_contfrac(
        const typename Kernel_::Algebraic_kernel_with_sqrt::Polynomial_1 &polynomial,
        OutputIterator_ output);
} // namespace CS

#include "Isolate_polynomial_roots.ipp"

#endif // LIBCS_ISOLATE_POLYNOMIAL_ROOTS_H
