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
#ifndef LIBCS_CONVERT_POLYNOMIAL_H
#define LIBCS_CONVERT_POLYNOMIAL_H

#include <CGAL/Gmpz.h>
#include <CGAL/Sqrt_extension.h>
#include <malloc.h>
#include <gmp.h>
#include <assert.h>

template<class Kernel_>
typename Kernel_::Algebraic_kernel_with_sqrt::Polynomial_1 convert_polynomial(
        const typename Kernel_::Hom_polynomial_with_sqrt &hom_polynomial);

template<class Kernel_>
typename Kernel_::Algebraic_kernel_with_sqrt::Polynomial_1 extend_simple_polynomial(
        const typename Kernel_::Algebraic_kernel::Polynomial_1 &polynomial);

template<class Kernel_>
void apart_polynomial(
        const typename Kernel_::Algebraic_kernel_with_sqrt::Polynomial_1 polynomial,
        typename Kernel_::Algebraic_kernel::Polynomial_1 &left_polynomial,
        CGAL::Gmpz &d,
        typename Kernel_::Algebraic_kernel::Polynomial_1 &right_polynomial);

template<class Kernel_>
bool is_polynomial_extended(
        const typename Kernel_::Algebraic_kernel_with_sqrt::Polynomial_1 polynomial);

template<class Kernel_>
typename Kernel_::Algebraic_kernel::Polynomial_1 remove_polynomial_extension(
        const typename Kernel_::Algebraic_kernel_with_sqrt::Polynomial_1 polynomial);

template<class Kernel_>
void raw_polynomial_alloc(
        const typename Kernel_::Algebraic_kernel::Polynomial_1 polynomial,
        mpz_t *&coeffs,
        unsigned long &degree);

template<class Kernel_>
void raw_polynomial_free(
        mpz_t *&coeffs,
        unsigned long &degree);

#include "Convert_polynomial.ipp"

#endif // LIBCS_CONVERT_POLYNOMIAL_H
