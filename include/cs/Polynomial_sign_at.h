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
#ifndef LIBCS_POLYNOMIAL_SIGN_AT_H
#define LIBCS_POLYNOMIAL_SIGN_AT_H

#include <CGAL/Uncertain.h>
#include <CGAL/number_utils.h>
#include <CGAL/polynomial_utils.h>
#include "Benchmark.h"

namespace CS
{
/*
 * New algorithm for polynomial interval sign determination
 *
 * It is faster than the algorithm proposed by Hemmer et. al.
 * The idea is to assume that the plynomial is locally mototonic.
 * By using polynomial sign at point, and polynomial differential it is possible
 * to deduce what the sign of whole interval is.
 */
template<class Polynomial_1_, class RT_>
CGAL::Uncertain<CGAL::Sign> polynomial_interval_sign_at(
        const std::vector<Polynomial_1_> &polynomials,
        const size_t &level,
        const RT_ &low, const RT_ &high);

template<class Polynomial_1_, class Algebraic_real_1_>
CGAL::Uncertain<CGAL::Sign> polynomial_sign_at(
        const std::vector<Polynomial_1_> &polynomials,
        const Algebraic_real_1_ &r);

template<class Polynomial_1_>
void polynomial_sign_at_prepare(
        const Polynomial_1_ &polynomial,
        std::vector<Polynomial_1_> &polynomials);

} // namespace CS

#include "Polynomial_sign_at.ipp"

#endif // LIBCS_POLYNOMIAL_SIGN_AT_H
