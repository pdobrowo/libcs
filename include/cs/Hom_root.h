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
#ifndef LIBCS_HOM_ROOT_H
#define LIBCS_HOM_ROOT_H

#include "Benchmark.h"
#include "Isolate_polynomial_roots.h"
#include "Convert_polynomial.h"
#include "Polynomial_interval.h"
#include "Polynomial_sign_at.h"
#include <CGAL/Gmpz.h>
#include <CGAL/Sqrt_extension.h>
#include <CGAL/polynomial_utils.h>
#include <libqi.h>

namespace CS
{
// this class represents an homogeneous root of a homogeneous polynomial
// the value is represented as a isolated root of a univariate polynomial
// or a rational number in some cases
template<class Kernel_>
class Hom_root
{
    typedef typename Kernel_::Algebraic_kernel_with_sqrt            Algebraic_kernel_with_sqrt;
    typedef typename Algebraic_kernel_with_sqrt::Algebraic_real_1   Algebraic_real_1;

public:
    Hom_root(bool infinity,
             const Algebraic_real_1 &parameter);

    bool is_infinite() const;
    bool is_finite() const;

    const Algebraic_real_1 &root() const;

private:
    bool                m_infinity;
    Algebraic_real_1    m_root;
};

template<class Kernel_, typename OutputIterator_>
void extract_isolated_hom_roots(
        const typename Kernel_::Hom_polynomial_with_sqrt &hom_polynomial_,
        OutputIterator_ output);

bigint to_bigint(const CGAL::Gmpz &g);

// return sign of hom polynomial at value at given by hom root
template<class Kernel_>
CGAL::Sign hom_polynomial_with_sqrt_sign_at(
        const typename Kernel_::Hom_polynomial_with_sqrt &hom_polynomial,
        const typename Kernel_::Hom_root &hom_root);

} // namespace CS

#include "Hom_root.ipp"

#endif // LIBCS_HOM_ROOT_H
