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
#ifndef LIBCS_SPIN_QSIP_3_H
#define LIBCS_SPIN_QSIP_3_H

#include "Benchmark.h"
#include "Hom_root.h"
#include <CGAL/Kernel/global_functions.h>
#include <CGAL/Sqrt_extension.h>
#include <string>
#include <vector>
#include <cassert>

namespace CS
{
// Spin_qsip_3:
//
//    P - oriented Spin_quadric_3
//    C - Spin_qsic_3
//
//    qsip = c(g_i)
//    g_i - zeros on curve C
template<class Kernel_>
class Spin_qsip_3
{
    typedef typename Kernel_::RT                RT;

    typedef typename Kernel_::Qsic_component    Qsic_component;
    typedef typename Kernel_::Qsic_curve        Qsic_curve;
    typedef typename Kernel_::Qsic_surface      Qsic_surface;

    typedef typename Kernel_::Spin_quadric_3    Spin_quadric_3;
    typedef typename Kernel_::Spin_qsic_3       Spin_qsic_3;

    typedef typename Kernel_::Spin_qsip_point   Spin_qsip_point;
    typedef typename Kernel_::Hom_root          Hom_root;

    typedef typename Kernel_::Hom_polynomial        Hom_polynomial;
    typedef typename Kernel_::Hom_hom_polynomial    Hom_hom_polynomial;

    typedef std::vector<Spin_qsip_point>                    Spin_qsip_point_list;
    typedef typename Spin_qsip_point_list::const_iterator   Point_const_iterator;

    typedef typename Kernel_::Algebraic_kernel_with_sqrt    Algebraic_kernel_with_sqrt;

    typedef typename Algebraic_kernel_with_sqrt::Polynomial_1       Polynomial_1;
    typedef typename Algebraic_kernel_with_sqrt::Algebraic_real_1   Algebraic_real_1;
    typedef typename Algebraic_kernel_with_sqrt::Bound              Bound;
    typedef typename Algebraic_kernel_with_sqrt::Multiplicity_type  Multiplicity_type;
    typedef typename Algebraic_kernel_with_sqrt::Solve_1            Solve_1;

    typedef typename Kernel_::Sqrt_extension            Sqrt_extension;
    typedef typename Kernel_::Hom_polynomial_with_sqrt  Hom_polynomial_with_sqrt;

    typedef Spin_qsip_3                                                 Self;

public:
    typedef Kernel_                                                     Kernel;

    Spin_qsip_3(const Spin_quadric_3 &q, const Spin_qsic_3 &c);

    size_t                      size_of_points() const;

    Point_const_iterator        points_begin() const;
    Point_const_iterator        points_end() const;

    const Spin_qsip_point &     point_at(size_t index) const;

private:
    Spin_qsip_point_list        m_points;

    // generic case
    void                        intersect_smooth_qsic(
                                    const Spin_quadric_3 &quadric,
                                    const Spin_qsic_3 &qsic);

    // special cases
    void                        intersect_rational_component(const Spin_quadric_3 &quadric,
                                    const Qsic_component &component);

    // extract intersections on a curve as a zeroes of a characteristic polynomial
    void                        populate_component_intersections(
                                    const Hom_polynomial_with_sqrt &characteristic,
                                    const Qsic_component &component);

    // utilities
    void                        mul3(
                                    Hom_polynomial_with_sqrt &out,
                                    const Hom_polynomial_with_sqrt &a,
                                    const Hom_polynomial_with_sqrt &b,
                                    const Hom_polynomial_with_sqrt &c) const;

    void                        add10(
                                    Hom_polynomial_with_sqrt &out,
                                    const Hom_polynomial_with_sqrt a[10]) const;

    void                        calc_implicit_qsic_equation(
                                    const bigint_matrix &quadric,       // first quadric
                                    const Qsic_surface &s1,             // bilinear parametrization of second quadric: first component of Q[sqrt(xi)]
                                    const Qsic_surface &s2,             // bilinear parametrization of second quadric: second component of Q[sqrt(xi)]
                                    const bigint &delta,                // bilinear parametrization of second quadric: xi
                                    Hom_hom_polynomial &outImplicitS1,  // out hom_hom polynomial: first component of Q[sqrt(xi)]
                                    Hom_hom_polynomial &outImplicitS2); // out hom_hom polynomial: second component of Q[sqrt(xi)]

    Hom_polynomial_with_sqrt    get_coefficient(
                                    int index,
                                    const Hom_hom_polynomial &a,
                                    const Hom_hom_polynomial &b,
                                    const RT &root);

    Hom_polynomial_with_sqrt    extend_hom_polynomial(
                                    const Hom_polynomial &a0);

    Hom_polynomial_with_sqrt    extend_hom_polynomial(
                                    const Hom_polynomial &a0,
                                    const Hom_polynomial &a1,
                                    const RT &root);
};
} // namespace CS

#include "Spin_qsip_3.ipp"

#endif // LIBCS_SPIN_QSIP_3_H
