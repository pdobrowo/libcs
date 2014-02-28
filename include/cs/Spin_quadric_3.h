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
#ifndef LIBCS_SPIN_QUADRIC_3_H
#define LIBCS_SPIN_QUADRIC_3_H

#include <CGAL/Kernel/global_functions.h>
#include "Spin_3.h"
#include "Predicate_h_3.h"
#include "Predicate_s_3.h"
#include "Predicate_g_3.h"
#include <ostream>
#include <cassert>
#include <sstream>
#include <string>

namespace CS
{
template<class Kernel_>
class Spin_quadric_3;

template<class Kernel_>
bool operator <(const Spin_quadric_3<Kernel_> &lhs, const Spin_quadric_3<Kernel_> &rhs);

// Spin_quadric_3:
//    F(s) = a11 * s12^2 + a22 * s23^2 + a33 * s31^2 + a44 * s0^2 +
//           2 * a12 * s12 * s23 + 2 * a13 * s12 * s31 + 2 * a14 * s12 * s0 +
//           2 * a23 * s23 * s31 + 2 * a24 * s23 * s0 + 2 * a34 * s31 * s0
//
//         | a11 a12 a13 a14 |
//    Qs = | a12 a22 a23 a24 |
//         | a13 a23 a33 a34 |
//         | a14 a24 a34 a44 |
//
//    s = [s12; s23, s31; s0]^T
//
//    F(s) = s^T Qs s
//
//    s12^2 + s23^2 + s31^2 + s0^2 = 1
//
template<class Kernel_>
class Spin_quadric_3
{
    typedef typename Kernel_::Vector_3      Vector_3;

    typedef Spin_quadric_3                  Self;

    void construct(const Predicate_g_3<Kernel_> &g3);

public:    
    typedef typename Kernel_::RT            RT;
    //typedef typename Kernel_::FT          FT;

    typedef Kernel_                      R;
    typedef typename Kernel_::Matrix     Matrix;

    Spin_quadric_3();

    Spin_quadric_3(const Predicate_h_3<Kernel_> &h3);

    Spin_quadric_3(const Predicate_s_3<Kernel_> &s3);

    Spin_quadric_3(const Predicate_g_3<Kernel_> &g3);

    const RT &      a11() const;
    const RT &      a22() const;
    const RT &      a33() const;
    const RT &      a44() const;
    const RT &      a12() const;
    const RT &      a13() const;
    const RT &      a14() const;
    const RT &      a23() const;
    const RT &      a24() const;
    const RT &      a34() const;

    Matrix          matrix() const;

    // four dimensional ellipsoid according to ellipsoid theorem
    Matrix          ellipsoid_matrix() const;

    std::string     to_string() const;

    // maximum of coefficients absolute values
    RT              max_abs_coefficient() const;

    // scale up the predicate
    // resulting predicate is equivalent
    // requires: scale > 0
    void            scale(const RT &s);

    void            inverse();
    Spin_quadric_3  inversed() const;

    template<class NT>
    NT              evaluate(const CS::Spin_3<NT> &spin) const;

    template<class NT>
    CGAL::Sign      evaluate_sign(const Spin_3<NT> &spin) const;

    friend bool operator < <>(const Spin_quadric_3<Kernel_> &lhs, const Spin_quadric_3<Kernel_> &rhs);

private:
    RT  m_a11, m_a22, m_a33, m_a44, m_a12, m_a13, m_a14, m_a23, m_a24, m_a34;
};
} // namespace CS

#include "Spin_quadric_3.ipp"

#endif // LIBCS_SPIN_QUADRIC_3_H
