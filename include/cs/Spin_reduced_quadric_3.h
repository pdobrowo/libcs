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
#ifndef LIBCS_SPIN_REDUCED_QUADRIC_3_H
#define LIBCS_SPIN_REDUCED_QUADRIC_3_H

#include <CGAL/Kernel/global_functions.h>
#include "Spin_quadric_3.h"

namespace CS
{
// Spin_reduced_quadric_3:
//    reduction of a spin quadric to a quartic in R^3
//
template<class Kernel_>
class Spin_reduced_quadric_3
{
    typedef Spin_reduced_quadric_3              Self;

public:    
    typedef typename Kernel_::FT                FT;
    typedef typename Kernel_::Spin_quadric_3    Spin_quadric_3;

    typedef Kernel_                             R;

    Spin_reduced_quadric_3();

    Spin_reduced_quadric_3(FT a11, FT a22, FT a33, FT a44, FT a12, FT a13, FT a14, FT a23, FT a24, FT a34);

    Spin_reduced_quadric_3(const Spin_quadric_3 &spin_quadric);

    //const RT &      a11() const;

    void            gradient(FT s12, FT s23, FT s31, FT &g12, FT &g23, FT &g31) const;
    FT              evaluate(FT s12, FT s23, FT s31) const;
    void            squared_s0_and_sign(FT s12, FT s23, FT s31, FT &squared_s0, CGAL::Sign &sign_s0) const;

private:
    // 22 coeffs
    FT m_c1233, m_c1223, m_c1123;
    FT m_c1122, m_c1133, m_c2233;
    FT m_c1113, m_c1112, m_c2223, m_c2221, m_c3332, m_c3331;
    FT m_c1111, m_c2222, m_c3333;
    FT m_c12, m_c13, m_c23;
    FT m_c11, m_c22, m_c33;
    FT m_c;

    void construct(FT a11, FT a22, FT a33, FT a44, FT a12, FT a13, FT a14, FT a23, FT a24, FT a34);
};
} // namespace CS

#include "Spin_reduced_quadric_3.ipp"

#endif // LIBCS_SPIN_REDUCED_QUADRIC_3_H
