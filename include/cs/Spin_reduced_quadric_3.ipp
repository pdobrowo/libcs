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
#include "Spin_reduced_quadric_3.h"

namespace CS
{
template<class Kernel_>
Spin_reduced_quadric_3<Kernel_>::Spin_reduced_quadric_3()
    : m_c1233(FT(0)),
      m_c1223(FT(0)),
      m_c1123(FT(0)),
      m_c1122(FT(0)),
      m_c1133(FT(0)),
      m_c2233(FT(0)),
      m_c1113(FT(0)),
      m_c1112(FT(0)),
      m_c2223(FT(0)),
      m_c2221(FT(0)),
      m_c3332(FT(0)),
      m_c3331(FT(0)),
      m_c1111(FT(0)),
      m_c2222(FT(0)),
      m_c3333(FT(0)),
      m_c12(FT(0)),
      m_c13(FT(0)),
      m_c23(FT(0)),
      m_c11(FT(0)),
      m_c22(FT(0)),
      m_c33(FT(0)),
      m_c(FT(0))
{
}

template<class Kernel_>
Spin_reduced_quadric_3<Kernel_>::Spin_reduced_quadric_3(const Spin_quadric_3 &spin_quadric)
{
    construct(spin_quadric.a11(), spin_quadric.a22(), spin_quadric.a33(), spin_quadric.a44(),
              spin_quadric.a12(), spin_quadric.a13(), spin_quadric.a14(), spin_quadric.a23(), spin_quadric.a24(), spin_quadric.a34());
}

template<class Kernel_>
void Spin_reduced_quadric_3<Kernel_>::construct(FT a11, FT a22, FT a33, FT a44, FT a12, FT a13, FT a14, FT a23, FT a24, FT a34)
{
    // coeffs
    m_c1233 = 8 * a13 * a23 + 4 * a12 * (a33 - a44) + 8 * a14 * a24;
    m_c1223 = 8 * a12 * a23 + 4 * a13 * (a22 - a44) + 8 * a14 * a34;
    m_c1123 = 8 * a12 * a13 + 4 * a23 * (a11 - a44) + 8 * a24 * a34;
    m_c1122 = 4 * a12 * a12 + 2 * (a11 - a44) * (a22 - a44) + 4 * a24 * a24 + 4 * a14 * a14;
    m_c1133 = 4 * a13 * a13 + 2 * (a11 - a44) * (a33 - a44) + 4 * a34 * a34 + 4 * a14 * a14;
    m_c2233 = 4 * a23 * a23 + 2 * (a22 - a44) * (a33 - a44) + 4 * a34 * a34 + 4 * a24 * a24;
    m_c1113 = 4 * a13 * (a11 - a44) + 8 * a14 * a34;
    m_c1112 = 4 * a12 * (a11 - a44) + 8 * a14 * a24;
    m_c2223 = 4 * a23 * (a22 - a44) + 8 * a24 * a34;
    m_c2221 = 4 * a12 * (a22 - a44) + 8 * a14 * a24;
    m_c3332 = 4 * a23 * (a33 - a44) + 8 * a24 * a34;
    m_c3331 = 4 * a13 * (a33 - a44) + 8 * a14 * a34;
    m_c1111 = (a11 - a44) * (a11 - a44) + 4 * a14 * a14;
    m_c2222 = (a22 - a44) * (a22 - a44) + 4 * a24 * a24;
    m_c3333 = (a33 - a44) * (a33 - a44) + 4 * a34 * a34;
    m_c12 = 4 * a12 * a44 - 8 * a14 * a24;
    m_c13 = 4 * a13 * a44 - 8 * a14 * a34;
    m_c23 = 4 * a23 * a44 - 8 * a24 * a34;
    m_c11 = 2 * (a11 - a44) * a44 - 4 * a14 * a14;
    m_c22 = 2 * (a22 - a44) * a44 - 4 * a24 * a24;
    m_c33 = 2 * (a33 - a44) * a44 - 4 * a34 * a34;
    m_c = a44 * a44;
}

template<class Kernel_>
void Spin_reduced_quadric_3<Kernel_>::gradient(FT s12, FT s23, FT s31, FT &g12, FT &g23, FT &g31) const
{
    g12 = m_c1233     * s23 * s31 * s31 +
          m_c1223     * s23 * s23 * s31 +
          m_c1123 * 2 * s12 * s23 * s31 +
          m_c1122 * 2 * s12 * s23 * s23 +
          m_c1133 * 2 * s12 * s31 * s31 +
          m_c1113 * 3 * s12 * s12 * s31 +
          m_c1112 * 2 * s12 * s12 * s23 +
          m_c2221     * s23 * s23 * s23 +
          m_c3331     * s31 * s31 * s31 +
          m_c1111 * 4 * s12 * s12 * s12 +
          m_c12       * s23 +
          m_c13       * s31 +
          m_c11   * 2 * s12;

    g23 = m_c1233     * s12 * s31 * s31 +
          m_c1223 * 2 * s12 * s23 * s31 +
          m_c1123     * s12 * s12 * s31 +
          m_c1122 * 2 * s12 * s12 * s23 +
          m_c2233 * 2 * s23 * s31 * s31 +
          m_c1112     * s12 * s12 * s12 +
          m_c2223 * 3 * s23 * s23 * s31 +
          m_c2221 * 3 * s23 * s23 * s12 +
          m_c3332     * s31 * s31 * s31 +
          m_c2222 * 4 * s23 * s23 * s23 +
          m_c12       * s12 +
          m_c23       * s31 +
          m_c22   * 2 * s23;

    g31 = m_c1233 * 2 * s12 * s23 * s31 +
          m_c1223     * s12 * s23 * s23 +
          m_c1123     * s12 * s12 * s23 +
          m_c1133 * 2 * s12 * s12 * s31 +
          m_c2233 * 2 * s23 * s23 * s31 +
          m_c1113     * s12 * s12 * s12 +
          m_c2223     * s23 * s23 * s23 +
          m_c3332 * 3 * s23 * s31 * s31 +
          m_c3331 * 3 * s12 * s31 * s31 +
          m_c3333 * 4 * s31 * s31 * s31 +
          m_c13       * s12 +
          m_c23       * s23 +
          m_c33   * 2 * s31;
}

template<class Kernel_>
typename Spin_reduced_quadric_3<Kernel_>::FT Spin_reduced_quadric_3<Kernel_>::evaluate(FT s12, FT s23, FT s31) const
{
    return m_c1233 * s12 * s23 * s31 * s31 +
           m_c1223 * s12 * s23 * s23 * s31 +
           m_c1123 * s12 * s12 * s23 * s31 +
           m_c1122 * s12 * s12 * s23 * s23 +
           m_c1133 * s12 * s12 * s31 * s31 +
           m_c2233 * s23 * s23 * s31 * s31 +
           m_c1113 * s12 * s12 * s12 * s31 +
           m_c1112 * s12 * s12 * s12 * s23 +
           m_c2223 * s23 * s23 * s23 * s31 +
           m_c2221 * s23 * s23 * s23 * s12 +
           m_c3332 * s31 * s31 * s31 * s23 +
           m_c3331 * s31 * s31 * s31 * s12 +
           m_c1111 * s12 * s12 * s12 * s12 +
           m_c2222 * s23 * s23 * s23 * s23 +
           m_c3333 * s31 * s31 * s31 * s31 +
           m_c12   * s12 * s23 +
           m_c13   * s12 * s31 +
           m_c23   * s23 * s31 +
           m_c11   * s12 * s12 +
           m_c22   * s23 * s23 +
           m_c33   * s31 * s31 +
           m_c;
}

template<class Kernel_>
void Spin_reduced_quadric_3<Kernel_>::squared_s0_and_sign(FT s12, FT s23, FT s31, FT &squared_s0, CGAL::Sign &sign_s0) const
{
    FT s0abs = std::sqrt(1 - s12 * s12 - s23 * s23 - s31 * s31);

    ; //
    (void)squared_s0;
    (void)sign_s0;

    return s0abs;
}
} // namespace CS
