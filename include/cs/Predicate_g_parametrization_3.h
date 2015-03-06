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
#ifndef LIBCS_PREDICATE_G_PARAMETRIZATION_3_H
#define LIBCS_PREDICATE_G_PARAMETRIZATION_3_H

#include "Spin_3.h"
#include "Matrix_44.h"
#include "Vector_4.h"
#include <CGAL/Kernel/global_functions.h>
#include <CGAL/Sqrt_extension.h>
#include <ostream>

namespace CS
{
//
// parametrization of a general predicate
//
template<class Kernel_>
class Predicate_g_parametrization_3
{
    typedef typename Kernel_::RT                RT;
    typedef typename Kernel_::Vector_3          Vector_3;
    typedef typename Kernel_::Predicate_g_3     Predicate_g_3;
    typedef typename Kernel_::Spin_quadric_3    Spin_quadric_3;
    typedef Matrix_44<RT>                       Matrix_RT;
    typedef CGAL::Sqrt_extension<RT, RT>        ERT;
    typedef Matrix_44<ERT>                      Matrix_ERT;
    typedef Vector_4<ERT>                       Vector_ERT;

public:
    typedef Kernel_                             Kernel;

    Predicate_g_parametrization_3(const Predicate_g_3 &g3);

    const Vector_3 &b() const;

    //  {u, v} in [0; 1]
    Spin_3<double> approx_evaluate(double u, double v) const;

private:
    Vector_3    m_b;
};
} // namespace CS

#include "Predicate_g_parametrization_3.ipp"

#endif // LIBCS_PREDICATE_G_PARAMETRIZATION_3_H
