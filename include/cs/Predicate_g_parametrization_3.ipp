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
#include "Predicate_g_parametrization_3.h"
#include <algorithm>

namespace CS
{
template<class Kernel_>
Predicate_g_parametrization_3<Kernel_>::Predicate_g_parametrization_3(const Predicate_g_3 &g3)
{
    Vector_3 k = g3.k();
    Vector_3 l = g3.l();
    Vector_3 a = g3.a();
    Vector_3 b = g3.b();
    RT c = g3.c();

    if (c < RT(0))
    {
        std::swap(k, l);
        c = -c;
    }

    RT p, d;

    RT i = CGAL::cross_product(k, l).squared_length();
    RT j = CGAL::cross_product(a, b).squared_length();

    if (i == RT(0))
    {
        p = j;
        d = RT(0);
    }
    else if (j == RT(0))
    {
        p = i;
        d = RT(0);
    }
    else
    {
        if (i < j)
            p = i;
        else
            p = j;

        d = i * j;
    }

    if (d == RT(0))
    {
        // no extension
    }
    else
    {
        // with extension
        RT pp = p * p;
        RT ppp = pp * p;
        RT cppp = c * ppp;
        ERT pd = p * ERT(RT(0), d);

        ERT e[4] = { cppp - ( pp + pd ),
                     cppp - ( pp - pd ),
                     cppp - (-pp + pd ) ,
                     cppp - (-pp - pd ) };

        Matrix_ERT m = Matrix_ERT(Spin_quadric_3(Predicate_g_3(k, l, a, b, c)).matrix());

        Matrix_ERT me[4];

    //    for (int z = 0; z < 4; ++z)
    //        me[z] = m - Matrix_ERT::diagonal(e[z]);
    }
}

template<class Kernel_>
const typename Predicate_g_parametrization_3<Kernel_>::Vector_3 &Predicate_g_parametrization_3<Kernel_>::b() const
{
    return m_b;
}
} // namespace CS
