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

#include <iostream>

namespace CS
{
template<class Kernel_>
Predicate_g_parametrization_3<Kernel_>::Predicate_g_parametrization_3(const Predicate_g_3 &g3)
{
    using std::cout;
    using std::endl;

    const RT zero = RT(0);
    const RT one = RT(1);

    Vector_3 k = g3.k();
    Vector_3 l = g3.l();
    Vector_3 a = g3.a();
    Vector_3 b = g3.b();
    RT c = g3.c();

    cout<<"_0k=["<<k<<"]"<<endl;cout<<"_0l=["<<l<<"]"<<endl;cout<<"_0a=["<<a<<"]"<<endl;cout<<"_0b=["<<b<<"]"<<endl;cout<<"_0c="<<c<<endl;

    if (c < zero)
    {
        std::swap(k, l);
        c = -c;
    }

    cout<<"_1k=["<<k<<"]"<<endl;cout<<"_1l=["<<l<<"]"<<endl;cout<<"_1a=["<<a<<"]"<<endl;cout<<"_1b=["<<b<<"]"<<endl;cout<<"_1c="<<c<<endl;

    RT i = CGAL::cross_product(k, l).squared_length();
    RT j = CGAL::cross_product(a, b).squared_length();

    cout<<"_2i="<<i<<endl;cout<<"_2j="<<j<<endl;

    RT r;
    ERT p, d;

    if (i == zero)
    {
        r = one;
        p = ERT(j, zero, r);
        d = ERT(zero, zero, r);
    }
    else if (j == zero)
    {
        r = one;
        p = ERT(i, zero, r);
        d = ERT(zero, zero, r);
    }
    else
    {
        r = i * j;

        if (i < j)
            p = ERT(i, zero, r);
        else
            p = ERT(j, zero, r);

        d = ERT(zero, one, r);
    }

    cout<<"_2p="<<p<<endl;cout<<"_2d="<<d<<endl;

    ERT pp = p * p;
    ERT ppp = pp * p;
    ERT cppp = c * ppp;
    ERT pd = p * d;

    cout<<"_3pp="<<pp<<endl;cout<<"_3ppp="<<ppp<<endl;cout<<"_3cppp="<<cppp<<endl;cout<<"_3pd="<<pd<<endl;

    ERT e[4] = { cppp - ( pp + pd ),
                 cppp - ( pp - pd ),
                 cppp - (-pp + pd ) ,
                 cppp - (-pp - pd ) };

    cout<<"_4e0="<<e[0]<<endl;cout<<"_4e1="<<e[1]<<endl;cout<<"_4e2="<<e[2]<<endl;cout<<"_4e3="<<e[3]<<endl;

    Matrix_ERT m(Spin_quadric_3(Predicate_g_3(k, l, a, b, c)).matrix());
    Matrix_ERT me[4];

//    for (int z = 0; z < 4; ++z)
//        me[z] = m - Matrix_ERT::diagonal(e[z]);

}

template<class Kernel_>
const typename Predicate_g_parametrization_3<Kernel_>::Vector_3 &Predicate_g_parametrization_3<Kernel_>::b() const
{
    return m_b;
}

template<class Kernel_>
Spin_3<double> Predicate_g_parametrization_3<Kernel_>::approx_evaluate(double u, double v) const
{
    return Spin_3<double>(1, 0, 0, 0);
}
} // namespace CS
