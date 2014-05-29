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

    const RT zero(0);
    const RT one(1);

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

    RT u = CGAL::cross_product(k, l).squared_length() * (a - b).squared_length();
    RT v = (k - l).squared_length() * CGAL::cross_product(a, b).squared_length();

    cout<<"_2u="<<u<<endl;cout<<"_2v="<<v<<endl;

    RT r;
    ERT p, d;

    if (u == zero)
    {
        r = one;
        p = ERT(v, zero, r);
        d = ERT(zero, zero, r);
    }
    else if (v == zero)
    {
        r = one;
        p = ERT(u, zero, r);
        d = ERT(zero, zero, r);
    }
    else
    {
        r = u * v;

        if (u < v)
            p = ERT(u, zero, r);
        else
            p = ERT(v, zero, r);

        d = ERT(zero, one, r);
    }

    cout<<"_2p="<<p<<endl;cout<<"_2d="<<d<<endl;

    ERT pp = p * p;
    ERT ppp = pp * p;
    ERT cppp = c * ppp;
    ERT pd = p * d;

    cout<<"_3pp="<<pp<<endl;cout<<"_3ppp="<<ppp<<endl;cout<<"_3cppp="<<cppp<<endl;cout<<"_3pd="<<pd<<endl;

    ERT eval[4] = { cppp - ( pp + pd ),
                    cppp - ( pp - pd ),
                    cppp - (-pp + pd ) ,
                    cppp - (-pp - pd ) };

    cout<<"_4e0="<<eval[0]<<endl;cout<<"_4e1="<<eval[1]<<endl;cout<<"_4e2="<<eval[2]<<endl;cout<<"_4e3="<<eval[3]<<endl;

    std::sort(eval, eval + 4);
    ERT *uni_eval = std::unique(eval, eval + 4);
    int uni_count = static_cast<int>(uni_eval - eval);

    cout<<"eu="<<uni_count<<endl;
    for(int z=0;z<uni_count;++z)cout<<"se"<<z<<":"<<eval[z]<<endl;

    Matrix_ERT mbase(Spin_quadric_3(Predicate_g_3(k, l, a, b, c)).matrix());
    Matrix_ERT q;

    int total_dim = 0;
    Vector_ERT vk[4];

    for (int uni_idx = 0; uni_idx < uni_count; ++uni_idx)
    {
        Matrix_ERT meval(mbase);

        for (int i = 0; i < 4; ++i)
            meval.set(i, i, meval.get(i, i) - eval[uni_idx]);

        int dim = meval.kernel(vk);

        cout<<"dim="<<dim<<endl;

        for (int i = 0; i < dim; ++i)
            vk[total_dim++] = vk[i];
    }

    //assert(total_dim == 4);
    cout<<"total_dim="<<total_dim<<endl;
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
