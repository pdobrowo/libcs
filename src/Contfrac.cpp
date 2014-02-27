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
#include <cs/Contfrac.h>
#include <cs/Benchmark.h>
#include <realroot/GMP.hpp>
#include <realroot/polynomial_tensor.hpp>
#include <realroot/solver_uv_continued_fraction.hpp>

contfrac_interval::contfrac_interval(const CGAL::Gmpq &low_, const CGAL::Gmpq &high_)
    : low(low_),
      high(high_)
{
}

void contfrac(mpz_t *coeffs, unsigned int degree, std::vector<contfrac_interval> &out)
{
    typedef mmx::GMP::integer   Z;
    typedef mmx::GMP::rational  Q;
    typedef mmx::polynomial<Z, mmx::with<mmx::MonomialTensor> > P;

    P p;

    // Conversion
    {
        CS_BENCHMARK_POINT();

        P m("x"), w("1");
        Z c;

        for (unsigned int i= 0; i <= degree; i++)
        {
            mpz_set(&c.rep(), coeffs[i]);
            p += w * c;
            w = w * m;
        }
    }

    // Isolation
    {
        CS_BENCHMARK_POINT();

        mmx::Seq<mmx::Interval<Q> > roots;
        mmx::solver<Q, mmx::ContFrac<mmx::Isolate> >::solve(roots, p);

        int num_roots = roots.size();
        //std::cerr << "roots:" << std::endl;

        for (int i=0; i<num_roots; ++i)
        {
            //std::cerr << "root: " << roots[i] << std::endl;
            out.push_back(contfrac_interval(CGAL::Gmpq(&roots[i].lower().rep()),
                                            CGAL::Gmpq(&roots[i].upper().rep())));
        }
    }
}
