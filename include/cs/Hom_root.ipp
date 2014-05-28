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
#include "Hom_root.h"

namespace CS
{
template<class Kernel_>
Hom_root<Kernel_>::Hom_root(bool infinity, const Algebraic_real_1 &parameter)
    : m_infinity(infinity),
      m_root(parameter)
{
}

template<class Kernel_>
bool Hom_root<Kernel_>::is_infinite() const
{
    return m_infinity;
}

template<class Kernel_>
bool Hom_root<Kernel_>::is_finite() const
{
    return !is_infinite();
}

template<class Kernel_>
const typename Hom_root<Kernel_>::Algebraic_real_1 &Hom_root<Kernel_>::root() const
{
    return m_root;
}

template<class Kernel_, typename OutputIterator_>
void extract_isolated_hom_roots(
        const typename Kernel_::Hom_polynomial_with_sqrt &hom_polynomial_,
        OutputIterator_ output)
{
    CS_BENCHMARK_POINT();

    //typedef typename Kernel_::Hom_polynomial                             Hom_polynomial;
    typedef typename Kernel_::Hom_polynomial_with_sqrt                   Hom_polynomial_with_sqrt;

    typedef typename Kernel_::Algebraic_kernel_with_sqrt                 Algebraic_kernel_with_sqrt;

    //typedef typename Algebraic_kernel_with_sqrt::Polynomial_1       Polynomial_1;
    typedef typename Algebraic_kernel_with_sqrt::Algebraic_real_1   Algebraic_real_1;
    //typedef typename Algebraic_kernel_with_sqrt::Bound              Bound;
    //typedef typename Algebraic_kernel_with_sqrt::Multiplicity_type  Multiplicity_type;
    //typedef typename Algebraic_kernel_with_sqrt::Solve_1            Solve_1;

    // minimal polynomial
    Hom_polynomial_with_sqrt hom_polynomial = hom_polynomial_;

    // check if the roots are isolated
    if (hom_polynomial.is_zero())
        return;

    // check whether infinity is a solution
    // coefficient x^degree is zero <=> infinity is a solution
    bool infinity_is_a_root = hom_polynomial[hom_polynomial.degree()].is_zero();

    if (infinity_is_a_root)
        hom_polynomial.set_degree(hom_polynomial.degree() - 1);

    // do root isolation
    std::vector<Algebraic_real_1> roots;
    isolate_polynomial_roots<Kernel_>(convert_polynomial<Kernel_>(hom_polynomial), std::back_inserter(roots));

    if (infinity_is_a_root)
    {
        //std::cerr << "  infinity root" << std::endl;

        (*output)++ = Hom_root<Kernel_>(true, Algebraic_real_1());
    }

    //std::cerr << "  number of roots is: " << roots.size() << std::endl;

    for (size_t i = 0; i < roots.size(); ++i)
    {
        //std::cerr << "    root: " << CGAL::to_double(roots[i].first) << " (multi: " << roots[i].second << ")" << std::endl;

        // add new points: a conjugate pair
        (*output)++ = Hom_root<Kernel_>(false, roots[i]);
    }
}

template<class Kernel_>
CGAL::Sign hom_polynomial_with_sqrt_sign_at(
        const typename Kernel_::Hom_polynomial_with_sqrt &hom_polynomial,
        const typename Kernel_::Hom_root &hom_root)
{
    CS_BENCHMARK_POINT();

    typedef typename Kernel_::Algebraic_kernel_with_sqrt              Algebraic_kernel_with_sqrt;
    typedef typename Algebraic_kernel_with_sqrt::Polynomial_1         Polynomial_1;

    // zero polynomial has always sign equal to zero at any point
    if (hom_polynomial.degree() == -1)
        return CGAL::ZERO;

    // if the hom root is infinte or rational it is a simple case
    if (hom_root.is_infinite())
    {
        // the sign at infinity can be deduced from coefficients
        return hom_polynomial[hom_polynomial.degree()].sign();
    }
    else
    {
        Polynomial_1 polynomial = convert_polynomial<Kernel_>(hom_polynomial);

        // prepare polynomials for fast sign computation
        std::vector<Polynomial_1> polynomials;
        polynomial_sign_at_prepare(polynomial, polynomials);

        // maximum number of refinements in fast method
        const size_t MAX_REFINE_COUNT       = 24;
        size_t refine_count = 0;

        // iterate fast method
        do
        {
            CS_BENCHMARK_POINT();

            CGAL::Uncertain<CGAL::Sign> sign = polynomial_sign_at(polynomials, hom_root.root());

            // if the evaluation yield a certain value, return it immediately
            if (CGAL::is_certain(sign))
                return CGAL::get_certain(sign);

            // refine root more
            hom_root.root().refine();
        }
        while (++refine_count <= MAX_REFINE_COUNT);

        // if still uncertain the value may be zero
        {
            CS_BENCHMARK_POINT();

            // bad, we have to use the algorithm
            Algebraic_kernel_with_sqrt algebraic_kernel;
            typename Algebraic_kernel_with_sqrt::Sign_at_1 sign_at_1 = algebraic_kernel.sign_at_1_object();
            return sign_at_1(polynomial, hom_root.root());
        }
    }
}
} // namespace CS
