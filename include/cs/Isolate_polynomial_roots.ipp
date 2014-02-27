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
#include "Isolate_polynomial_roots.h"

namespace CS
{
template<class Kernel_, typename OutputIterator>
void isolate_polynomial_roots_cgal(
        const typename Kernel_::Algebraic_kernel_with_sqrt::Polynomial_1 &polynomial,
        OutputIterator output)
{
    //CS_BENCHMARK_POINT();

    typedef typename Kernel_::Algebraic_kernel_with_sqrt                 Algebraic_kernel_with_sqrt;
    //typedef typename Algebraic_kernel_with_sqrt::Polynomial_1       Polynomial_1;
    typedef typename Algebraic_kernel_with_sqrt::Algebraic_real_1   Algebraic_real_1;
    typedef typename Algebraic_kernel_with_sqrt::Multiplicity_type  Multiplicity_type;
    typedef typename Algebraic_kernel_with_sqrt::Solve_1            Solve_1;
    typedef typename Kernel_::Algebraic_kernel                 Algebraic_kernel;
    //typedef typename Algebraic_kernel::Polynomial_1       Simple_polynomial_1;
    typedef typename Algebraic_kernel::Algebraic_real_1   Simple_algebraic_real_1;
    typedef typename Algebraic_kernel::Multiplicity_type  Simple_multiplicity_type;
    typedef typename Algebraic_kernel::Solve_1            Simple_solve_1;

    if (is_polynomial_extended<Kernel_>(polynomial))
    {
        CS_BENCHMARK_POINT();

        //
        // Random H3 predicate test: 32% of all isolations are extended
        //
        Algebraic_kernel_with_sqrt ak;
        typedef std::vector<std::pair<Algebraic_real_1, Multiplicity_type> > Roots;
        Roots roots;

        Solve_1 solve_1 = ak.solve_1_object();
        solve_1(polynomial, std::back_inserter(roots));

        // copy roots
        for (typename Roots::const_iterator iterator = roots.begin(); iterator != roots.end(); ++iterator)
            (*output)++ = iterator->first;
    }
    else
    {
        CS_BENCHMARK_POINT();

        //
        // Random H3 predicate test: 68% of all isolations are simple
        //

        // do not go into sqrt computations which are 5 times more expensive
        Algebraic_kernel ak;
        typedef std::vector<std::pair<Simple_algebraic_real_1, Simple_multiplicity_type> > Roots;
        Roots roots;

        Simple_solve_1 solve_1 = ak.solve_1_object();
        solve_1(remove_polynomial_extension<Kernel_>(polynomial), std::back_inserter(roots));

        // copy roots
        for (typename Roots::const_iterator iterator = roots.begin(); iterator != roots.end(); ++iterator)
        {
#if 1
            (*output)++ = Algebraic_real_1(polynomial,
                                           iterator->first.low(),
                                           iterator->first.high()); // convert back to extended real
#else
            (*output)++ = Algebraic_real_1(polynomial,
                                           iterator->first.inf(),
                                           iterator->first.sup()); // convert back to extended real
#endif
        }
    }
}

template<class Kernel_>
typename Kernel_::Algebraic_kernel_with_sqrt::Algebraic_real_1 interval_to_algebraic_real(
        const typename Kernel_::Algebraic_kernel_with_sqrt::Polynomial_1 &polynomial,
        const interval &z)
{
    CS_BENCHMARK_POINT();

    typedef typename Kernel_::Algebraic_kernel_with_sqrt            Algebraic_kernel_with_sqrt;

    typedef typename Algebraic_kernel_with_sqrt::Algebraic_real_1   Algebraic_real_1;

    CGAL::Gmpq low;
    CGAL::Gmpq high;

    if (z.k <= 0)
    {
        low = CGAL::Gmpz(z.c);
    }
    else
    {
        low = CGAL::Gmpq(
            CGAL::Gmpz(z.c),
            CGAL::Gmpz(1) << z.k);
    }

    if (z.isexact == 1)
    {
        // just the same
        high = low;
    }
    else
    {
        if (z.k <= 0)
        {
            high = (CGAL::Gmpz(1) << -z.k) + CGAL::Gmpz(z.c);
        }
        else
        {
            high = CGAL::Gmpq(
                CGAL::Gmpz(z.c) + 1,
                CGAL::Gmpz(1) << z.k);
        }
    }

    return Algebraic_real_1(polynomial, low, high);
}

#if 0

template<class Kernel_, typename OutputIterator>
void isolate_polynomial_roots_hanrot(
        const typename Kernel_::Algebraic_kernel_with_sqrt::Polynomial_1 &polynomial,
        OutputIterator output)
{
    CS_BENCHMARK_POINT();

    typedef typename Kernel_::Hom_polynomial                             Hom_polynomial;
    typedef typename Kernel_::Hom_polynomial_with_sqrt                   Hom_polynomial_with_sqrt;

    typedef typename Kernel_::Algebraic_kernel_with_sqrt                 Algebraic_kernel_with_sqrt;

    typedef typename Algebraic_kernel_with_sqrt::Polynomial_1       Polynomial_1;
    typedef typename Algebraic_kernel_with_sqrt::Algebraic_real_1   Algebraic_real_1;

    //
    // polynomial = left_polynomial + sqrt(d) * right_polynomial
    // big_polynomial = (left_polynomial + sqrt(d) * right_polynomial) * (left_polynomial - sqrt(d) * right_polynomial) =
    //                = (left_polynomial^2 - d * right_polynomial^2)
    //
    // num_zeroes(polynomial) = k
    // num_zeroes(big_polynomial) = l
    //
    // refine zeroes intervals that are not-zeroes of polynomial until there are k left zeroes, which must be
    // zeroes of polynomial
    //
    // FIXME: proof
    typedef typename Kernel_::Algebraic_kernel            Algebraic_kernel;

    typedef typename Algebraic_kernel::Polynomial_1       Simple_polynomial_1;

    Simple_polynomial_1 left_polynomial, right_polynomial;
    CGAL::Gmpz d;

    apart_polynomial<Kernel_>(polynomial, left_polynomial, d, right_polynomial);

    // raw poly
    mpz_t *coeffs;
    unsigned long degree;

    // is right polynomial zero ?
    if (d.sign() == CGAL::ZERO)
    {
        // simple case - right polynomial is zero
        raw_polynomial_alloc<Kernel_>(left_polynomial, coeffs, degree);

        // run Hanrot's Uspensky algorithm
        unsigned int nroots = 0;
        interval *roots;

        roots = Uspensky(coeffs, degree, &nroots);

        for (unsigned int i = 0; i < nroots; ++i)
            (*output)++ = interval_to_algebraic_real<Kernel_>(polynomial, roots[i]);

        // free raw poly
        raw_polynomial_free(coeffs, degree);
    }
    else
    {
        // general case - both polynomials
        Simple_polynomial_1 big_polynomial = left_polynomial * left_polynomial - d * right_polynomial * right_polynomial;
        Polynomial_1 extended_big_polynomial = extend_simple_polynomial<Kernel_>(big_polynomial);

        raw_polynomial_alloc<Kernel_>(big_polynomial, coeffs, degree);

        // run Hanrot's Uspensky algorithm
        unsigned int nroots = 0;
        interval *roots;

        roots = Uspensky(coeffs, degree, &nroots);

        // take all roots
        std::list<Algebraic_real_1> all_roots;

        for (unsigned int i = 0; i < nroots; ++i)
            all_roots.push_back(
                interval_to_algebraic_real<Kernel_>(extended_big_polynomial, roots[i]));

        // free raw poly
        raw_polynomial_free(coeffs, degree);

        // definition of isolated square root of d
        //Algebraic_real_1 sqrt(Simple_polynomial_1)

        // now, refine until all foreign roots are discarded
        for (typename std::list<Algebraic_real_1>::const_iterator rootIterator = all_roots.begin(); rootIterator != all_roots.end(); ++rootIterator)
        {
            const Algebraic_real_1 &root = *rootIterator;

            CGAL::Gmpfi left_value = polynomial_interval_value(left_polynomial, root);
            CGAL::Gmpfi right_value = polynomial_interval_value(right_polynomial, root);

            CGAL::Gmpfi sqrt_d = CGAL::sqrt(CGAL::Gmpfi(2));

            CGAL::Gmpfi value_this = left_value + sqrt_d * right_value;
            CGAL::Gmpfi value_other = left_value - sqrt_d * right_value;

            if (CGAL::certainly_not(value_this.is_zero()))
            {
                // the zero is certainly from the other polynomial
//                std::cerr << "certainly other poly zero" << std::endl;
            }
            else if (CGAL::certainly_not(value_other.is_zero()))
            {
                // the zero is certainly from the this polynomial
//                std::cerr << "certainly this poly zero" << std::endl;
            }
            else
            {
                // unknown, refine intervals
//                std::cerr << "dunno poly zero" << std::endl;
            }
        }
    }
}

#endif

template<class Kernel_, typename OutputIterator>
void isolate_polynomial_roots_contfrac(
        const typename Kernel_::Algebraic_kernel_with_sqrt::Polynomial_1 &polynomial,
        OutputIterator output)
{
    //CS_BENCHMARK_POINT();

    typedef typename Kernel_::Algebraic_kernel_with_sqrt            Algebraic_kernel_with_sqrt;
    //typedef typename Algebraic_kernel_with_sqrt::Polynomial_1       Polynomial_1;
    typedef typename Algebraic_kernel_with_sqrt::Algebraic_real_1   Algebraic_real_1;
    typedef typename Algebraic_kernel_with_sqrt::Multiplicity_type  Multiplicity_type;
    typedef typename Algebraic_kernel_with_sqrt::Solve_1            Solve_1;
    //typedef typename Kernel_::Algebraic_kernel                      Algebraic_kernel;
    //typedef typename Algebraic_kernel::Polynomial_1                 Simple_polynomial_1;
    //typedef typename Algebraic_kernel::Algebraic_real_1             Simple_algebraic_real_1;
    //typedef typename Algebraic_kernel::Multiplicity_type            Simple_multiplicity_type;
    //typedef typename Algebraic_kernel::Solve_1                      Simple_solve_1;

    //
    // WARNING: sometimes the contfrac implementation hangs up!
    //          this is clearly a bug in realroot library
    //
    //          update the library and/or report a bug
    //
#if 1
    if (true)
#else
    if (is_polynomial_extended<Kernel_>(polynomial))
#endif
    {
        CS_BENCHMARK_POINT();

        //
        // Random H3 predicate test: 32% of all isolations are extended
        //
        Algebraic_kernel_with_sqrt ak;
        typedef std::vector<std::pair<Algebraic_real_1, Multiplicity_type> > Roots;
        Roots roots;

        Solve_1 solve_1 = ak.solve_1_object();
        solve_1(polynomial, std::back_inserter(roots));

        // copy roots
        for (typename Roots::const_iterator iterator = roots.begin(); iterator != roots.end(); ++iterator)
            (*output)++ = iterator->first;
    }
    else
    {
        //
        // Random H3 predicate test: 68% of all isolations are simple
        //
        typedef std::vector<contfrac_interval> Roots;
        Roots roots;

        // Isolation
        {
            CS_BENCHMARK_POINT();

            int degree = polynomial.degree();

            // read coefficients
            mpz_t *coeffs = new mpz_t[degree + 1];

            for (int i = 0; i <= degree; ++i)
                mpz_init_set(coeffs[i], polynomial[i].a0().mpz());

            contfrac(coeffs, degree, roots);

            for (int i = 0; i <= degree; ++i)
                mpz_clear(coeffs[i]);

            delete [] coeffs;
        }

        // Root copying
        {
            CS_BENCHMARK_POINT();

            // FIXME: This shit is more than 14 times slower than confrac root isolation!!!!!

            for (typename Roots::const_iterator iterator = roots.begin(); iterator != roots.end(); ++iterator)
                (*output)++ = Algebraic_real_1(polynomial, iterator->low, iterator->high);
        }
    }
}
} // namespace CS
