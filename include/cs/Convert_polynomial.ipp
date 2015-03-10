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
#include "Convert_polynomial.h"

template<class Kernel_>
typename Kernel_::Algebraic_kernel_with_sqrt::Polynomial_1 convert_polynomial(
        const typename Kernel_::Hom_polynomial_with_sqrt &hom_polynomial)
{
    CS_BENCHMARK_POINT();

    typedef typename Kernel_::Algebraic_kernel_with_sqrt              Algebraic_kernel_with_sqrt;
    typedef typename Kernel_::RT                                      RT;

    typedef typename Algebraic_kernel_with_sqrt::Polynomial_1         Polynomial_1;

    // if the polynomial is zero, do not do anything
    if (hom_polynomial.degree() == -1)
        return Polynomial_1();

    // read coefficients
    std::vector<CGAL::Sqrt_extension<CGAL::Gmpz, CGAL::Gmpz> > coeffs;
    coeffs.reserve(hom_polynomial.degree() + 1);

    // read root extended or not
    if (hom_polynomial[0].is_extended())
    {
        RT root = hom_polynomial[0].root();
        mpz_t root_mpz;
        root_mpz[0] = root.bigint_rep();

        for (int c = 0; c <= hom_polynomial.degree(); ++c)
        {
            RT a0 = hom_polynomial[c].a0();
            RT a1 = hom_polynomial[c].a1();

            mpz_t a0_mpz, a1_mpz;
            a0_mpz[0] = a0.bigint_rep();
            a1_mpz[0] = a1.bigint_rep();

            coeffs.push_back(
                CGAL::Sqrt_extension<CGAL::Gmpz, CGAL::Gmpz>(
                    CGAL::Gmpz(a0_mpz),
                    CGAL::Gmpz(a1_mpz),
                    CGAL::Gmpz(root_mpz)));
        }
    }
    else
    {
        for (int c = 0; c <= hom_polynomial.degree(); ++c)
        {
            RT a0 = hom_polynomial[c].a0();

            mpz_t a0_mpz;
            a0_mpz[0] = a0.bigint_rep();

            coeffs.push_back(
                    CGAL::Sqrt_extension<CGAL::Gmpz, CGAL::Gmpz>(
                        CGAL::Gmpz(a0_mpz)));
       }
    }

    return Polynomial_1(coeffs.begin(), coeffs.end());
}

template<class Kernel_>
typename Kernel_::Algebraic_kernel_with_sqrt::Polynomial_1 extend_simple_polynomial(
        const typename Kernel_::Algebraic_kernel::Polynomial_1 &polynomial)
{
    CS_BENCHMARK_POINT();

    typedef typename Kernel_::Algebraic_kernel_with_sqrt                   Algebraic_kernel_with_sqrt;
    typedef typename Algebraic_kernel_with_sqrt::Polynomial_1         Polynomial_1;

    if (polynomial.degree() == -1)
        return Polynomial_1();

    // read coefficients
    std::vector<CGAL::Sqrt_extension<CGAL::Gmpz, CGAL::Gmpz> > coeffs;
    coeffs.reserve(polynomial.degree() + 1);

    for (int i = 0; i <= polynomial.degree(); ++i)
    {
        coeffs.push_back(
                CGAL::Sqrt_extension<CGAL::Gmpz, CGAL::Gmpz>(
                        polynomial[i]));
    }

    return Polynomial_1(coeffs.begin(), coeffs.end());
}

template<class Kernel_>
void apart_polynomial(
        const typename Kernel_::Algebraic_kernel_with_sqrt::Polynomial_1 polynomial,
        typename Kernel_::Algebraic_kernel::Polynomial_1 &left_polynomial,
        CGAL::Gmpz &d,
        typename Kernel_::Algebraic_kernel::Polynomial_1 &right_polynomial)
{
    CS_BENCHMARK_POINT();

    typedef typename Kernel_::Algebraic_kernel                             Algebraic_kernel;
    //typedef typename Kernel_::Algebraic_kernel_with_sqrt                   Algebraic_kernel_with_sqrt;

    typedef typename Algebraic_kernel::Polynomial_1                   Simple_polynomial_1;
    //typedef typename Algebraic_kernel_with_sqrt::Polynomial_1         Polynomial_1;

    // a common polynomial degree
    int degree = polynomial.degree();

    if (degree == -1)
    {
        left_polynomial = Simple_polynomial_1();
        right_polynomial = Simple_polynomial_1();
        return;
    }

    // read coefficients
    std::vector<CGAL::Gmpz> left_coeffs, right_coeffs;
    left_coeffs.reserve(degree + 1);
    right_coeffs.reserve(degree + 1);

    // make apart
    for (int i = 0; i <= degree; ++i)
    {
        CGAL::Sqrt_extension<CGAL::Gmpz, CGAL::Gmpz> coeff = polynomial[i];

        left_coeffs.push_back(coeff.a0());
        right_coeffs.push_back(coeff.a1());
        d = coeff.root();
    }

    left_polynomial = Simple_polynomial_1(left_coeffs.begin(), left_coeffs.end());
    right_polynomial = Simple_polynomial_1(right_coeffs.begin(), right_coeffs.end());
}

template<class Kernel_>
bool is_polynomial_extended(
        const typename Kernel_::Algebraic_kernel_with_sqrt::Polynomial_1 polynomial)
{
    CS_BENCHMARK_POINT();

    int degree = polynomial.degree();

    if (degree == -1)
        return false;

    return polynomial[0].root().sign() != CGAL::ZERO;
}

template<class Kernel_>
typename Kernel_::Algebraic_kernel::Polynomial_1 remove_polynomial_extension(
        const typename Kernel_::Algebraic_kernel_with_sqrt::Polynomial_1 polynomial)
{
    CS_BENCHMARK_POINT();

    typedef typename Kernel_::Algebraic_kernel                        Algebraic_kernel;
    typedef typename Algebraic_kernel::Polynomial_1                   Simple_polynomial_1;

    int degree = polynomial.degree();

    if (degree == -1)
        return Simple_polynomial_1();

    assert(polynomial[0].root().sign() == CGAL::ZERO);

    // read coefficients
    std::vector<CGAL::Gmpz> coeffs, right_coeffs;
    coeffs.reserve(degree + 1);

    // make apart
    for (int i = 0; i <= degree; ++i)
    {
        CGAL::Sqrt_extension<CGAL::Gmpz, CGAL::Gmpz> coeff = polynomial[i];

        coeffs.push_back(coeff.a0());
    }

    return Simple_polynomial_1(coeffs.begin(), coeffs.end());
}

template<class Kernel_>
void raw_polynomial_alloc(
        const typename Kernel_::Algebraic_kernel::Polynomial_1 polynomial,
        mpz_t *&coeffs,
        unsigned long &degree)
{
    CS_BENCHMARK_POINT();

    assert(polynomial.degree() >= 0);

    degree = static_cast<unsigned long>(polynomial.degree());
    coeffs = static_cast<mpz_t *>(malloc((degree + 1) * sizeof(mpz_t)));

    for (unsigned long i = 0; i <= degree; ++i)
    {
        mpz_init(coeffs[i]);
        mpz_set(coeffs[i], polynomial[i].mpz());
    }
}

template<class Kernel_>
void raw_polynomial_free(
        mpz_t *&coeffs,
        unsigned long &degree)
{
    CS_BENCHMARK_POINT();

    for (unsigned long i = 0; i <= degree; ++i)
        mpz_clear(coeffs[i]);

    free(coeffs);

    degree = 0;
    coeffs = 0;
}
