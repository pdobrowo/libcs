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
#include "Polynomial_sign_at.h"

namespace CS
{
template<class Polynomial_1_, class RT_>
CGAL::Uncertain<CGAL::Sign> polynomial_interval_sign_at(
        const std::vector<Polynomial_1_> &polynomials,
        const size_t &level,
        const RT_ &low, const RT_ &high)
{
    // current polynomial
    const Polynomial_1_ &polynomial = polynomials[level];

    // if the polynomial degree is zero, there is a certain answer
    if (polynomial.degree() == 0)
        return CGAL::make_certain(CGAL::sign(polynomial[0]));

    // find out interval sign by using polynomial differential
    CGAL::Uncertain<CGAL::Sign> uncertain_differential_sign = polynomial_interval_sign_at(polynomials, level + 1, low, high);

    if (!uncertain_differential_sign.is_certain())
    {
        // if the dp sign is not certain, it is impossible to determine polynomial interval sign: return uncertain value
        return CGAL::make_uncertain(CGAL::ZERO);
    }

    CGAL::Sign polynomial_sign_low = CGAL::sign(polynomial.evaluate(low));
    CGAL::Sign polynomial_sign_high = CGAL::sign(polynomial.evaluate(high));

    CGAL::Sign differential_sign = CGAL::get_certain(uncertain_differential_sign);

    bool differential_sign_not_negative = (differential_sign != CGAL::NEGATIVE);
    bool differential_sign_not_positive = (differential_sign != CGAL::POSITIVE);

    // the result is certainly positive
    if ((polynomial_sign_low == CGAL::POSITIVE && differential_sign_not_negative) ||
        (polynomial_sign_high == CGAL::POSITIVE && differential_sign_not_positive))
    {
        return CGAL::make_certain(CGAL::POSITIVE);
    }

    // the result is certainly negative
    if ((polynomial_sign_low == CGAL::NEGATIVE && differential_sign_not_positive) ||
        (polynomial_sign_high == CGAL::NEGATIVE && differential_sign_not_negative))
    {
        return CGAL::make_certain(CGAL::NEGATIVE);
    }

    // the result is certainly zero
    if (polynomial_sign_low == CGAL::ZERO && differential_sign == CGAL::ZERO) // this can be also constructed with polynomial_sign_high
    {
        return CGAL::make_certain(CGAL::ZERO);
    }

    // the result is possibly zero
    return CGAL::make_uncertain(CGAL::ZERO);
}

template<class Polynomial_1_, class Algebraic_real_1_>
CGAL::Uncertain<CGAL::Sign> polynomial_sign_at(
        const std::vector<Polynomial_1_> &polynomials,
        const Algebraic_real_1_ &r)
{
    CS_BENCHMARK_POINT();

    return polynomial_interval_sign_at(polynomials, 0, r.low(), r.high());
}

template<class Polynomial_1_>
void polynomial_sign_at_prepare(
        const Polynomial_1_ &polynomial,
        std::vector<Polynomial_1_> &polynomials)
{
    CS_BENCHMARK_POINT();

    // prepare buffer
    polynomials.clear();
    polynomials.reserve(polynomial.degree() + 1); // for degrees: N, N-1, ..., 0

    // always push initial polynomial
    polynomials.push_back(polynomial);

    Polynomial_1_ current_polynomial = polynomial;

    do
    {
        current_polynomial = CGAL::differentiate(current_polynomial);
        polynomials.push_back(current_polynomial);
    }
    while (current_polynomial.degree() > 0);
}
} // namespace CS
