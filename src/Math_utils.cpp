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
#include <cs/Math_utils.h>

namespace CS
{
namespace Math
{
mpfr_rnd_t gmp_rounding_mode(std::float_round_style r)
{
    switch (r)
    {
        case std::round_toward_infinity: return MPFR_RNDU;
        case std::round_toward_neg_infinity: return MPFR_RNDD;
        case std::round_toward_zero: return MPFR_RNDZ;
        case std::round_to_nearest: return MPFR_RNDN;
        default:
            assert(0);
            return MPFR_RNDN;
    }
}

CGAL::Gmpfr::Precision_type gmp_result_precision(const CGAL::Gmpfr &x)
{
    return (x.get_precision() > CGAL::Gmpfr::get_default_precision() ?
            x.get_precision() :
            CGAL::Gmpfr::get_default_precision());
}

CGAL::Gmpfr sin(const CGAL::Gmpfr &x)
{
    CGAL::Gmpfr result(0, gmp_result_precision(x));
    mpfr_sin(result.fr(), x.fr(), gmp_rounding_mode(CGAL::Gmpfr::get_default_rndmode()));
    return result;
}

CGAL::Gmpfr cos(const CGAL::Gmpfr &x)
{
    CGAL::Gmpfr result(0, gmp_result_precision(x));
    mpfr_cos(result.fr(), x.fr(), gmp_rounding_mode(CGAL::Gmpfr::get_default_rndmode()));
    return result;
}

CGAL::Gmpfr acos(const CGAL::Gmpfr &x)
{
    CGAL::Gmpfr result(0, gmp_result_precision(x));
    mpfr_acos(result.fr(), x.fr(), gmp_rounding_mode(CGAL::Gmpfr::get_default_rndmode()));
    return result;
}

CGAL::Gmpfr atan(const CGAL::Gmpfr &x)
{
    CGAL::Gmpfr result(0, gmp_result_precision(x));
    mpfr_atan(result.fr(), x.fr(), gmp_rounding_mode(CGAL::Gmpfr::get_default_rndmode()));
    return result;
}

CGAL::Gmpfr pow(const CGAL::Gmpfr &x, const CGAL::Gmpfr &y)
{
    CGAL::Gmpfr result(0, std::max(gmp_result_precision(x), gmp_result_precision(y)));
    mpfr_pow(result.fr(), x.fr(), y.fr(), gmp_rounding_mode(CGAL::Gmpfr::get_default_rndmode()));
    return result;
}

CGAL::Gmpfr atan2(const CGAL::Gmpfr &y, const CGAL::Gmpfr &x)
{
    CGAL::Gmpfr result(0, std::max(gmp_result_precision(y), gmp_result_precision(x)));
    mpfr_atan2(result.fr(), y.fr(), x.fr(), gmp_rounding_mode(CGAL::Gmpfr::get_default_rndmode()));
    return result;
}

CGAL::Gmpfr sqrt(const CGAL::Gmpfr &x)
{
    return x.sqrt();
}

CGAL::Gmpfr fabs(const CGAL::Gmpfr &x)
{
    return x.abs();
}

template<>
float pi<float>()
{
    return static_cast<float>(M_PI);
}

template<>
double pi<double>()
{
    return M_PI;
}

template<>
CGAL::Gmpfr pi<CGAL::Gmpfr>()
{
    CGAL::Gmpfr result(0, CGAL::Gmpfr::get_default_precision());
    mpfr_const_pi(result.fr(), gmp_rounding_mode(CGAL::Gmpfr::get_default_rndmode()));
    return result;
}

void print(const CGAL::Gmpfr &x)
{
    mpfr_out_str(stderr, 10, 0, x.fr(), GMP_RNDD);
}
} // namespace Math
} // namespace CS
