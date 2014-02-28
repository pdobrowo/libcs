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
#ifndef LIBCS_BIGINT_H
#define LIBCS_BIGINT_H

#include <LiDIA/bigint.h>

namespace CGAL
{
// LiDIA bigint
typedef LiDIA::bigint Bigint;

// Algebraic structure traits
template <> class Algebraic_structure_traits< Bigint >
    : public Algebraic_structure_traits_base< Bigint, Euclidean_ring_tag > {
public:
    typedef Tag_true            Is_exact;
    typedef Tag_false           Is_numerical_sensitive;

    typedef INTERN_AST::Is_square_per_sqrt< Type > Is_square;

    class Integral_division
        : public std::binary_function< Type, Type, Type > {
    public:
        Type operator()( const Type& x, const Type& y ) const {
            Bigint result;
            mpz_t mpz_result, mpz_x, mpz_y;
            mpz_result[0] = result.bigint_rep();
            mpz_x[0] = x.bigint_rep();
            mpz_y[0] = y.bigint_rep();
            mpz_divexact(mpz_result, mpz_x, mpz_y);
            CGAL_postcondition_msg(result * y == x, "exact_division failed\n");
            return result;
        }
    };

    class Gcd
        : public std::binary_function< Type, Type, Type > {
    public:
        Type operator()( const Type& x, const Type& y ) const {
            Bigint result;
            mpz_t mpz_result, mpz_x, mpz_y;
            mpz_result[0] = result.bigint_rep();
            mpz_x[0] = x.bigint_rep();
            mpz_y[0] = y.bigint_rep();
            mpz_gcd(mpz_result, mpz_x, mpz_y);
            return result;
        }

        Type operator()( const Type& x, const int& y ) const {
          if (y > 0)
          {
              Bigint result;
              mpz_t mpz_result, mpz_x;
              mpz_result[0] = result.bigint_rep();
              mpz_x[0] = x.bigint_rep();
              mpz_gcd_ui(mpz_result, mpz_x, y);
              return result;
          }
          return CGAL_NTS gcd(x, Bigint(y));
        }

        Type operator()( const int& x, const Type& y ) const {
          return CGAL_NTS gcd(Bigint(x), y );
        }
    };

    typedef INTERN_AST::Div_per_operator< Type > Div;
    typedef INTERN_AST::Mod_per_operator< Type > Mod;

    class Sqrt
        : public std::unary_function< Type, Type > {
    public:
        Type operator()( const Type& x ) const {
            Bigint result;
            mpz_t mpz_result, mpz_x;
            mpz_result[0] = result.bigint_rep();
            mpz_x[0] = x.bigint_rep();
            mpz_sqrt(mpz_result, mpz_x);
            return result;
        }
    };
};

template <> class Real_embeddable_traits< Bigint >
    : public INTERN_RET::Real_embeddable_traits_base< Bigint, CGAL::Tag_true > {
public:
    class Sgn
        : public std::unary_function< Type, ::CGAL::Sign > {
    public:
        ::CGAL::Sign operator()( const Type& x ) const {
            if (x.is_gt_zero())
                return CGAL::POSITIVE;
            else if (x.is_lt_zero())
                return CGAL::NEGATIVE;
            else
                return CGAL::ZERO;
        }
    };

    class To_double
        : public std::unary_function< Type, double > {
    public:
        double operator()( const Type& x ) const {
            return x.dbl();
        }
    };

    class To_interval
        : public std::unary_function< Type, std::pair< double, double > > {
    public:
        std::pair<double, double> operator()( const Type& x ) const {
            mpz_t mpz;
            mpz[0] = x.bigint_rep();
            mpfr_t y;
            mpfr_init2 (y, 53); /* Assume IEEE-754 */
            mpfr_set_z (y, mpz, GMP_RNDD);
            double i = mpfr_get_d (y, GMP_RNDD); /* EXACT but can overflow */
            mpfr_set_z (y, mpz, GMP_RNDU);
            double s = mpfr_get_d (y, GMP_RNDU); /* EXACT but can overflow */
            mpfr_clear (y);
            return std::pair<double, double>(i, s);
        }
    };
};

template<> class Algebraic_structure_traits< Quotient<Bigint> >
    : public INTERN_QUOTIENT::Algebraic_structure_traits_quotient_base<Quotient<Bigint> >{
    // specialization of to double functor
public:
    typedef Quotient<Bigint> Type;

    struct To_double: public std::unary_function<Quotient<Bigint>, double>{
        double operator()(const Quotient<Bigint>& quot){
            mpq_t  mpQ;
            mpq_init(mpQ);
            const Bigint& n = quot.numerator();
            const Bigint& d = quot.denominator();
            mpz_t mpz_n, mpz_d;
            mpz_n[0] = n.bigint_rep();
            mpz_d[0] = d.bigint_rep();
            mpz_set(mpq_numref(mpQ), mpz_n);
            mpz_set(mpq_denref(mpQ), mpz_d);

            mpq_canonicalize(mpQ);

            double ret = mpq_get_d(mpQ);
            mpq_clear(mpQ);
            return ret;
        }
    };
};

//
// Needs_parens_as_product
//
template <>
struct Needs_parens_as_product<Bigint> {
  bool operator()(const Bigint& x) {
    return CGAL_NTS is_negative(x);
  }
};

/*! \ingroup NiX_Modular_traits_spec
 *  \brief a model of concept ModularTraits,
 *  specialization of NiX::Modular_traits.
 */
template<>
class Modular_traits< Bigint > {
  typedef Residue RES;
 public:
    typedef Bigint NT;
    typedef CGAL::Tag_true Is_modularizable;
    typedef Residue Residue_type;

    struct Modular_image{
        Residue_type operator()(const NT& a){
          NT tmp_1(a % NT(RES::get_current_prime()));
          mpz_t mpz_tmp_1;
          mpz_tmp_1[0] = tmp_1.bigint_rep();
          return CGAL::Residue(int(mpz_get_si(mpz_tmp_1)));
        }
    };
    struct Modular_image_representative{
        NT operator()(const Residue_type& x){
          return NT(x.get_value());
        }
    };
};
} // namespace CGAL

namespace CS
{
inline to_double(const CGAL::Bigint &x)
{
    return LiDIA::dbl(x);
}
} // namespace CS

#include "Bigint.ipp"

#endif // LIBCS_BIGINT_H
