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
#ifndef LIBCS_PREDICATE_TT_3_H
#define LIBCS_PREDICATE_TT_3_H

#include <ostream>

namespace CS
{
// TT3: triangle-triangle collision predicate
//
//      ABC - oriented rotating triangle
//      KLM - oriented stable triangle
//
// predicate evaluates to true if and only if
// there is a collision of two triangles
//
template<class Kernel_>
class Predicate_tt_3
{
    typedef typename Kernel_::RT            RT;
    typedef typename Kernel_::FT            FT;
    typedef typename Kernel_::Vector_3      Vector_3;
    typedef typename Kernel_::Predicate_s_3 Predicate_s_3;

public:
    typedef Kernel_                         Kernel;

    Predicate_tt_3();

    Predicate_tt_3(
            const Vector_3 &k,
            const Vector_3 &l,
            const Vector_3 &m,
            const Vector_3 &a,
            const Vector_3 &b,
            const Vector_3 &c);

    const Vector_3 &k() const;
    const Vector_3 &l() const;
    const Vector_3 &m() const;
    const Vector_3 &a() const;
    const Vector_3 &b() const;
    const Vector_3 &c() const;

    // access to sub-predicates
    //
    // index:
    // i = { 'kl', 'lm', 'mk' }
    // j = { 'ab', 'bc', 'ca' }
    //
    // [0:kl/ab] [1:kl/bc] [2:kl/ca]
    // [3:lm/ab] [4:lm/bc] [5:lm/ca]
    // [6:mk/ab] [7:mk/bc] [8:mk/ca]
    //
    // return value is a pointer to a table with nine elements
    static const size_t SUB_PREDICATE_COUNT = 9;
    typedef Predicate_s_3 Sub_predicate;

    const Predicate_s_3 *sub_predicates() const;
    bool evaluate(const bool *signs) const;

private:
    Vector_3    m_k;
    Vector_3    m_l;
    Vector_3    m_m;
    Vector_3    m_a;
    Vector_3    m_b;
    Vector_3    m_c;

    Predicate_s_3 m_predicates[9];
};
} // namespace CS

#include "Predicate_tt_3.ipp"

#endif // LIBCS_PREDICATE_TT_3_H
