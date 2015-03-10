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
#include "Predicate_tt_3.h"

namespace CS
{
template<class Kernel_>
Predicate_tt_3<Kernel_>::Predicate_tt_3()
    : m_k(Vector_3(RT(0), RT(0), RT(0))),
      m_l(Vector_3(RT(0), RT(0), RT(0))),
      m_m(Vector_3(RT(0), RT(0), RT(0))),
      m_a(Vector_3(RT(0), RT(0), RT(0))),
      m_b(Vector_3(RT(0), RT(0), RT(0))),
      m_c(Vector_3(RT(0), RT(0), RT(0)))
{
}

template<class Kernel_>
Predicate_tt_3<Kernel_>::Predicate_tt_3(
        const Vector_3  &k,
        const Vector_3  &l,
        const Vector_3  &m,
        const Vector_3  &a,
        const Vector_3  &b,
        const Vector_3  &c)
    : m_k(k),
      m_l(l),
      m_m(m),
      m_a(a),
      m_b(b),
      m_c(c)
{
    // construct the matrix of s3 predicates
    m_predicates[0] = Predicate_s_3(m_k, m_l, m_a, m_b);
    m_predicates[1] = Predicate_s_3(m_k, m_l, m_b, m_c);
    m_predicates[2] = Predicate_s_3(m_k, m_l, m_c, m_a);
    m_predicates[3] = Predicate_s_3(m_l, m_m, m_a, m_b);
    m_predicates[4] = Predicate_s_3(m_l, m_m, m_b, m_c);
    m_predicates[5] = Predicate_s_3(m_l, m_m, m_c, m_a);
    m_predicates[6] = Predicate_s_3(m_m, m_k, m_a, m_b);
    m_predicates[7] = Predicate_s_3(m_m, m_k, m_b, m_c);
    m_predicates[8] = Predicate_s_3(m_m, m_k, m_c, m_a);
}

template<class Kernel_>
const typename Predicate_tt_3<Kernel_>::Vector_3 &Predicate_tt_3<Kernel_>::k() const
{
    return m_k;
}

template<class Kernel_>
const typename Predicate_tt_3<Kernel_>::Vector_3 &Predicate_tt_3<Kernel_>::l() const
{
    return m_l;
}

template<class Kernel_>
const typename Predicate_tt_3<Kernel_>::Vector_3 &Predicate_tt_3<Kernel_>::m() const
{
    return m_m;
}

template<class Kernel_>
const typename Predicate_tt_3<Kernel_>::Vector_3 &Predicate_tt_3<Kernel_>::a() const
{
    return m_a;
}

template<class Kernel_>
const typename Predicate_tt_3<Kernel_>::Vector_3 &Predicate_tt_3<Kernel_>::b() const
{
    return m_b;
}

template<class Kernel_>
const typename Predicate_tt_3<Kernel_>::Vector_3 &Predicate_tt_3<Kernel_>::c() const
{
    return m_c;
}

template<class Kernel_>
const typename Predicate_tt_3<Kernel_>::Predicate_s_3 *Predicate_tt_3<Kernel_>::sub_predicates() const
{
    return m_predicates;
}

template<class Kernel_>
bool Predicate_tt_3<Kernel_>::evaluate(const bool *signs) const
{
#if 1
    // note: sign is boolean

    // evaluate collision predicate
    // rule: any row or column has the same signs for each cell
    return (signs[0] == signs[1] && signs[1] == signs[2]) ||
           (signs[3] == signs[4] && signs[4] == signs[5]) ||
           (signs[6] == signs[7] && signs[7] == signs[8]) ||
           (signs[0] == signs[3] && signs[3] == signs[6]) ||
           (signs[1] == signs[4] && signs[4] == signs[7]) ||
           (signs[2] == signs[5] && signs[5] == signs[8]);

#else
    // previous implementation: for integer type

    return (signs[0] >= 0 && signs[1] >= 0 && signs[2] >= 0) ||
           (signs[0] <= 0 && signs[1] <= 0 && signs[2] <= 0) ||
           (signs[3] >= 0 && signs[4] >= 0 && signs[5] >= 0) ||
           (signs[3] <= 0 && signs[4] <= 0 && signs[5] <= 0) ||
           (signs[6] >= 0 && signs[7] >= 0 && signs[8] >= 0) ||
           (signs[6] <= 0 && signs[7] <= 0 && signs[8] <= 0) ||
           (signs[0] >= 0 && signs[3] >= 0 && signs[6] >= 0) ||
           (signs[0] <= 0 && signs[3] <= 0 && signs[6] <= 0) ||
           (signs[1] >= 0 && signs[4] >= 0 && signs[7] >= 0) ||
           (signs[1] <= 0 && signs[4] <= 0 && signs[7] <= 0) ||
           (signs[2] >= 0 && signs[5] >= 0 && signs[8] >= 0) ||
           (signs[2] <= 0 && signs[5] <= 0 && signs[8] <= 0);

#endif
}
} // namespace CS
