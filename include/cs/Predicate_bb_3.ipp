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
#include "Predicate_bb_3.h"

namespace CS
{
template<class Kernel>
Predicate_bb_3<Kernel>::Predicate_bb_3()
    : m_b(Vector_3(RT(0), RT(0), RT(0))),
      m_l(RT(0)),
      m_a(Vector_3(RT(0), RT(0), RT(0))),
      m_r(RT(0))
{
}

template<class Kernel>
Predicate_bb_3<Kernel>::Predicate_bb_3(
        const Vector_3 &b,
        const RT &l,
        const Vector_3 &a,
        const RT &r)
    : m_b(b),
      m_l(l),
      m_a(a),
      m_r(r)
{
    // construct corresponding H predicate
    construct(b, l, a, r);
}

template<class Kernel>
Predicate_bb_3<Kernel>::Predicate_bb_3(
        const Ball_3 &b,
        const Ball_3 &a)
{
    // construct corresponding H predicate
    construct(b.center(), b.radius(), a.center(), a.radius());
}

template<class Kernel>
void Predicate_bb_3<Kernel>::construct(
        const Vector_3 &b,
        const RT &l,
        const Vector_3 &a,
        const RT &r)
{
    m_b = b;
    m_l = l;
    m_a = a;
    m_r = r;

    // init predicate
    RT aa = m_a.squared_length();
    RT bb = m_b.squared_length();
    Vector_3 normal = RT(2) * m_b;

    m_predicate = Predicate_h_3(m_a,
                                Plane_3(normal.x(), normal.y(), normal.z(),
                                        (m_r + m_l) * (m_r + m_l) - aa - bb));
}

template<class Kernel>
const typename Predicate_bb_3<Kernel>::Vector_3 &Predicate_bb_3<Kernel>::b() const
{
    return m_b;
}

template<class Kernel>
const typename Predicate_bb_3<Kernel>::RT &Predicate_bb_3<Kernel>::l() const
{
    return m_l;
}

template<class Kernel>
const typename Predicate_bb_3<Kernel>::Vector_3 &Predicate_bb_3<Kernel>::a() const
{
    return m_a;
}

template<class Kernel>
const typename Predicate_bb_3<Kernel>::RT &Predicate_bb_3<Kernel>::r() const
{
    return m_r;
}

template<class Kernel>
const typename Predicate_bb_3<Kernel>::Predicate_h_3 *Predicate_bb_3<Kernel>::sub_predicates() const
{
    return &m_predicate;
}

template<class Kernel>
bool Predicate_bb_3<Kernel>::evaluate(const bool *signs) const
{
    // evaluate collision predicate
    return signs[0];
}
} // namespace CS
