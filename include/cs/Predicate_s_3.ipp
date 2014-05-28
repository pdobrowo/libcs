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
#include "Predicate_s_3.h"

namespace CS
{
template<class Kernel_>
Predicate_s_3<Kernel_>::Predicate_s_3()
    : m_k(Vector_3(RT(0), RT(0), RT(0))),
      m_l(Vector_3(RT(0), RT(0), RT(0))),
      m_a(Vector_3(RT(0), RT(0), RT(0))),
      m_b(Vector_3(RT(0), RT(0), RT(0)))
{
}

template<class Kernel_>
Predicate_s_3<Kernel_>::Predicate_s_3(const Vector_3  &k,
        const Vector_3  &l,
        const Vector_3  &a,
        const Vector_3  &b)
    : m_k(k),
      m_l(l),
      m_a(a),
      m_b(b)
{
}

template<class Kernel_>
const typename Predicate_s_3<Kernel_>::Vector_3 &Predicate_s_3<Kernel_>::k() const
{
    return m_k;
}

template<class Kernel_>
const typename Predicate_s_3<Kernel_>::Vector_3 &Predicate_s_3<Kernel_>::l() const
{
    return m_l;
}

template<class Kernel_>
const typename Predicate_s_3<Kernel_>::Vector_3 &Predicate_s_3<Kernel_>::a() const
{
    return m_a;
}

template<class Kernel_>
const typename Predicate_s_3<Kernel_>::Vector_3 &Predicate_s_3<Kernel_>::b() const
{
    return m_b;
}

template<class Kernel_>
Predicate_s_3<Kernel_> Predicate_s_3<Kernel_>::opposite() const
{
    // Swap k & l components or a & b components
    return Predicate_s_3<Kernel_>(l(), k(), a(), b());
}

template<class Kernel_>
std::ostream &operator <<(std::ostream &os, const Predicate_s_3<Kernel_> &predicate)
{
    return (os << "[" << predicate.k() << ";" << predicate.l() << ";"
                      << predicate.a() << ";" << predicate.b() << "]");
}
} // namespace CS
