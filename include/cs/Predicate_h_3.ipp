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
#include "Predicate_h_3.h"

namespace CS
{
template<class Kernel_>
Predicate_h_3<Kernel_>::Predicate_h_3()
    : m_b(Vector_3(RT(0), RT(0), RT(0))),
      m_p(Plane_3(RT(0), RT(0), RT(0), RT(0)))
{
}

template<class Kernel_>
Predicate_h_3<Kernel_>::Predicate_h_3(const Vector_3 &b,
                                const Plane_3 &p)
    : m_b(b),
      m_p(p)
{
}

template<class Kernel_>
const typename Predicate_h_3<Kernel_>::Vector_3 &Predicate_h_3<Kernel_>::b() const
{
    return m_b;
}

template<class Kernel_>
const typename Predicate_h_3<Kernel_>::Plane_3 &Predicate_h_3<Kernel_>::p() const
{
    return m_p;
}

template<class Kernel_>
Predicate_h_3<Kernel_> Predicate_h_3<Kernel_>::opposite() const
{
    return Predicate_h_3<Kernel_>(b(), p().opposite());
}

template<class Kernel_>
std::ostream &operator <<(std::ostream &os, const Predicate_h_3<Kernel_> &predicate)
{
    return (os << "[" << predicate.b() << ";" << predicate.p() << "]");
}

template<class Kernel_>
const typename Predicate_h_3<Kernel_>::Predicate_h_3 *Predicate_h_3<Kernel_>::sub_predicates() const
{
    return this;
}

template<class Kernel_>
bool Predicate_h_3<Kernel_>::evaluate(const bool *signs) const
{
    // evaluate collision predicate
    return signs[0];
}
} // namespace CS
