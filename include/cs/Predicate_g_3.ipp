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
#include "Predicate_g_3.h"

namespace CS
{
template<class R>
Predicate_g_3<R>::Predicate_g_3()
    : m_k(Vector_3(RT(0), RT(0), RT(0))),
      m_l(Vector_3(RT(0), RT(0), RT(0))),
      m_a(Vector_3(RT(0), RT(0), RT(0))),
      m_b(Vector_3(RT(0), RT(0), RT(0))),
      m_c(RT(0))
{
}

template<class R>
Predicate_g_3<R>::Predicate_g_3(
        const Vector_3  &k,
        const Vector_3  &l,
        const Vector_3  &a,
        const Vector_3  &b,
        const RT &c)
    : m_k(k),
      m_l(l),
      m_a(a),
      m_b(b),
      m_c(c)
{
}

template<class R>
Predicate_g_3<R>::Predicate_g_3(const Predicate_h_3 &h3)
{
    Vector_3 v(h3.p().a(), h3.p().b(), h3.p().c());
    Vector_3 r = choose_r_vector(v);
    Vector_3 vr = CGAL::cross_product(v, r);
    RT d = h3.p().d();

    RT sq = vr.x() * vr.x() + vr.y() * vr.y() + vr.z() * vr.z();

    m_k = vr;
    m_l = CGAL::cross_product(v, vr);
    m_a = h3.b();
    m_b = -h3.b();
    m_c = RT(2) * d * sq;
}

template<class R>
Predicate_g_3<R>::Predicate_g_3(const Predicate_s_3 &s3)
    : m_k(s3.k()),
      m_l(s3.l()),
      m_a(s3.a()),
      m_b(s3.b()),
      m_c(0)
{
}

template<class R>
const typename Predicate_g_3<R>::Vector_3 &Predicate_g_3<R>::k() const
{
    return m_k;
}

template<class R>
const typename Predicate_g_3<R>::Vector_3 &Predicate_g_3<R>::l() const
{
    return m_l;
}

template<class R>
const typename Predicate_g_3<R>::Vector_3 &Predicate_g_3<R>::a() const
{
    return m_a;
}

template<class R>
const typename Predicate_g_3<R>::Vector_3 &Predicate_g_3<R>::b() const
{
    return m_b;
}

template<class R>
const typename Predicate_g_3<R>::RT &Predicate_g_3<R>::c() const
{
    return m_c;
}

template<class R>
Predicate_g_3<R> Predicate_g_3<R>::opposite() const
{
    // Swap k & l components or a & b components
    return Predicate_g_3<R>(l(), k(), a(), b());
}

template<class R>
std::ostream &operator <<(std::ostream &os, const Predicate_g_3<R> &predicate)
{
    return (os << "[" << predicate.k() << ";" << predicate.l() << ";"
                      << predicate.a() << ";" << predicate.b() << ";"
                      << predicate.c() << "]");
}
} // namespace CS

