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
#ifndef LIBCS_PREDICATE_G_3_H
#define LIBCS_PREDICATE_G_3_H

#include <CGAL/Kernel/global_functions.h>
#include <ostream>

namespace CS
{
// G3: (k x l) * rot(a - b) + (k - l) * rot(a x b) + c
//
// general predicate
//
template<class Kernel_>
class Predicate_g_3
{
    typedef typename Kernel_::RT             RT;
    typedef typename Kernel_::Vector_3       Vector_3;

    typedef typename Kernel_::Predicate_h_3  Predicate_h_3;
    typedef typename Kernel_::Predicate_s_3  Predicate_s_3;

    Vector_3 choose_r_vector(const Vector_3 &v)
    {
        // choose a vector that is not parallel to v
        // project the vector on a plane and take one orthogonal to it
        if (v.x() != RT(0))
            return Vector_3(-v.y(), v.x(), v.z()); // take x-y plane
        else if (v.y() != RT(0))
            return Vector_3(v.x(), -v.z(), v.y()); // take y-z plane
        else
            return Vector_3(v.z(), v.y(), -v.x()); // take z-x plane
    }

public:
    typedef Kernel_                             Kernel;

    Predicate_g_3();

    Predicate_g_3(
        const Vector_3 &k,
        const Vector_3 &l,
        const Vector_3 &a,
        const Vector_3 &b,
        const RT &c);

    Predicate_g_3(const Predicate_h_3 &h3);

    Predicate_g_3(const Predicate_s_3 &s3);

    const Vector_3 &k() const;
    const Vector_3 &l() const;
    const Vector_3 &a() const;
    const Vector_3 &b() const;
    const RT &c() const;

    Predicate_g_3 opposite() const;

private:
    Vector_3    m_k;
    Vector_3    m_l;
    Vector_3    m_a;
    Vector_3    m_b;
    RT          m_c;
};
} // namespace CS

#include "Predicate_g_3.ipp"

#endif // LIBCS_PREDICATE_G_3_H
