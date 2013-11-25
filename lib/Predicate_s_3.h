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
#ifndef LIBCS_PREDICATE_S_3_H
#define LIBCS_PREDICATE_S_3_H

#include <ostream>

namespace CS
{
// S3: (k x l) * rot(a - b) + (k - l) * rot(a x b)
template<class Kernel_>
class Predicate_s_3
{
    typedef typename Kernel_::RT         RT;
    typedef typename Kernel_::Vector_3   Vector_3;

public:
    typedef Kernel_                      Kernel;

    Predicate_s_3();

    Predicate_s_3(const Vector_3 &k,
       const Vector_3 &l,
       const Vector_3 &a,
       const Vector_3 &b);

    const Vector_3 &k() const;
    const Vector_3 &l() const;
    const Vector_3 &a() const;
    const Vector_3 &b() const;

    Predicate_s_3 opposite() const;

private:
    Vector_3    m_k;
    Vector_3    m_l;
    Vector_3    m_a;
    Vector_3    m_b;
};
} // namespace CS

#include "Predicate_s_3.ipp"

#endif // LIBCS_PREDICATE_S_3_H
