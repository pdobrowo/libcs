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
#ifndef LIBCS_PREDICATE_H_3_H
#define LIBCS_PREDICATE_H_3_H

#include <ostream>

namespace CS
{
// H3: p.n * rot(b) + p.d
template<class Kernel_>
class Predicate_h_3
{
    typedef typename Kernel_::RT         RT;
    typedef typename Kernel_::Vector_3   Vector_3;
    typedef typename Kernel_::Plane_3    Plane_3;

public:
    typedef Kernel_                      Kernel;

    Predicate_h_3();

    Predicate_h_3(const Vector_3 &b,
                  const Plane_3 &p);

    const Vector_3 &b() const;
    const Plane_3 &p() const;

    Predicate_h_3 opposite() const;

private:
    Vector_3  m_b;  ///< Base vector
    Plane_3   m_p;  ///< Plane
};
} // namespace CS

#include "Predicate_h_3.ipp"

#endif // LIBCS_PREDICATE_H_3_H
