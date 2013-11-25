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
#include "Smooth_triangle_3.h"

namespace CS
{
template <class R_>
Smooth_triangle_3<R_>::Smooth_triangle_3()
{
}

template <class R_>
Smooth_triangle_3<R_>::Smooth_triangle_3(const CGAL::Triangle_3<R_> &triangle)
    : m_triangle(triangle)
{
}

template <class R_>
Smooth_triangle_3<R_>::Smooth_triangle_3(
        const CGAL::Point_3<R_> &v0,
        const CGAL::Point_3<R_> &v1,
        const CGAL::Point_3<R_> &v2)
    : m_triangle(CGAL::Triangle_3<R_>(v0, v1, v2))
{
}

template <class R_>
Smooth_triangle_3<R_>::Smooth_triangle_3(
        const CGAL::Triangle_3<R_> &triangle,
        const CGAL::Vector_3<R_> &normal_0,
        const CGAL::Vector_3<R_> &normal_1,
        const CGAL::Vector_3<R_> &normal_2)
    : m_triangle(triangle),
      m_normal_0(normal_0),
      m_normal_1(normal_1),
      m_normal_2(normal_2)
{
}

template <class R_>
CGAL::Point_3<R_> Smooth_triangle_3<R_>::vertex(int i) const
{
    return m_triangle.vertex(i);
}

template <class R_>
CGAL::Point_3<R_> Smooth_triangle_3<R_>::operator[](int i) const
{
    return vertex(i);
}

template <class R_>
const CGAL::Triangle_3<R_> &Smooth_triangle_3<R_>::triangle() const
{
    return m_triangle;
}

template <class R_>
const CGAL::Vector_3<R_> &Smooth_triangle_3<R_>::normal_0() const
{
    return m_normal_0;
}

template <class R_>
const CGAL::Vector_3<R_> &Smooth_triangle_3<R_>::normal_1() const
{
    return m_normal_1;
}

template <class R_>
const CGAL::Vector_3<R_> &Smooth_triangle_3<R_>::normal_2() const
{
    return m_normal_2;
}
} // namespace CS
