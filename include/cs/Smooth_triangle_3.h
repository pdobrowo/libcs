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
#ifndef LIBCS_SMOOTH_TRIANGLE_3_H
#define LIBCS_SMOOTH_TRIANGLE_3_H

#include <CGAL/Triangle_3.h>
#include <CGAL/Vector_3.h>
#include <CGAL/Point_3.h>

namespace CS
{
template <class Kernel_>
class Smooth_triangle_3
{
public:
    typedef CGAL::Point_3<Kernel_> Point_3;
    typedef CGAL::Vector_3<Kernel_> Vector_3;
    typedef CGAL::Triangle_3<Kernel_> Triangle_3;

    Smooth_triangle_3();

    Smooth_triangle_3(const CGAL::Triangle_3<Kernel_> &triangle);

    Smooth_triangle_3(const CGAL::Point_3<Kernel_> &v0,
                      const CGAL::Point_3<Kernel_> &v1,
                      const CGAL::Point_3<Kernel_> &v2);

    Smooth_triangle_3(const CGAL::Triangle_3<Kernel_> &triangle,
                      const CGAL::Vector_3<Kernel_> &normal_0,
                      const CGAL::Vector_3<Kernel_> &normal_1,
                      const CGAL::Vector_3<Kernel_> &normal_2);

    const CGAL::Triangle_3<Kernel_>    &triangle() const;
    const CGAL::Vector_3<Kernel_>      &normal_0() const;
    const CGAL::Vector_3<Kernel_>      &normal_1() const;
    const CGAL::Vector_3<Kernel_>      &normal_2() const;

    CGAL::Point_3<Kernel_> vertex(int i) const;
    CGAL::Point_3<Kernel_> operator[](int i) const;

private:
    CGAL::Triangle_3<Kernel_>    m_triangle;
    CGAL::Vector_3<Kernel_>      m_normal_0;
    CGAL::Vector_3<Kernel_>      m_normal_1;
    CGAL::Vector_3<Kernel_>      m_normal_2;
};
} // namespace CS

#include "Smooth_triangle_3.ipp"

#endif // LIBCS_SMOOTH_TRIANGLE_3_H
