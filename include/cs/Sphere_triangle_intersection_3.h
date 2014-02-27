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
#ifndef LIBCS_SPHERE_TRIANGLE_INTERSECTION_3_H
#define LIBCS_SPHERE_TRIANGLE_INTERSECTION_3_H

#include <CGAL/Kernel/global_functions.h>

namespace CS
{
template<class K>
typename K::RT dot_product(const CGAL::Vector_3<K> &v, const CGAL::Vector_3<K> &w);

template<class K>
bool sphere_triangle_intersection_3(const CGAL::Vector_3<K> &A_, const CGAL::Vector_3<K> &B_, const CGAL::Vector_3<K> &C_,
                                    const typename K::RT &rr, const CGAL::Vector_3<K> &P);

} // namespace CS

#include "Sphere_triangle_intersection_3.ipp"

#endif // LIBCS_SPHERE_TRIANGLE_INTERSECTION_3_H
