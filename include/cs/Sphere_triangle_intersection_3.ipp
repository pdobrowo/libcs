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
#include "Sphere_triangle_intersection_3.h"

namespace CS
{
template<class Kernel_>
typename Kernel_::RT dot_product(const CGAL::Vector_3<Kernel_> &v, const CGAL::Vector_3<Kernel_> &w)
{
    return v * w;
}

template<class Kernel_>
bool sphere_triangle_intersection_3(const CGAL::Vector_3<Kernel_> &A_, const CGAL::Vector_3<Kernel_> &B_, const CGAL::Vector_3<Kernel_> &C_,
                                    const typename Kernel_::RT &rr, const CGAL::Vector_3<Kernel_> &P)
{
    typedef typename Kernel_::RT      RT;
    typedef CGAL::Vector_3<Kernel_>   Vector_3;

    // Taken from: http://realtimecollisiondetection.net/blog/?p=103
    //    by Christer Ericson, December 30, 2010
    //
    // Adaptated and generalized
    Vector_3 A = A_ - P;
    Vector_3 B = B_ - P;
    Vector_3 C = C_ - P;
    Vector_3 V = CGAL::cross_product(B - A, C - A);
    RT d = dot_product(A, V);
    RT e = dot_product(V, V);
    bool sep1 = d * d > rr * e;
    RT aa = dot_product(A, A);
    RT ab = dot_product(A, B);
    RT ac = dot_product(A, C);
    RT bb = dot_product(B, B);
    RT bc = dot_product(B, C);
    RT cc = dot_product(C, C);
    bool sep2 = (aa > rr) && (ab > aa) && (ac > aa);
    bool sep3 = (bb > rr) && (ab > bb) && (bc > bb);
    bool sep4 = (cc > rr) && (ac > cc) && (bc > cc);
    Vector_3 AB = B - A;
    Vector_3 BC = C - B;
    Vector_3 CA = A - C;
    RT d1 = ab - aa;
    RT d2 = bc - bb;
    RT d3 = ac - cc;
    RT e1 = dot_product(AB, AB);
    RT e2 = dot_product(BC, BC);
    RT e3 = dot_product(CA, CA);
    Vector_3 Q1 = A * e1 - d1 * AB;
    Vector_3 Q2 = B * e2 - d2 * BC;
    Vector_3 Q3 = C * e3 - d3 * CA;
    Vector_3 QC = C * e1 - Q1;
    Vector_3 QA = A * e2 - Q2;
    Vector_3 QB = B * e3 - Q3;
    bool sep5 = (dot_product(Q1, Q1) > rr * e1 * e1) && (dot_product(Q1, QC) > 0);
    bool sep6 = (dot_product(Q2, Q2) > rr * e2 * e2) && (dot_product(Q2, QA) > 0);
    bool sep7 = (dot_product(Q3, Q3) > rr * e3 * e3) && (dot_product(Q3, QB) > 0);
    bool separated = sep1 || sep2 || sep3 || sep4 || sep5 || sep6 || sep7;
    return !separated;
}
} // namespace CS
