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
#ifndef LIBCS_SPIN_QUADRIC_MESH_3_H
#define LIBCS_SPIN_QUADRIC_MESH_3_H

#include "Smooth_triangle_3.h"
#include "Spin_3.h"
#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/Complex_2_in_triangulation_3.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Implicit_surface_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/output_surface_facets_to_polyhedron.h>
#include <CGAL/IO/output_surface_facets_to_triangle_soup.h>
#include <boost/function.hpp>
#include <boost/bind.hpp>

namespace CS
{
// Spin_quadric_3:
//     mesher. s0 is dropped.
//
// This is a decorator. It does not take ownership of the object.
template<class Kernel_>
class Spin_quadric_mesh_3
{
public:
    typedef CGAL::Surface_mesh_default_triangulation_3                          Surface_mesh_triangulation_3;
    typedef CGAL::Complex_2_in_triangulation_3<Surface_mesh_triangulation_3>    Complex_2_in_triangulation_3;

    typedef Surface_mesh_triangulation_3::Geom_traits                           Geom_traits;

    typedef Geom_traits::Point_3                                                Point_3;

    typedef Geom_traits::K                                                      Geom_kernel;
    typedef Geom_traits::FT                                                     FT;
    typedef CS::Smooth_triangle_3<Geom_kernel>                                  Smooth_triangle_3;

private:
    typedef Geom_traits::Sphere_3                                               Sphere_3;
    typedef Geom_traits::Triangle_3                                             Triangle_3;

    typedef Geom_traits::Vector_3                                               Vector_3;
    typedef CS::Spin_3<FT>                                                      Spin_3;
    typedef CS::Diff_spin_3<FT>                                                 Diff_spin_3;

    typedef boost::function<FT (Point_3)>                                       Function;

    typedef CGAL::Implicit_surface_3<Geom_traits, Function>                     Implicit_surface_3;

    typedef Spin_quadric_mesh_3                                                 Self;

    FT m_a11, m_a12, m_a13, m_a14;
    FT m_a22, m_a23, m_a24;
    FT m_a33, m_a34;
    FT m_a44;

    FT evaluate_with_sign(Point_3 p, FT sign);

    // negative (left) signed set (leaf)
    FT evaluate_left(Point_3 p);

    // positive (right) signed set (leaf)
    FT evaluate_right(Point_3 p);

    Vector_3 evaluate_gradient(Point_3 p, FT sign);

public:
    typedef Kernel_                             R;
    typedef CGAL::Polyhedron_3<CGAL::Epick>     Polyhedron_3;

    /*
     * Construct a mesher for a spin quadric
     */
    Spin_quadric_mesh_3(const Spin_quadric_3<R> &spin_quadric);

    /*
     * Mesh to a polyhedron
     * Valid iff the spin quadric is a manifold
     */
    void mesh_polyhedron(Polyhedron_3 &polyhedron_left,
            Polyhedron_3 &polyhedron_right,
            double angular_bound = 30.0,
            double radius_bound = 0.1,
            double distance_bound = 0.1);

    /*
     * Mesh to a triangulation
     */
    void mesh_surface_triangulation(
            Surface_mesh_triangulation_3 &triangulation_left,
            Surface_mesh_triangulation_3 &triangulation_right,
            double angular_bound = 30.0,
            double radius_bound = 0.1,
            double distance_bound = 0.1);

    /*
     * Mesh to triangle soup
     *
     * @param OutputIterator value_type must be convertible from Smooth_triangle_3<Kernel>
     */
    template<typename OutputIterator>
    void mesh_triangle_soup(
            OutputIterator output_iterator_left,
            OutputIterator output_iterator_right,
            double angular_bound = 30.0,
            double radius_bound = 0.1,
            double distance_bound = 0.1);

private:
    template<typename OutputIterator>
    void mesh_internal_triangle_soup(
            OutputIterator output_iterator,
            Surface_mesh_triangulation_3 &triangulation,
            Complex_2_in_triangulation_3 &c2t3,
            const FT &sign);

    void mesh_internal_complex_in_triangulation(
            Complex_2_in_triangulation_3 &c2t3_left,
            Complex_2_in_triangulation_3 &c2t3_right,
            double angular_bound,
            double radius_bound,
            double distance_bound);

    const Spin_quadric_3<R> &       m_spin_quadric;
};
} // namespace CS

#include "Spin_quadric_mesh_3.ipp"

#endif // LIBCS_SPIN_QUADRIC_MESH_3_H
