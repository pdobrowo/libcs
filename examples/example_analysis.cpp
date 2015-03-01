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
#include "example_analysis.h"
#include <cs/Spin_inexact_kernel_3.h>
#include <CGAL/Cartesian.h>
#include <eigen3/Eigen/Dense>
#include <iostream>
#include <cmath>

// kernel settings
typedef double RT;
typedef CGAL::Cartesian<RT> Kernel_base;
typedef CS::Spin_inexact_kernel_3<Kernel_base> Kernel;

// scene type
typedef Kernel::Predicate_tt_3 Predicate_tt_3;
typedef Kernel::Spin_null_configuration_space_3<Predicate_tt_3>::Type Configuration_space_3;
typedef Configuration_space_3::Parameters Parameters;

// objects
typedef Kernel::Point_3 Point_3;
typedef Kernel::Vector_3 Vector_3;
typedef Kernel::Triangle_3 Triangle_3;
typedef Kernel::Spin_quadric_3 Spin_quadric_3;
typedef Kernel::Predicate_g_3 Predicate_g_3;

void example_analysis()
{
    // 4-sided rectangular pyramid
    const Point_3 knob_vertices[5] = {
        Point_3(0, 0, 0),
        Point_3(2, 1, 8),
        Point_3(-2, 1, 8),
        Point_3(-2, -1, 8),
        Point_3(2, -1, 8)
    };

    const Triangle_3 knob_triangles[6] = {
        Triangle_3(knob_vertices[0], knob_vertices[2], knob_vertices[1]),
        Triangle_3(knob_vertices[0], knob_vertices[3], knob_vertices[2]),
        Triangle_3(knob_vertices[0], knob_vertices[4], knob_vertices[3]),
        Triangle_3(knob_vertices[0], knob_vertices[1], knob_vertices[4]),
        Triangle_3(knob_vertices[1], knob_vertices[2], knob_vertices[3]),
        Triangle_3(knob_vertices[3], knob_vertices[4], knob_vertices[1])
    };

    // rectangular stamp with tooth
    const Point_3 body_vertices[12] = {
        Point_3(7, -3, 5),
        Point_3(9, -3, 5),
        Point_3(7, 6, 5),
        Point_3(7, 3, 5),
        Point_3(-9, 3, 5),
        Point_3(-7, 3, 5),
        Point_3(-7, -6, 5),
        Point_3(-7, -3, 5),
        Point_3(2, 3, 5),
        Point_3(-2, 3, 5),
        Point_3(-2, -1, 5),
        Point_3(2, -1, 5)
    };

    const Triangle_3 body_triangles[6] = {
        Triangle_3(body_vertices[0], body_vertices[1], body_vertices[2]),
        Triangle_3(body_vertices[3], body_vertices[2], body_vertices[4]),
        Triangle_3(body_vertices[5], body_vertices[4], body_vertices[6]),
        Triangle_3(body_vertices[7], body_vertices[6], body_vertices[1]),
        Triangle_3(body_vertices[8], body_vertices[9], body_vertices[10]),
        Triangle_3(body_vertices[10], body_vertices[11], body_vertices[8])
    };

    // configuration space
    Configuration_space_3 cs;

    cs.create_from_scene(knob_triangles, knob_triangles + 6,
                         body_triangles, body_triangles + 6,
                         Parameters());

    // analysis

    // print spin-surface associated ellipsoids
    std::vector<Spin_quadric_3>::const_iterator spin_quadrics_iterator = cs.spin_quadrics_begin();
    std::vector<Predicate_g_3>::const_iterator general_predicate_iterator = cs.general_predicates_begin();
    assert(cs.size_of_spin_quadrics() == cs.size_of_general_predicates());

    typedef Configuration_space_3::General_predicate_size_type Count_type;
    Count_type count = cs.size_of_general_predicates();

    for (Count_type idx = 0; idx < count; ++idx)
    {
        // current
        const Spin_quadric_3 &spin_quadric = *spin_quadrics_iterator;
        const Predicate_g_3 &general_predicate = *general_predicate_iterator;

        std::cout << "--------------------------------------------------------" << std::endl;

        // matrix
        Kernel::Matrix ellipsoid_matrix = spin_quadric.matrix();

        Eigen::Matrix4d eigen_matrix;

        for (int row = 0; row < 4; ++row)
            for (int column = 0; column < 4; ++column)
                eigen_matrix(row, column) = ellipsoid_matrix.get(row, column);

        std::cout << "spin-quadric associated matrix:" << std::endl << eigen_matrix << std::endl;

        // eigensystem
        Eigen::EigenSolver<Eigen::Matrix4d> eigensolver(eigen_matrix);

        if (eigensolver.info() == Eigen::Success)
        {
            std::cout << "eigenvalues:" << std::endl << eigensolver.eigenvalues() << std::endl;
            std::cout << "eigenvectors:" << std::endl << eigensolver.eigenvectors() << std::endl;
        }
        else
            std::cout << "failed to calculate eigenvalues!" << std::endl;

        // analysis
        Vector_3 k = general_predicate.k();
        Vector_3 l = general_predicate.l();
        Vector_3 a = general_predicate.a();
        Vector_3 b = general_predicate.b();
        RT c = general_predicate.c();

        Vector_3 p = CGAL::cross_product(k, l);
        Vector_3 q = a - b;
        Vector_3 u = k - l;
        Vector_3 v = CGAL::cross_product(a, b);

        RT pp = p.squared_length();
        RT qq = q.squared_length();
        RT uu = u.squared_length();
        RT vv = v.squared_length();

        RT sqrt_ppqq = std::sqrt(pp * qq);
        RT sqrt_uuvv = std::sqrt(uu * vv);

        RT l1 = c - (+ sqrt_ppqq + sqrt_uuvv);
        RT l2 = c - (+ sqrt_ppqq - sqrt_uuvv);
        RT l3 = c - (- sqrt_ppqq + sqrt_uuvv);
        RT l4 = c - (- sqrt_ppqq - sqrt_uuvv);

        std::cout << "analysis values:" << std::endl << l1 << std::endl << l2 << std::endl << l3 << std::endl << l4 << std::endl;

        // Z and z

        // next spin quadrics and general predicate
        ++spin_quadrics_iterator;
        ++general_predicate_iterator;
    }
}
