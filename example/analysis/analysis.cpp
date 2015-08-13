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
#include "analysis.h"
#include <cs/Spin_inexact_kernel_3.h>
#include <cs/Spin_3.h>
#include <CGAL/Cartesian.h>
#include <eigen3/Eigen/Dense>
#include <cstdlib>
#include <iostream>
#include <iomanip>
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

typedef CS::Spin_3<double> Spin;

Spin calculate_analysis_vector(const Vector_3 &p, const Vector_3 &q, const Vector_3 &u, const Vector_3 &v, RT alpha, RT beta)
{
    Vector_3 pxq = CGAL::cross_product(p, q);
    Vector_3 uxv = CGAL::cross_product(u, v);
    Vector_3 pxv = CGAL::cross_product(p, v);
    Vector_3 qxu = CGAL::cross_product(q, u);
    Vector_3 r = pxq + uxv;
    RT trr = r.x() + r.y() + r.z();
    RT rr = r * r;
    RT pq = p * q;
    RT uv = u * v;
    RT qu = q * u;
    RT pv = p * v;
    RT pp = p * p;
    RT qq = q * q;
    RT uu = u * u;
    RT vv = v * v;
    RT trp = p.x() + p.y() + p.z();
    RT trq = q.x() + q.y() + q.z();
    RT tru = u.x() + u.y() + u.z();
    RT trv = v.x() + v.y() + v.z();
    Vector_3 pxu = CGAL::cross_product(p, u);
    Vector_3 qxv = CGAL::cross_product(q, v);
    RT trpxu = pxu.x() + pxu.y() + pxu.z();
    RT trqxv = qxv.x() + qxv.y() + qxv.z();
    RT l = alpha * std::sqrt(pp * qq) + beta * std::sqrt(uu * vv);
    RT ll = l * l;
    RT a44 = pq + uv + l;
    
    Vector_3 big_z = a44 * (qxv * trpxu
                            + pxu * trqxv
                            + p * qq * trp
                            + q * pp * trq
                            + u * vv * tru
                            + v * uu * trv
                            - l * (p * trq
                                   + q * trp
                                   + u * trv
                                   + v * tru
                                  )
                           )
                     - (p * tru + u * trp) * (pv * qq + qu * vv)
                     - (q * trv + v * trq) * (pv * uu + qu * pp)
                     + (p * trq + q * trp) * (pp * qq + pv * qu - pq * uv)
                     + (u * trv + v * tru) * (uu * vv + pv * qu - pq * uv)
                     + p * uv * trp * qq
                     + q * uv * trq * pp
                     - u * uv * tru * vv
                     - v * uv * trv * uu
                     - p * pq * trp * qq
                     - q * pq * trq * pp
                     + u * pq * tru * vv
                     + v * pq * trv * uu
                     + l * trr * r;
  
    RT small_z = - l * rr + a44 * (ll - pp * qq - uu * vv);

    RT e1 = a44 * (big_z.z() + small_z);
    RT e2 = a44 * (big_z.x() + small_z);
    RT e3 = a44 * (big_z.y() + small_z);
    RT e4 = -big_z * r - small_z * trr;

    Vector_3 vs_e4 = pxq * (+ 4 * uv * pv * qu)
                     + uxv * (+ 4 * pq * pv * qu)
                     + pxu * (- 2 * qq * uv * pv + 2 * vv * pq * qu)
                     + qxv * (+ 2 * pp * uv * qu - 2 * uu * pq * pv)
                     + pxv * (- 2 * qu * pq * uv + 2 * pv * (qq * uu - qu * qu))
                     + qxu * (+ 2 * pv * pq * uv - 2 * qu * (pp * vv - pv * pv))
                     + (pq + uv + l) * (
                         + 2 * (uu * vv * pxq + pp * qq * uxv)
                         + (pxq + uxv) * (l * l - pp * qq - uu * vv)
                         - l * ( (uxv * q) * p + (uxv * p) * q + (pxq * v) * u + (pxq * u) * v) );

    Vector_3 test = pxq * (+ 4 * uv * pv * qu)
            + uxv * (+ 4 * pq * pv * qu)
            + pxu * (- 2 * qq * uv * pv + 2 * vv * pq * qu)
            + qxv * (+ 2 * pp * uv * qu - 2 * uu * pq * pv)
            + pxv * (- 2 * qu * pq * uv + 2 * pv * (qq * uu - qu * qu))
            + qxu * (+ 2 * pv * pq * uv - 2 * qu * (pp * vv - pv * pv));

    std::cout << "$$$ TEST: " << test << std::endl;

    RT s_e4 = - vs_e4.x() - vs_e4.y() - vs_e4.z();
    std::cout << std::setprecision(10) << std::fixed << "e4: " << e4 << ", s_e4: " << s_e4 << std::endl;

    RT serr = e4 - s_e4;

    if (std::fabs(serr) > 10e-3)
    {
        std::cout << "MISMATCH S/E4; serr = " << serr << std::endl;
        exit(0);
    }

    RT ee = e1 * e1 + e2 * e2 + e3 * e3 + e4 * e4;
    RT en = std::sqrt(ee);
    
    std::cout << "norm: " << en << std::endl;

    if (std::fabs(en) < 0.0000001)
    {
        std::cout << "ZERO NORM" << std::endl;
        return Spin();
    }
    
    return Spin(e1 / en, e2 / en, e3 / en, e4 / en);
}

void analysis()
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
        Kernel::Matrix matrix = spin_quadric.matrix();
        Eigen::Matrix4d eigen_matrix;

        for (int row = 0; row < 4; ++row)
            for (int column = 0; column < 4; ++column)
                eigen_matrix(row, column) = matrix.get(row, column);

        //std::cout << "spin-quadric associated matrix:" << std::endl << eigen_matrix << std::endl;

        // eigensystem
        Eigen::EigenSolver<Eigen::Matrix4d> eigensolver(eigen_matrix);

        if (eigensolver.info() == Eigen::Success)
        {
            std::cout << "eigen values:" << std::endl
                << eigensolver.eigenvalues()[0].real() << std::endl
                << eigensolver.eigenvalues()[1].real() << std::endl
                << eigensolver.eigenvalues()[2].real() << std::endl
                << eigensolver.eigenvalues()[3].real() << std::endl;
                
            std::cout << "eigen vector 1: [" << eigensolver.eigenvectors()(0, 0).real() << ", " << eigensolver.eigenvectors()(1, 0).real() << ", " << eigensolver.eigenvectors()(2, 0).real() << ", " << eigensolver.eigenvectors()(3, 0).real() << "]" << std::endl;
            std::cout << "eigen vector 2: [" << eigensolver.eigenvectors()(0, 1).real() << ", " << eigensolver.eigenvectors()(1, 1).real() << ", " << eigensolver.eigenvectors()(2, 1).real() << ", " << eigensolver.eigenvectors()(3, 1).real() << "]" << std::endl;
            std::cout << "eigen vector 3: [" << eigensolver.eigenvectors()(0, 2).real() << ", " << eigensolver.eigenvectors()(1, 2).real() << ", " << eigensolver.eigenvectors()(2, 2).real() << ", " << eigensolver.eigenvectors()(3, 2).real() << "]" << std::endl;
            std::cout << "eigen vector 4: [" << eigensolver.eigenvectors()(0, 3).real() << ", " << eigensolver.eigenvectors()(1, 3).real() << ", " << eigensolver.eigenvectors()(2, 3).real() << ", " << eigensolver.eigenvectors()(3, 3).real() << "]" << std::endl;
        }
        else
            std::cout << "failed to calculate eigenvalues!" << std::endl;

        std::cout << std::endl;
        
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
        Spin e1 = calculate_analysis_vector(p, q, u, v, +1, +1);
        Spin e2 = calculate_analysis_vector(p, q, u, v, +1, -1);
        Spin e3 = calculate_analysis_vector(p, q, u, v, -1, +1);
        Spin e4 = calculate_analysis_vector(p, q, u, v, -1, -1);
        
        std::cout << "analysis vector 1: [" << e1.s12() << ", " << e1.s23() << ", " << e1.s31() << ", " << e1.s0() << "]" << std::endl;
        std::cout << "analysis vector 2: [" << e2.s12() << ", " << e2.s23() << ", " << e2.s31() << ", " << e2.s0() << "]" << std::endl;
        std::cout << "analysis vector 3: [" << e3.s12() << ", " << e3.s23() << ", " << e3.s31() << ", " << e3.s0() << "]" << std::endl;
        std::cout << "analysis vector 4: [" << e4.s12() << ", " << e4.s23() << ", " << e4.s31() << ", " << e4.s0() << "]" << std::endl;

        // next spin quadrics and general predicate
        ++spin_quadrics_iterator;
        ++general_predicate_iterator;
    }
}
