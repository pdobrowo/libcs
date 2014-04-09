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
#include <cs/Spin_inexact_kernel_3.h>
#include <cs/Spin_qsic_trim_3.h>
#include <CGAL/Cartesian.h>
#include <log4cxx/logger.h>

//#define DISPLAY_SPIN_SURFACE_ASSOCIATED_ELLIPSOID

#ifdef DISPLAY_SPIN_SURFACE_ASSOCIATED_ELLIPSOID
    #include <iostream>
    #include <eigen3/Eigen/Dense>
#endif // DISPLAY_SPIN_SURFACE_ASSOCIATED_ELLIPSOID

typedef CS::Default_inexact_kernel              Kernel;
typedef Kernel::RT                              RT;

typedef Kernel::Predicate_tt_3                  Predicate_tt_3;
typedef Kernel::Predicate_bb_3                  Predicate_bb_3;

typedef Kernel::Predicate_g_3                   Predicate_g_3;

typedef Kernel::Vector_3                        Vector_3;
typedef Kernel::Point_3                         Point_3;

typedef Kernel::Ball_3                          Ball_3;
typedef Kernel::Triangle_3                      Triangle_3;

typedef Kernel::Spin_cell_configuration_space_3<Predicate_tt_3>::Type CS_TT_C;
typedef Kernel::Spin_cell_configuration_space_3<Predicate_bb_3>::Type CS_BB_C;
typedef Kernel::Spin_raster_configuration_space_3<Predicate_tt_3>::Type CS_TT_R;
typedef Kernel::Spin_raster_configuration_space_3<Predicate_bb_3>::Type CS_BB_R;

typedef CS::Spin_qsic_trim_3                    Spin_qsic_trim_3;

typedef Kernel::Spin_reduced_quadric_3          Spin_reduced_quadric_3;

static void test_inexact_perf_single(int number_of_predicates)
{
    unsigned long long start = CS::get_tick_count();

#if 0

    Ball_3 rotating_ball(Vector_3(-5, 0, 0), 5);
    std::vector<Ball_3> obstacle_balls;

    for (int p = 0; p < number_of_predicates; ++p)
    {
        double r = CS::random_double(1, 5);
        double x = CS::random_double(-5, 5);
        double y = CS::random_double(-5, 5);
        double z = CS::random_double(-5, 5);

        obstacle_balls.push_back(Kernel::Ball_3(Kernel::Vector_3(x, y, z), r));
    }

    typedef CS_BB_C ConfigurationSpace;
    ConfigurationSpace cs;

    cs.create_from_scene(&rotating_ball, 1 + &rotating_ball,
                         obstacle_balls.begin(), obstacle_balls.end());

#else

    Triangle_3 rotating_triangle(Point_3(0, 0, 0),
                                 Point_3(5, 0, 0),
                                 Point_3(0, 5, 0));

    Triangle_3 stationary_triangle(Point_3(0, 0, 1),
                                   Point_3(0, 0, 4),
                                   Point_3(4, 0, 4));

    typedef CS_TT_C ConfigurationSpace;
    ConfigurationSpace cs;

    cs.create_from_scene(&rotating_triangle, 1 + &rotating_triangle,
                         &stationary_triangle, 1 + &stationary_triangle);

#endif

    unsigned long long end = CS::get_tick_count();
    unsigned long long delta = end - start;

    std::cout << number_of_predicates << " complex predicate(s) in " << delta << " ms" << std::endl;

    typedef ConfigurationSpace::Spin_quadric_const_iterator Spin_quadric_const_iterator;
    typedef Kernel::Matrix Matrix;

    typedef ConfigurationSpace::Generic_predicate_const_iterator Generic_predicate_const_iterator;

#ifdef DISPLAY_SPIN_SURFACE_ASSOCIATED_ELLIPSOID

    // print spin-surface associated ellipsoids
    for (Spin_quadric_const_iterator it = cs.spin_quadrics_begin(); it != cs.spin_quadrics_end(); ++it)
    {
        Matrix matrix = it->matrix();
        Eigen::Matrix4d eigen_matrix;

        for (int row = 0; row < 4; ++row)
            for (int column = 0; column < 4; ++column)
                eigen_matrix(row, column) = matrix.get(row, column);

        //std::cout << "spin-quadric associated matrix:" << std::endl << eigen_matrix << std::endl;

        Eigen::EigenSolver<Eigen::Matrix4d> eigensolver(eigen_matrix);

        if (eigensolver.info() == Eigen::Success)
            std::cout << "eigenvalues:" << std::endl << eigensolver.eigenvalues() << std::endl;
        else
            std::cout << "failed to calculate eigenvalues!" << std::endl;
    }


    // print predicates
    for (Generic_predicate_const_iterator it = cs.generic_predicates_begin(); it != cs.generic_predicates_end(); ++it)
    {
        Predicate_g_3 gp = *it;
        //std::cout << "generic predicate: " << generic_predicate << std::endl;

        Vector_3 k = gp.k();
        Vector_3 l = gp.l();
        Vector_3 a = gp.a();
        Vector_3 b = gp.b();
        RT c = gp.c();

        RT p = std::sqrt(CGAL::cross_product(k, l).squared_length() * (a - b).squared_length());
        RT q = std::sqrt(CGAL::cross_product(a, b).squared_length() * (k - l).squared_length());

        std::cout << "p: " << p << ", q: " << q << ", c: " << c << std::endl;

        RT e1 = c - p - q;
        RT e2 = c - p + q;
        RT e3 = c + p - q;
        RT e4 = c + p + q;

        std::cout << "values: " << e1 << ", " << e2 << ", " << e3 << ", " << e4 << std::endl;
    }

#endif // DISPLAY_SPIN_SURFACE_ASSOCIATED_ELLIPSOID

    for (Spin_quadric_const_iterator it = cs.spin_quadrics_begin(); it != cs.spin_quadrics_end(); ++it)
    {
        for (Spin_quadric_const_iterator it2 = it; it2 != cs.spin_quadrics_end(); ++it2)
        {
            Spin_qsic_trim_3 qsic_trim = Spin_qsic_trim_3(Spin_reduced_quadric_3(*it), Spin_reduced_quadric_3(*it2));

            double s12, s23, s31;

            if (!qsic_trim.begin(s12, s23, s31))
            {
                std::cout << "FAILURE!!!" << std::endl;
                continue;
            }

            for (int i = 0; i<1; ++i)
            {
                std::cout << "s12: " << s12 << ", s23: " << s23 << ", s31: " << s31 << std::endl;

                if (!qsic_trim.next(s12, s23, s31))
                {
                    std::cout << "FAILURE!!!" << std::endl;
                    break;
                }
            }

        }
    }
}

void test_inexact_perf()
{
    for (int p = 1; p <= 1; ++p)
        test_inexact_perf_single(1);
}
