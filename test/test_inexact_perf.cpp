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
#include <CGAL/Cartesian.h>
#include <log4cxx/logger.h>

#define DISPLAY_SPIN_SURFACE_ASSOCIATED_ELLIPSOID

#ifdef DISPLAY_SPIN_SURFACE_ASSOCIATED_ELLIPSOID
    #include <iostream>
    #include <eigen3/Eigen/Dense>
#endif // DISPLAY_SPIN_SURFACE_ASSOCIATED_ELLIPSOID

typedef double                                  RT;
typedef CGAL::Cartesian<RT>                     Base_kernel;
typedef CS::Spin_inexact_kernel_3<Base_kernel>  Kernel;

typedef Kernel::Predicate_tt_3 Predicate_tt_3;
typedef Kernel::Predicate_bb_3 Predicate_bb_3;

typedef Kernel::Predicate_g_3 Predicate_g_3;

typedef Kernel::Spin_cell_configuration_space_3<Predicate_tt_3>::Type CS_TT_C;
typedef Kernel::Spin_cell_configuration_space_3<Predicate_bb_3>::Type CS_BB_C;
typedef Kernel::Spin_raster_configuration_space_3<Predicate_tt_3>::Type CS_TT_R;
typedef Kernel::Spin_raster_configuration_space_3<Predicate_bb_3>::Type CS_BB_R;

static void test_inexact_perf_single(int number_of_predicates)
{
    unsigned long long start = CS::get_tick_count();

    Kernel::Ball_3 rotating_ball(Kernel::Vector_3(-50, 0, 0), 50);
    std::vector<Kernel::Ball_3> obstacle_balls;

    for (int p = 0; p < number_of_predicates; ++p)
    {
        double r = CS::random_double(1, 3);
        double x = CS::random_double(-5, 5);
        double y = CS::random_double(-5, 5);
        double z = CS::random_double(-5, 5);

        obstacle_balls.push_back(Kernel::Ball_3(Kernel::Vector_3(x, y, z), r));
    }

    typedef CS_BB_C ConfigurationSpace;

    ConfigurationSpace cs;
    cs.create_from_scene(&rotating_ball, 1 + &rotating_ball,
                         obstacle_balls.begin(), obstacle_balls.end());

    unsigned long long end = CS::get_tick_count();
    unsigned long long delta = end - start;

    printf("%i %llu\n", number_of_predicates, delta);

#ifdef DISPLAY_SPIN_SURFACE_ASSOCIATED_ELLIPSOID

    typedef ConfigurationSpace::Spin_quadric_const_iterator Spin_quadric_const_iterator;
    typedef Kernel::Matrix Matrix;

    // print spin-surface associated ellipsoids
    for (Spin_quadric_const_iterator it = cs.spin_quadrics_begin(); it != cs.spin_quadrics_end(); ++it)
    {
        Matrix matrix = it->matrix();
        Eigen::Matrix4d eigen_matrix;

        for (int row = 0; row < 4; ++row)
            for (int column = 0; column < 4; ++column)
                eigen_matrix(row, column) = matrix.get(row, column);

        std::cout << "spin-quadric associated matrix:" << std::endl << eigen_matrix << std::endl;

        Eigen::EigenSolver<Eigen::Matrix4d> eigensolver(eigen_matrix);

        if (eigensolver.info() == Eigen::Success)
            std::cout << "eigenvalues:" << std::endl << eigensolver.eigenvalues() << std::endl;
        else
            std::cout << "failed to calculate eigenvalues!" << std::endl;
    }

    typedef ConfigurationSpace::Generic_predicate_const_iterator Generic_predicate_const_iterator;

    // print predicates
    for (Generic_predicate_const_iterator it = cs.generic_predicates_begin(); it != cs.generic_predicates_end(); ++it)
    {
        Predicate_g_3 generic_predicate = *it;
        std::cout << "generic predicate:" << generic_predicate << std::endl;
    }

#endif // DISPLAY_SPIN_SURFACE_ASSOCIATED_ELLIPSOID

}

void test_inexact_perf()
{
    for (int p = 1; p <= 3; ++p)
        test_inexact_perf_single(p);
}
