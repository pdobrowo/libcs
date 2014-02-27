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
#include <cs/Spin_kernel_3.h>
#include <cs/Spin_configuration_space_3.h>
#include <cs/Spin_cell_graph_3.h>
#include <cs/Spin_exact_graph_3.h>
#include <cs/Spin_raster_graph_3.h>
#include <cs/Predicate_tt_3.h>
#include <cs/Predicate_bb_3.h>
#include <cs/Benchmark.h>
#include <cs/Bigint.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Filtered_kernel.h>
#include <log4cxx/logger.h>

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("Test.test_exact_perf"));

typedef CGAL::Bigint                                RT;
typedef CGAL::Filtered_kernel<CGAL::Cartesian<RT> > Base_kernel;
typedef CS::Spin_kernel_3<Base_kernel>              Kernel;

typedef Kernel::Predicate_tt_3 Predicate_tt_3;
typedef Kernel::Predicate_bb_3 Predicate_bb_3;

// cell
typedef Kernel::Spin_cell_configuration_space_3<Predicate_tt_3>::Type CS_TT_C;
typedef Kernel::Spin_cell_configuration_space_3<Predicate_bb_3>::Type CS_BB_C;

// raster
typedef Kernel::Spin_raster_configuration_space_3<Predicate_tt_3>::Type CS_TT_R;
typedef Kernel::Spin_raster_configuration_space_3<Predicate_bb_3>::Type CS_BB_R;

// exact
typedef Kernel::Spin_exact_configuration_space_3<Predicate_tt_3>::Type CS_TT_E;
typedef Kernel::Spin_exact_configuration_space_3<Predicate_bb_3>::Type CS_BB_E;

static void test_exact_perf_single(int number_of_predicates)
{
    unsigned long long start = CS::get_tick_count();

    Kernel::Ball_3 rotating_ball(Kernel::Vector_3(-50, 0, 0), 50);
    std::vector<Kernel::Ball_3> obstacle_balls;

    for (int p = 0; p < number_of_predicates; ++p)
    {
        int r = CS::random_int(1, 25);
        int x = CS::random_int(-50, 50);
        int y = CS::random_int(-50, 50);
        int z = CS::random_int(-50, 50);

        obstacle_balls.push_back(Kernel::Ball_3(Kernel::Vector_3(x, y, z), r));
    }

    CS_BB_E cs;
    cs.create_from_scene(&rotating_ball, 1 + &rotating_ball,
                         obstacle_balls.begin(), obstacle_balls.end(),
                         CS_BB_E::Parameters());

    unsigned long long end = CS::get_tick_count();
    unsigned long long delta = end - start;

    printf("%i %llu\n", number_of_predicates, delta);
}

void test_exact_perf()
{
    for (int p = 1; p <= 20; ++p)
        test_exact_perf_single(p);

    //for (int p = 10; p <= 100; p += 10)
    //    test_perfornamce_single(p);
}
