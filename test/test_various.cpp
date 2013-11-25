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
#include "Spin_kernel_3.h"
#include "Spin_qsic_mesh_3.h"
#include "Spin_qsip_mesh_3.h"
#include "Spin_quadric_mesh_3.h"
#include "Spin_quadric_tree_3.h"
#include "Spin_quadric_tree_generator.h"
#include "Spin_configuration_space_3.h"
#include "Spin_cell_graph_3.h"
#include "Build.h"
#include "Loader_scene.h"
#include "Loader_sphere_tree.h"
#include "Predicate_tt_3.h"
#include "Predicate_bb_3.h"
#include "Benchmark.h"
#include "Bigint.h"
#include "Uniform_random_spin_3.h"
#include <CGAL/Cartesian.h>
#include <CGAL/Filtered_kernel.h>
#include <CGAL/Root_of_traits.h>
#include <CGAL/Gmpfr.h>
#include <iostream>
#include <iomanip>
#include <string>
#include <boost/scoped_ptr.hpp>
#include <log4cxx/logger.h>
#include <log4cxx/basicconfigurator.h>
#include <log4cxx/propertyconfigurator.h>
#include <log4cxx/helpers/exception.h>

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("Test.test_various"));

#if 1

    // Exact kernel
    typedef CGAL::Bigint                                RT;
    typedef CGAL::Filtered_kernel<CGAL::Cartesian<RT> > Base_kernel;
    typedef CS::Spin_kernel_3<Base_kernel>              Kernel;

    // need to scale scene in double coordinates
    const double GLOBAL_SCENE_SCALE = 1000000.0;

#else

    // Floating point kernel
    #if 1
        typedef double                                  RT;
    #else
        typedef CGAL::Gmpfr                             RT;
        // default CGAL::Gmpfr precision
        CGAL::Gmpfr::Precision_type g_old_gmpfr_precision = CGAL::Gmpfr::set_default_precision(512);
    #endif

    typedef CGAL::Cartesian<RT>                     Kernel_base;
    typedef CS::Spin_inexact_kernel_3<Kernel_base>  Kernel;

    // no need to scale scene - it is already in double coordinates
    const double GLOBAL_SCENE_SCALE = 1.0;

#endif

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

void test_various(size_t param)
{
    (void)param;

    std::vector<Kernel::Triangle_3> robot_faces;
    std::vector<Kernel::Triangle_3> obstacle_faces;

#if 0
    // Scene from file
    const char *BENCHMARK_FILE = "benchmark/Test_Clip_Med";

    const Kernel::Vector_3 scaledRobotTranslation(RT(GLOBAL_SCENE_SCALE * 0.38),
                                                  RT(GLOBAL_SCENE_SCALE * 0.90),
                                                  RT(GLOBAL_SCENE_SCALE * 0.90));

    // WARN/TODO: This case is very difficult: it contains many degenerate TT3 predicates (parallel segments are involved)
//    const Kernel::Vector_3 scaledRobotTranslation(RT(GLOBAL_SCENE_SCALE * 0.00),
//                                                  RT(GLOBAL_SCENE_SCALE * 0.50),
//                                                  RT(GLOBAL_SCENE_SCALE * 0.90));

    if (!CS::load_scene(BENCHMARK_FILE, std::back_inserter(robot_faces), std::back_inserter(obstacle_faces), GLOBAL_SCENE_SCALE, scaledRobotTranslation))
    {
        LOG4CXX_FATAL(g_logger, "Cannot load benchmark file: " << BENCHMARK_FILE);
        return;
    }

#else
    // Scene defined
    robot_faces.push_back( // Robot should be placed near zero
        Kernel::Triangle_3(
            Kernel::Point_3(0, 0, 0),
            Kernel::Point_3(1, 0, 0),
            Kernel::Point_3(0, 0, 1)));

    obstacle_faces.push_back(
        Kernel::Triangle_3(
            Kernel::Point_3(-10, 0.99, -10),
            Kernel::Point_3(70, 0.99, -10),
            Kernel::Point_3(-10, 0.99, 70)));


    const double scale = 1.0;
    const Kernel::Vector_3 scaledRobotTranslation(scale * 0.0, scale * 0.0, scale * 0.0);

#endif

    LOG4CXX_INFO(logger, "------------------------------ TT");

#if 0
    {
        CS_TT_C cs;

        cs.create_from_scene(robot_faces.begin(), robot_faces.end(),
                             obstacle_faces.begin(), obstacle_faces.end(),
                             CS_TT_C::Parameters(param));

        CS_TT_C::Sample begin(0, 0, 0, 1);
        CS_TT_C::Sample end(1, 0, 0, 0);

        // find route
        cs.find_route(begin, end);
    }
#endif

#if 0
    {
        CS_TT_R cs;

        cs.create_from_scene(robot_faces.begin(), robot_faces.end(),
                             obstacle_faces.begin(), obstacle_faces.end(),
                             CS_TT_R::Parameters(55));

        CS_TT_R::Sample begin(0, 0, 0, 1);
        CS_TT_R::Sample end(1, 0, 0, 0);

        // find route
        cs.find_route(begin, end);
    }

    LOG4CXX_INFO(g_logger, "------------------------------ BB");
#endif

#if 0
    {
        CS::Loader_sphere_tree<Kernel> robot, obstacle;

        if (!robot.load_from_file("models/lamp600.sph") ||
            !obstacle.load_from_file("models/snakie.sph"))
        {
            LOG4CXX_ERROR(g_logger, "Failed to load sphere scene");
            return;
        }

        size_t maximum_number_of_levels = std::min(robot.number_of_levels(), obstacle.number_of_levels());
        size_t level = maximum_number_of_levels - 1;

        level = 2;

        CS_BB_C cs;
        cs.create_from_scene(robot.level_begin(level), robot.level_end(level),
                             obstacle.level_begin(level), obstacle.level_end(level),
                             CS_BB_C::Parameters(param));

        CS_BB_C::Sample begin(0, 0, 0, 1);
        CS_BB_C::Sample end(1, 0, 0, 0);

        // find route
        cs.find_route(begin, end);
    }
#endif

#if 0
    {
        CS::Loader_sphere_tree<Kernel> robot, obstacle;

        if (!robot.load_from_file("models/lamp600.sph") ||
            !obstacle.load_from_file("models/snakie.sph"))
        {
            LOG4CXX_ERROR(g_logger, "Failed to load sphere scene");
            return;
        }

        size_t maximum_number_of_levels = std::min(robot.number_of_levels(), obstacle.number_of_levels());
        size_t level = maximum_number_of_levels - 1;

        level = 2;

        CS_BB_R cs;

        cs.create_from_scene(robot.level_begin(level), robot.level_end(level),
                             obstacle.level_begin(level), obstacle.level_end(level),
                             CS_BB_R::Parameters(55));

        CS_BB_R::Sample begin(0, 0, 0, 1);
        CS_BB_R::Sample end(1, 0, 0, 0);

        // find route
        cs.find_route(begin, end);
    }
#endif

#if 0
    {
        Kernel::Ball_3 robot(Kernel::Vector_3(200, 0, 0), 1);
        Kernel::Ball_3 obstacle(Kernel::Vector_3(-200, 0, 0), 200);

        CS_BB_C cs;
        cs.create_from_scene(&robot, 1 + &robot,
                             &obstacle, 1 + &obstacle,
                             CS_BB_C::Parameters(param));

        CS_BB_C::Sample begin(0, 0, 0, 1);
        CS_BB_C::Sample end(1, 0, 0, 0);

        // find route
        cs.find_route(begin, end);
    }
#endif

#if 0
    {
        Kernel::Ball_3 robot(Kernel::Vector_3(2, 0, 0), 1);
        Kernel::Ball_3 obstacle(Kernel::Vector_3(-1, 0, 0), 1);

        CS_BB_R cs;
        cs.create_from_scene(&robot, 1 + &robot,
                             &obstacle, 1 + &obstacle,
                             CS_BB_R::Parameters(param));

        CS_BB_R::Sample begin(0, 0, 0, 1);
        CS_BB_R::Sample end(1, 0, 0, 0);

        // find route
        cs.find_route(begin, end);
    }
#endif

#if 1
    {
        Kernel::Ball_3 robot(Kernel::Vector_3(2, 0, 0), 1);
        Kernel::Ball_3 obstacle(Kernel::Vector_3(-1, 0, 0), 1);

        typedef CS_BB_E CS;

        CS cs;
        cs.create_from_scene(&robot, 1 + &robot,
                             &obstacle, 1 + &obstacle,
                             CS::Parameters());

        CS::Sample begin(0, 0, 0, 1);
        CS::Sample end(1, 0, 0, 0);

        // find route
        //cs.find_route(begin, end);
    }
#endif

#if 0
    {
        Kernel::Ball_3 robot(Kernel::Vector_3(2, 0, 0), 1);

        Kernel::Ball_3 obstacles[3] = { Kernel::Ball_3(Kernel::Vector_3(-1,  1,  0), 1),
                                        Kernel::Ball_3(Kernel::Vector_3(-1,  0,  0), 1),
                                        Kernel::Ball_3(Kernel::Vector_3(-1, -1,  0), 1) };

        typedef CS_BB_R CS;

        CS cs;
        cs.create_from_scene(&robot, 1 + &robot,
                             obstacles, 3 + obstacles,
                             CS::Parameters(param));

        CS::Sample begin(0, 0, 0, 1);
        CS::Sample end(1, 0, 0, 0);

        // find route
        cs.find_route(begin, end);
    }
#endif

#if 0
    {
        Kernel::Ball_3 robot(Kernel::Vector_3(2, 0, 0), 1);

        Kernel::Ball_3 obstacles[4] = { Kernel::Ball_3(Kernel::Vector_3(-1,  0,  0), 1),
                                        Kernel::Ball_3(Kernel::Vector_3(-1,  1,  0), 1),
                                        Kernel::Ball_3(Kernel::Vector_3(-2,  0,  0), 1),
                                        Kernel::Ball_3(Kernel::Vector_3(-1, -1,  0), 1) };

        CS_BB_E cs;
        cs.create_from_scene(&robot, 1 + &robot,
                             obstacles, 4 + obstacles,
                             CS_BB_E::Parameters());

        CS_BB_E::Sample begin(0, 0, 0, 1);
        CS_BB_E::Sample end(1, 0, 0, 0);

        // find route
        cs.find_route(begin, end);
    }
#endif
}
