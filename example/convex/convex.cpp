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
#include "convex.h"
#include <cs/Spin_inexact_kernel_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Gmpz.h>
#include <fstream>

// kernel
typedef CGAL::Gmpz RT;
typedef CGAL::Cartesian<RT> Kernel_base;
typedef CS::Spin_inexact_kernel_3<Kernel_base> Kernel;

// scene
typedef Kernel::Predicate_h_3 Predicate_h_3;
typedef Kernel::Spin_null_configuration_space_3<Predicate_h_3>::Type Configuration_space_3;
typedef Configuration_space_3::Parameters Parameters;

// object
typedef Kernel::Point_3 Point_3;
typedef Kernel::Polyhedron_3 Polyhedron_3;

void convex()
{
    // robot
    Polyhedron_3 robot;

    /*
    robot.make_tetrahedron(Point_3(-1, -1, -1),
                           Point_3(+1, -1, -1),
                           Point_3(+0, +1, -1),
                           Point_3(+0, +0, +2));
    */

    std::ifstream robot_file("models/cube.off");

    if (!robot_file.is_open())
        exit(-1);

    robot_file >> robot;

    // obstacle
    Polyhedron_3 obstacle;

    /*
    obstacle.make_tetrahedron(Point_3(+1, -1, -1),
                              Point_3(+3, -1, -1),
                              Point_3(+2, +1, -1),
                              Point_3(+2, +0, +1));
    */

    std::ifstream obstacle_file("models/cube.off");

    if (!obstacle_file.is_open())
        exit(-1);

    obstacle_file >> obstacle;

    // configuration space
    Configuration_space_3 cs;

    cs.create_from_scene(&robot, &robot + 1,
                         &obstacle, &obstacle + 1,
                         Parameters());
}
