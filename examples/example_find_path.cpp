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
#include "example_find_path.h"

#include <cs/Spin_kernel_3.h>
#include <CGAL/Cartesian.h>

// kernel settings
typedef double RT;
typedef CGAL::Cartesian<RT> Kernel_base;
typedef CS::Spin_inexact_kernel_3<Kernel_base> Kernel;

// scene type
typedef Kernel::Predicate_tt_3 Predicate_tt_3;
typedef Kernel::Spin_raster_configuration_space_3<Predicate_tt_3>::Type Configuration_space_3;

// objects
typedef Kernel::Triangle_3 Triangle_3;
typedef Kernel::Point_3 Point_3;

// representation
typedef CS::Spin_raster_graph_3<Kernel, Predicate_tt_3> Spin_raster_graph_3;

// path searching
typedef Configuration_space_3::Sample Sample;
typedef Configuration_space_3::Route Route;

void example_find_path()
{
    // rotating object
    Triangle_3 rotating(Point_3(0, 0, 0), Point_3(0, 1, 0), Point_3(1, 0, 0));

    // obstacle
    Triangle_3 obstacle(Point_3(0, 0, 0.5), Point_3(0, 1, 0.5), Point_3(1, 0, 0.5));

    // configuration space
    Configuration_space_3 cs;
    cs.create_from_scene(&rotating, &rotating + 1,
                         &obstacle, &obstacle + 1,
                         Configuration_space_3::Parameters(64));

    // find a path
    Sample begin(0, 0, 0, 1); // s0 = 1
    Sample end(1, 0, 0, 0); // s1 = e12

    Route route = cs.find_route(begin, end);

    // check if path was found
    if (!route.is_valid())
        return;

    // read all path nodes
    std::vector<Sample> path_nodes = route.nodes();
    std::size_t number_of_path_nodes = path_nodes.size();

    // calculate configurations along the path
    for (int i = 0; i <= 100; ++i)
    {
        double t = 0.01 * i;
        Sample sample = route.evaluate(t);

        // sample is the rotation of the object at the moment t, where t is in [0; 1]
    }
}
