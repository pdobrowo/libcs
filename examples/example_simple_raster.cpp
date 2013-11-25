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
#include "example_simple_raster.h"

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

void example_simple_raster()
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

    // get access to configuration space representation
    const Spin_raster_graph_3 &raster_graph = cs.rep();

    // access individual voxels in raster
    const Spin_raster_graph_3::Voxel &voxel = raster_graph.voxel(0, 0, 0);
}
