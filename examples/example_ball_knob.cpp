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
#include "example_ball_knob.h"

#include <cs/Spin_kernel_3.h>
#include <cs/Bigint.h>
#include <CGAL/Cartesian.h>

struct Floating_point_inexact_kernel_raster_graph
{
    // kernel settings
    typedef double RT;
    typedef CGAL::Cartesian<RT> Kernel_base;
    typedef CS::Spin_inexact_kernel_3<Kernel_base> Kernel;

    // scene type
    typedef Kernel::Predicate_bb_3 Predicate_bb_3;
    typedef Kernel::Spin_raster_configuration_space_3<Predicate_bb_3>::Type Configuration_space_3;

    static Configuration_space_3::Parameters parameters()
    {
        return Configuration_space_3::Parameters(128);
    }

    // objects
    typedef Kernel::Ball_3 Ball_3;
    typedef Kernel::Point_3 Point_3;
};

struct Floating_point_inexact_kernel_cell_graph
{
    // kernel settings
    typedef double RT;
    typedef CGAL::Cartesian<RT> Kernel_base;
    typedef CS::Spin_inexact_kernel_3<Kernel_base> Kernel;

    // scene type
    typedef Kernel::Predicate_bb_3 Predicate_bb_3;
    typedef Kernel::Spin_cell_configuration_space_3<Predicate_bb_3>::Type Configuration_space_3;

    static Configuration_space_3::Parameters parameters()
    {
        return Configuration_space_3::Parameters(1000000);
    }

    // objects
    typedef Kernel::Ball_3 Ball_3;
    typedef Kernel::Point_3 Point_3;
};

struct Bigint_exact_kernel_exact_graph
{
    // kernel settings
    typedef CGAL::Bigint RT;
    typedef CGAL::Cartesian<RT> Kernel_base;
    typedef CS::Spin_kernel_3<Kernel_base> Kernel;

    // scene type
    typedef Kernel::Predicate_bb_3 Predicate_bb_3;
    typedef Kernel::Spin_exact_configuration_space_3<Predicate_bb_3>::Type Configuration_space_3;

    static Configuration_space_3::Parameters parameters()
    {
        return Configuration_space_3::Parameters();
    }

    // objects
    typedef Kernel::Ball_3 Ball_3;
    typedef Kernel::Point_3 Point_3;
};

template<typename Traits>
void example_ball_knob_generic()
{
    // bring types here
    typedef typename Traits::Point_3 Point_3;
    typedef typename Traits::Ball_3 Ball_3;

    typedef typename Traits::Configuration_space_3 Configuration_space_3;

    // ball knob
    const Ball_3 knob_balls[4] = {
        Ball_3(0, 0, 0, 1),
        Ball_3(0, 0, 2, 1),
        Ball_3(0, 0, 4, 1),
        Ball_3(0, 2, 4, 1)
    };

    // ball body
    const Ball_3 body_balls[8] = {
        Ball_3(4, 0, 0, 1),
        Ball_3(-4, 0, 0, 1),
        Ball_3(0, 4, 0, 1),
        Ball_3(0, -4, 0, 1),
        Ball_3(3, 3, 0, 1),
        Ball_3(3, -3, 0, 1),
        Ball_3(-3, -3, 0, 1),
        Ball_3(-3, 3, 0, 1),
    };

    // configuration space
    Configuration_space_3 cs;

    cs.create_from_scene(knob_balls, knob_balls + 4,
                         body_balls, body_balls + 8,
                         Traits::parameters());
}

void example_ball_knob()
{
    //example_ball_knob_generic<Floating_point_inexact_kernel_raster_graph>();
    //example_ball_knob_generic<Floating_point_inexact_kernel_cell_graph>();
    example_ball_knob_generic<Bigint_exact_kernel_exact_graph>();
}
