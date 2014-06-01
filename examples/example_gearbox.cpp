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
#include "example_gearbox.h"

#define EXACT 0

#if EXACT
#include <cs/Spin_kernel_3.h>
#include <cs/Bigint.h>
#else // EXACT
#include <cs/Spin_inexact_kernel_3.h>
#include <cs/Coordinate.h>
#endif // EXACT
#include <CGAL/Cartesian.h>

#include <iostream>

struct Floating_point_inexact_kernel_raster_graph
{
    // kernel settings
    typedef double RT;
    typedef CGAL::Cartesian<RT> Kernel_base;
    typedef CS::Spin_inexact_kernel_3<Kernel_base> Kernel;

    // scene type
    typedef Kernel::Predicate_tt_3 Predicate_tt_3;
    typedef Kernel::Spin_raster_configuration_space_3<Predicate_tt_3>::Type Configuration_space_3;

    static Configuration_space_3::Parameters parameters()
    {
        return Configuration_space_3::Parameters(128);
    }

    // objects
    typedef Kernel::Triangle_3 Triangle_3;
    typedef Kernel::Point_3 Point_3;
};

struct Floating_point_inexact_kernel_cell_graph
{
    // kernel settings
    typedef double RT;
    typedef CGAL::Cartesian<RT> Kernel_base;
    typedef CS::Spin_inexact_kernel_3<Kernel_base> Kernel;

    // scene type
    typedef Kernel::Predicate_tt_3 Predicate_tt_3;
    typedef Kernel::Spin_cell_configuration_space_3<Predicate_tt_3>::Type Configuration_space_3;

    static Configuration_space_3::Parameters parameters()
    {
        return Configuration_space_3::Parameters(1000000);
    }

    // objects
    typedef Kernel::Triangle_3 Triangle_3;
    typedef Kernel::Point_3 Point_3;
};

#if EXACT
struct Bigint_exact_kernel_exact_graph
{
    // kernel settings
    typedef CGAL::Bigint RT;
    typedef CGAL::Cartesian<RT> Kernel_base;
    typedef CS::Spin_kernel_3<Kernel_base> Kernel;

    // scene type
    typedef Kernel::Predicate_tt_3 Predicate_tt_3;
    typedef Kernel::Spin_exact_configuration_space_3<Predicate_tt_3>::Type Configuration_space_3;

    static Configuration_space_3::Parameters parameters()
    {
        return Configuration_space_3::Parameters();
    }

    // objects
    typedef Kernel::Triangle_3 Triangle_3;
    typedef Kernel::Point_3 Point_3;
};
#endif // EXACT

template<typename Traits>
void example_gearbox_generic()
{
    // bring types here
    typedef typename Traits::Point_3 Point_3;
    typedef typename Traits::Triangle_3 Triangle_3;

    typedef typename Traits::Configuration_space_3 Configuration_space_3;

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
                         Traits::parameters());

#if 1

    typedef typename Configuration_space_3::Representation Representation;
    const Representation &rep = cs.rep();

    typedef typename Representation::Cell_const_iterator Cell_const_iterator;
    typedef typename Representation::Cell Cell;

    size_t index = 0;

    for (Cell_const_iterator cell_iterator = rep.cells_begin(); cell_iterator != rep.cells_end(); ++cell_iterator)
    {
        const Cell &cell = *cell_iterator;
        CS::Coordinate coordinate = cell.coordinate();

        std::cout << index << ": [";

        for (CS::Coordinate::const_iterator bit_iterator = coordinate.begin(); bit_iterator != coordinate.end(); ++bit_iterator)
        {
            if (bit_iterator == coordinate.begin())
                std::cout << (*bit_iterator ? "1" : "0");
            else
                std::cout << (*bit_iterator ? ", 1" : ", 0");
        }

        std::cout << "]" << std::endl;

        ++index;
    }

#endif
}

void example_gearbox()
{
    //example_gearbox_generic<Floating_point_inexact_kernel_raster_graph>();
    example_gearbox_generic<Floating_point_inexact_kernel_cell_graph>();
    //example_gearbox_generic<Bigint_exact_kernel_exact_graph>();
}
