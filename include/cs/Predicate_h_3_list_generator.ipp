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
#include "Predicate_h_3_list_generator.h"

namespace CS
{
template<class Kernel_>
template<typename Polyhedron_, typename OutputIterator_>
void Predicate_list_generator<Kernel_, Predicate_h_3<Kernel_> >::convex_decomposition(const Polyhedron_ &polyhedron, OutputIterator_ convex_polyhedron_iterator)
{
    typedef Polyhedron_ Polyhedron;
    typedef CGAL::Nef_polyhedron_3<Kernel_/*, CGAL::SNC_indexed_items*/> Nef_polyhedron_3;
    typedef typename Nef_polyhedron_3::Volume_const_iterator Volume_const_iterator;

    // need to make a copy of polyhedron
    Polyhedron polyhedron_copy((polyhedron));

    // convert source to nef
    Nef_polyhedron_3 nef_polyhedron((polyhedron_copy));

    // decompose nef
    CGAL::convex_decomposition_3(nef_polyhedron);

    // the first volume is the outer volume, which is ignored in the decomposition
    Volume_const_iterator ci = ++nef_polyhedron.volumes_begin();

    for (; ci != nef_polyhedron.volumes_end(); ++ci)
    {
        if (ci->mark())
        {
            // convert nef to polyhedron
            Polyhedron convex_polyhedron;
            nef_polyhedron.convert_inner_shell_to_polyhedron(ci->shells_begin(), convex_polyhedron);
            *(convex_polyhedron_iterator++) = convex_polyhedron;
        }
    }
}

template<class Kernel_>
template<typename RobotInputIterator_, typename ObstacleInputIterator_, typename OutputIterator_>
void Predicate_list_generator<Kernel_, Predicate_h_3<Kernel_> >::create_predicate_list(
        RobotInputIterator_ robot_begin, RobotInputIterator_ robot_end,
        ObstacleInputIterator_ obstacle_begin, ObstacleInputIterator_ obstacle_end,
        OutputIterator_ predicates_iterator)
{
    typedef typename Kernel_::Polyhedron_3              Polyhedron_3;
    typedef typename Kernel_::Vector_3                  Vector_3;

    typedef typename Polyhedron_3::Plane_const_iterator Plane_const_iterator;
    typedef typename Polyhedron_3::Plane_3              Plane_3;

    typedef typename Polyhedron_3::Point_const_iterator Point_const_iterator;
    typedef typename Polyhedron_3::Point_3              Point_3;

    // decompose obstacles
    std::list<Polyhedron_3> obstacle_parts;

    for (ObstacleInputIterator_ obstacle_iterator = obstacle_begin; obstacle_iterator != obstacle_end; ++obstacle_iterator)
    {
        const Polyhedron_3 &obstacle = *obstacle_iterator;
        convex_decomposition(obstacle, std::back_inserter(obstacle_parts));
    }

    // decompose robot
    std::list<Polyhedron_3> robot_parts;

    for (RobotInputIterator_ robot_iterator = robot_begin; robot_iterator != robot_end; ++robot_iterator)
    {
        const Polyhedron_3 &robot = *robot_iterator;
        convex_decomposition(robot, std::back_inserter(robot_parts));
    }

    // create H predicate list
    for (typename std::list<Polyhedron_3>::const_iterator obstacle_part_iterator = obstacle_parts.begin(); obstacle_part_iterator != obstacle_parts.end(); ++obstacle_part_iterator)
    {
        const Polyhedron_3 &obstacle_part = *obstacle_part_iterator;

        for (typename std::list<Polyhedron_3>::const_iterator robot_part_iterator = obstacle_parts.begin(); robot_part_iterator != obstacle_parts.end(); ++robot_part_iterator)
        {
            const Polyhedron_3 &robot_part = *robot_part_iterator;

            // stationary obstacle plane and rotational robot vertex
            for (Plane_const_iterator obstacle_part_plane_iterator = obstacle_part.planes_begin(); obstacle_part_plane_iterator != obstacle_part.planes_end(); ++obstacle_part_plane_iterator)
            {
                for (Point_const_iterator robot_part_point_iterator = robot_part.points_begin(); robot_part_point_iterator != robot_part.points_end(); ++robot_part_point_iterator)
                {
                    const Plane_3 &obstacle_plane = *obstacle_part_plane_iterator;
                    const Point_3 &robot_point = *robot_part_point_iterator;

                    // include this predicate
                    *(predicates_iterator++) = Predicate_h_3<Kernel_>(Vector_3(robot_point.x(), robot_point.y(), robot_point.z()), obstacle_plane);
                }
            }

            // stationary obstacle vertex and rotational robot plane
            for (Plane_const_iterator robot_part_plane_iterator = robot_part.planes_begin(); robot_part_plane_iterator != robot_part.planes_end(); ++robot_part_plane_iterator)
            {
                for (Point_const_iterator obstacle_part_point_iterator = obstacle_part.points_begin(); obstacle_part_point_iterator != obstacle_part.points_end(); ++obstacle_part_point_iterator)
                {
                    const Plane_3 &robot_plane = *robot_part_plane_iterator;
                    const Point_3 &obstacle_point = *obstacle_part_point_iterator;

                    // include this predicate
                    *(predicates_iterator++) = Predicate_h_3<Kernel_>(Vector_3(-obstacle_point.x(), -obstacle_point.y(), -obstacle_point.z()), robot_plane);
                }
            }
        }
    }
}
} // namespace CS
