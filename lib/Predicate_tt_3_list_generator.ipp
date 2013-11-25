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
#include "Predicate_tt_3_list_generator.h"

namespace CS
{
template<class Kernel_>
template<typename RobotInputIterator, typename ObstacleInputIterator, typename OutputIterator>
void Predicate_list_generator<Kernel_, Predicate_tt_3<Kernel_> >::create_predicate_list(
        RobotInputIterator robot_begin, RobotInputIterator robot_end,
        ObstacleInputIterator obstacle_begin, ObstacleInputIterator obstacle_end,
        OutputIterator predicates_iterator)
{
    typedef typename Kernel_::Point_3   Point_3;
    typedef typename Kernel_::Vector_3  Vector_3;
    typedef typename Kernel_::RT        RT;

    // create TT predicate list
    Point_3 ZERO(RT(0), RT(0), RT(0));

    for (ObstacleInputIterator i = obstacle_begin; i != obstacle_end; ++i)
    {
        for (RobotInputIterator j = robot_begin; j != robot_end; ++j)
        {
            // obstacle
            Vector_3 vk = i->vertex(0) - ZERO;
            Vector_3 vl = i->vertex(1) - ZERO;
            Vector_3 vm = i->vertex(2) - ZERO;

            // rotating
            Vector_3 va = j->vertex(0) - ZERO;
            Vector_3 vb = j->vertex(1) - ZERO;
            Vector_3 vc = j->vertex(2) - ZERO;

            // discard all predicates that are empty for sure
            RT robot_face_vertex_squared_distance[3] = { va.squared_length(), vb.squared_length(), vc.squared_length() };

            // get maximum rotating distance

            // FIXME: handle minimum distance too
            RT max_robot_face_vertex_squared_distance =
                    std::max(robot_face_vertex_squared_distance[0], std::max(robot_face_vertex_squared_distance[1], robot_face_vertex_squared_distance[2]));

            // discard empty collisions
            if (!sphere_triangle_intersection_3(vk, vl, vm, max_robot_face_vertex_squared_distance, Vector_3(0, 0, 0)))
                continue;

            // accumulate predicates
            *(predicates_iterator++) = Predicate_tt_3<Kernel_>(vk, vl, vm, va, vb, vc);
        }
    }
}
} // namespace CS
