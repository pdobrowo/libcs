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
#include "Predicate_bb_3_list_generator.h"

namespace CS
{
template<class Kernel_>
template<typename RobotInputIterator, typename ObstacleInputIterator, typename OutputIterator>
void Predicate_list_generator<Kernel_, Predicate_bb_3<Kernel_> >::create_predicate_list(
        RobotInputIterator robot_begin, RobotInputIterator robot_end,
        ObstacleInputIterator obstacle_begin, ObstacleInputIterator obstacle_end,
        OutputIterator predicates_iterator)
{
    //typedef typename Kernel_::Point_3   Point_3;
    typedef typename Kernel_::Vector_3  Vector_3;
    typedef typename Kernel_::RT        RT;

    // create BB predicate list
    for (ObstacleInputIterator i = obstacle_begin; i != obstacle_end; ++i)
    {
        for (RobotInputIterator j = robot_begin; j != robot_end; ++j)
        {
            // obstacle
            Vector_3 center_b = i->center();
            RT radius_b = i->radius();

            // rotating
            Vector_3 center_a = j->center();
            RT radius_a = j->radius();

            RT sqr_distance_b = center_b.squared_length();
            RT sqr_distance_a = center_a.squared_length();

            // discard all predicates that are empty for sure
            RT ra_p_rb = radius_a + radius_b;
            RT left = sqr_distance_a + sqr_distance_b - ra_p_rb * ra_p_rb;

            if (left <= RT(0) ||
                left * left <= RT(4) * sqr_distance_a * sqr_distance_b)
            {
                // include this predicate
                *(predicates_iterator++) = Predicate_bb_3<Kernel_>(center_b, radius_b, center_a, radius_a);
            }
        }
    }
}
} // namespace CS
