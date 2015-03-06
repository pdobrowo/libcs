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
#ifndef LIBCS_PREDICATE_BB_3_LIST_GENERATOR_H
#define LIBCS_PREDICATE_BB_3_LIST_GENERATOR_H

#include "Predicate_list_generator.h"
#include "Predicate_bb_3.h"

namespace CS
{
template<class Kernel_>
struct Predicate_list_generator<Kernel_, Predicate_bb_3<Kernel_> >
{
    template<typename RobotInputIterator_, typename ObstacleInputIterator_, typename OutputIterator_>
    void create_predicate_list(
        RobotInputIterator_ robot_begin, RobotInputIterator_ robot_end,
        ObstacleInputIterator_ obstacle_begin, ObstacleInputIterator_ obstacle_end,
        OutputIterator_ predicates_iterator);
};
} // namespace CS

#include "Predicate_bb_3_list_generator.ipp"

#endif // LIBCS_PREDICATE_BB_3_LIST_GENERATOR_H
