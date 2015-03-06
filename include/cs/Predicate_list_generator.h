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
#ifndef LIBCS_PREDICATE_LIST_GENERATOR_H
#define LIBCS_PREDICATE_LIST_GENERATOR_H

namespace CS
{
template<class Kernel_, class Predicate_>
struct Predicate_list_generator
{
    template<typename RobotInputIterator_, typename ObstacleInputIterator_, typename OutputIterator_>
    void create_predicate_list(
        RobotInputIterator_ robot_begin, RobotInputIterator_ robot_end,
        ObstacleInputIterator_ obstacle_begin, ObstacleInputIterator_ obstacle_end,
        OutputIterator_ predicates_iterator);

private:
    template<typename Polyhedron_, typename OutputIterator_>
    void convex_decomposition(const Polyhedron_ &polyhedron, OutputIterator_ convex_polyhedron_iterator);
};
} // namespace CS

#endif // LIBCS_PREDICATE_LIST_GENERATOR_H
