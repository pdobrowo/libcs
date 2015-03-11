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
#include "Spin_configuration_space_3.h"

namespace CS
{
template<class Kernel_, class Predicate_, class Representation_>
Spin_configuration_space_3<Kernel_, Predicate_, Representation_>::Spin_configuration_space_3()
{
}

template<class Kernel_, class Predicate_, class Representation_>
template<typename RobotInputIterator, typename ObstacleInputIterator>
void Spin_configuration_space_3<Kernel_, Predicate_, Representation_>::create_from_scene(
        RobotInputIterator robot_begin, RobotInputIterator robot_end,
        ObstacleInputIterator obstacle_begin, ObstacleInputIterator obstacle_end,
        const Parameters &parameters)
{
    size_t number_of_robot_primitives = std::distance(robot_begin, robot_end);
    size_t number_of_obstacle_primitives = std::distance(obstacle_begin, obstacle_end);

    CS_logger_info(MODULE, "Creating configuration space");
    CS_logger_info(MODULE, "Configuration: " << number_of_robot_primitives << " robot primitives; " << number_of_obstacle_primitives << " obstacle primitives");

    // timing
    unsigned long long timer_start, timer_end, timer_all_start, timer_all_end;

    // create TT predicate list
    CS_logger_debug(MODULE, "Creating predicate list");

    timer_all_start = get_tick_count();
    timer_start = get_tick_count();

    Generator generator;

    generator.create_predicate_list(robot_begin, robot_end,
                                    obstacle_begin, obstacle_end,
                                    std::back_inserter(m_predicates));

    // note: the following calculation is not valid for H-scenes created from polyhedrons
    size_t total_number_of_scene_predicates = number_of_obstacle_primitives * number_of_robot_primitives;

    timer_end = get_tick_count();

    CS_logger_info(MODULE, "Prepared " << m_predicates.size() << "/" << total_number_of_scene_predicates << " scene predicates, " << std::fixed << std::setprecision(2) << (float(100 * m_predicates.size()) / float(total_number_of_scene_predicates)) << "% [" << (timer_end - timer_start) << " ms]");

    // prepare general predicate list
    CS_logger_debug(MODULE, "Creating general predicate list");

    CS_logger_debug(MODULE, "Sub-predicate count is " << SUB_PREDICATE_COUNT);

    m_general_predicates.reserve(SUB_PREDICATE_COUNT * m_predicates.size());

    for (const Predicate &predicate: m_predicates)
        for (size_t j = 0; j < SUB_PREDICATE_COUNT; ++j)
            m_general_predicates.push_back(Predicate_g_3(predicate.sub_predicates()[j]));

    // create spin quadric list
    CS_logger_debug(MODULE, "Creating spin quadric list");

    m_spin_quadrics.reserve(m_general_predicates.size());

    for (const Predicate_g_3 &general_predicate: m_general_predicates)
        m_spin_quadrics.push_back(Spin_quadric_3(general_predicate));

    CS_logger_info(MODULE, "Prepared " << m_spin_quadrics.size() << " spin-quadrics [" << (timer_end - timer_start) << " ms]");

    // create representation
    CS_logger_debug(MODULE, "Creating representation");

    m_representation.reset(new Representation(m_predicates, m_spin_quadrics, parameters));

    timer_all_end = get_tick_count();
    CS_logger_info(MODULE, "Configuration space created [" << (timer_all_end - timer_all_start) << " ms]");
}

template<class Kernel_, class Predicate_, class Representation_>
typename Spin_configuration_space_3<Kernel_, Predicate_, Representation_>::Route Spin_configuration_space_3<Kernel_, Predicate_, Representation_>::find_route(const Sample &begin, const Sample &end)
{
    CS_logger_info(MODULE, "Find route: " << begin << " -> " << end);

    // select upper hemisphere of Spin(3) as the representation
    // relocate begin and end locations if needed
    Sample normalized_begin = begin.s0() >= 0 ? begin : -begin;
    Sample normalized_end = end.s0() >= 0 ? end : -end;

    CS_logger_info(MODULE, "Normalized find route: " << begin << " -> " << end);

    // use representation
    return rep().find_route(normalized_begin, normalized_end);
}

template<class Kernel_, class Predicate_, class Representation_>
typename Spin_configuration_space_3<Kernel_, Predicate_, Representation_>::Representation &Spin_configuration_space_3<Kernel_, Predicate_, Representation_>::rep()
{
    assert(!!m_representation);
    return *m_representation;
}

template<class Kernel_, class Predicate_, class Representation_>
const typename Spin_configuration_space_3<Kernel_, Predicate_, Representation_>::Representation &Spin_configuration_space_3<Kernel_, Predicate_, Representation_>::rep() const
{
    assert(!!m_representation);
    return *m_representation;
}

template<class Kernel_, class Predicate_, class Representation_>
typename Spin_configuration_space_3<Kernel_, Predicate_, Representation_>::Spin_quadric_const_iterator Spin_configuration_space_3<Kernel_, Predicate_, Representation_>::spin_quadrics_begin() const
{
    return m_spin_quadrics.begin();
}

template<class Kernel_, class Predicate_, class Representation_>
typename Spin_configuration_space_3<Kernel_, Predicate_, Representation_>::Spin_quadric_const_iterator Spin_configuration_space_3<Kernel_, Predicate_, Representation_>::spin_quadrics_end() const
{
    return m_spin_quadrics.end();
}

template<class Kernel_, class Predicate_, class Representation_>
typename Spin_configuration_space_3<Kernel_, Predicate_, Representation_>::Spin_quadric_size_type Spin_configuration_space_3<Kernel_, Predicate_, Representation_>::size_of_spin_quadrics() const
{
    return m_spin_quadrics.size();
}

template<class Kernel_, class Predicate_, class Representation_>
typename Spin_configuration_space_3<Kernel_, Predicate_, Representation_>::Predicate_const_iterator Spin_configuration_space_3<Kernel_, Predicate_, Representation_>::predicates_begin() const
{
    return m_predicates.begin();
}

template<class Kernel_, class Predicate_, class Representation_>
typename Spin_configuration_space_3<Kernel_, Predicate_, Representation_>::Predicate_const_iterator Spin_configuration_space_3<Kernel_, Predicate_, Representation_>::predicates_end() const
{
    return m_predicates.end();
}

template<class Kernel_, class Predicate_, class Representation_>
typename Spin_configuration_space_3<Kernel_, Predicate_, Representation_>::Predicate_size_type Spin_configuration_space_3<Kernel_, Predicate_, Representation_>::size_of_predicates() const
{
    return m_predicates.size();
}

template<class Kernel_, class Predicate_, class Representation_>
typename Spin_configuration_space_3<Kernel_, Predicate_, Representation_>::General_predicate_const_iterator Spin_configuration_space_3<Kernel_, Predicate_, Representation_>::general_predicates_begin() const
{
    return m_general_predicates.begin();
}

template<class Kernel_, class Predicate_, class Representation_>
typename Spin_configuration_space_3<Kernel_, Predicate_, Representation_>::General_predicate_const_iterator Spin_configuration_space_3<Kernel_, Predicate_, Representation_>::general_predicates_end() const
{
    return m_general_predicates.end();
}

template<class Kernel_, class Predicate_, class Representation_>
typename Spin_configuration_space_3<Kernel_, Predicate_, Representation_>::General_predicate_size_type Spin_configuration_space_3<Kernel_, Predicate_, Representation_>::size_of_general_predicates() const
{
    return m_general_predicates.size();
}
} // namespace CS
