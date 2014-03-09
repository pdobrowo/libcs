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
#include "Spin_configuration_space_3.h"

//#define DISPLAY_SPIN_SURFACE_ASSOCIATED_ELLIPSOID

#ifdef DISPLAY_SPIN_SURFACE_ASSOCIATED_ELLIPSOID
    #include <iostream>
    #include <eigen3/Eigen/Dense>
#endif // DISPLAY_SPIN_SURFACE_ASSOCIATED_ELLIPSOID

namespace CS
{
template<class K, class P, class R>
Spin_configuration_space_3<K, P, R>::Spin_configuration_space_3()
    : m_logger(log4cxx::Logger::getLogger("CS.Spin_configuration_space_3"))
{
}

template<class K, class P, class R>
template<typename RobotInputIterator, typename ObstacleInputIterator>
void Spin_configuration_space_3<K, P, R>::create_from_scene(
        RobotInputIterator robot_begin, RobotInputIterator robot_end,
        ObstacleInputIterator obstacle_begin, ObstacleInputIterator obstacle_end,
        const Parameters &parameters)
{
    size_t number_of_robot_faces = std::distance(robot_begin, robot_end);
    size_t number_of_obstacle_faces = std::distance(obstacle_begin, obstacle_end);

    LOG4CXX_INFO(m_logger, "Creating configuration space");
    LOG4CXX_INFO(m_logger, "Configuration: " << number_of_robot_faces << " robot primitives; " << number_of_obstacle_faces << " obstacle primitives");

    // timing
    unsigned long long timer_start, timer_end, timer_all_start, timer_all_end;

    // create TT predicate list
    LOG4CXX_DEBUG(m_logger, "Creating predicate list");

    timer_all_start = get_tick_count();
    timer_start = get_tick_count();

    Generator generator;

    generator.create_predicate_list(robot_begin, robot_end,
                                    obstacle_begin, obstacle_end,
                                    std::back_inserter(m_predicates));

    size_t total_number_of_predicates = number_of_obstacle_faces * number_of_robot_faces;

    timer_end = get_tick_count();

    LOG4CXX_INFO(m_logger, "Prepared " << m_predicates.size() << "/" << total_number_of_predicates << " predicates, " << std::fixed << std::setprecision(2) << (float(100 * m_predicates.size()) / float(total_number_of_predicates)) << "% [" << (timer_end - timer_start) << " ms]");

    // prepare generic predicate list
    LOG4CXX_DEBUG(m_logger, "Creating generic predicate list");

    LOG4CXX_DEBUG(m_logger, "Sub-predicate count is " << SUB_PREDICATE_COUNT);

    m_generic_predicates.reserve(SUB_PREDICATE_COUNT * m_predicates.size());

    BOOST_FOREACH(const Predicate &predicate, m_predicates)
        for (size_t j = 0; j < SUB_PREDICATE_COUNT; ++j)
            m_generic_predicates.push_back(Predicate_g_3(predicate.sub_predicates()[j]));

    // create spin quadric list
    LOG4CXX_DEBUG(m_logger, "Creating spin quadric list");

    m_spin_quadrics.reserve(m_generic_predicates.size());

    BOOST_FOREACH(const Predicate_g_3 &generic_predicate, m_generic_predicates)
        m_spin_quadrics.push_back(Spin_quadric_3(generic_predicate));

#ifdef DISPLAY_SPIN_SURFACE_ASSOCIATED_ELLIPSOID

    // print spin-surface associated ellipsoids
    for (typename std::vector<Spin_quadric_3>::const_iterator it = m_quadrics.begin(); it != m_quadrics.end(); ++it)
    {
        typename K::Matrix ellipsoid_matrix = it->ellipsoid_matrix();

        Eigen::Matrix4d eigen_matrix;

        for (int row = 0; row < 4; ++row)
            for (int column = 0; column < 4; ++column)
                eigen_matrix(row, column) = ellipsoid_matrix(row, column);

        std::cout << "spin-quadric associated ellipsoid:" << std::endl << eigen_matrix << std::endl;

        Eigen::EigenSolver<Eigen::Matrix4d> eigensolver(eigen_matrix);

        if (eigensolver.info() == Eigen::Success)
            std::cout << "eigenvalues:" << std::endl << eigensolver.eigenvalues() << std::endl;
        else
            std::cout << "failed to calculate eigenvalues!" << std::endl;
    }

#endif // DISPLAY_SPIN_SURFACE_ASSOCIATED_ELLIPSOID

    LOG4CXX_INFO(m_logger, "Prepared " << m_spin_quadrics.size() << " spin-quadrics [" << (timer_end - timer_start) << " ms]");

    // create representation
    LOG4CXX_DEBUG(m_logger, "Creating representation");

    m_representation.reset(new Representation(m_predicates, m_spin_quadrics, parameters));

    timer_all_end = get_tick_count();
    LOG4CXX_INFO(m_logger, "Configuration space created [" << (timer_all_end - timer_all_start) << " ms]");
}

template<class K, class P, class R>
typename Spin_configuration_space_3<K, P, R>::Route Spin_configuration_space_3<K, P, R>::find_route(const Sample &begin, const Sample &end)
{
    LOG4CXX_INFO(m_logger, "Find route: " << begin << " -> " << end);

    // select upper hemisphere of Spin(3) as the representation
    // relocate begin and end locations if needed
    Sample normalized_begin = begin.s0() >= 0 ? begin : -begin;
    Sample normalized_end = end.s0() >= 0 ? end : -end;

    LOG4CXX_INFO(m_logger, "Normalized find route: " << begin << " -> " << end);

    // use representation
    return rep().find_route(normalized_begin, normalized_end);
}

template<class K, class P, class R>
typename Spin_configuration_space_3<K, P, R>::Representation &Spin_configuration_space_3<K, P, R>::rep()
{
    assert(!!m_representation);
    return *m_representation;
}

template<class K, class P, class R>
const typename Spin_configuration_space_3<K, P, R>::Representation &Spin_configuration_space_3<K, P, R>::rep() const
{
    assert(!!m_representation);
    return *m_representation;
}

template<class K, class P, class R>
typename Spin_configuration_space_3<K, P, R>::Spin_quadric_const_iterator Spin_configuration_space_3<K, P, R>::spin_quadrics_begin() const
{
    return m_spin_quadrics.begin();
}

template<class K, class P, class R>
typename Spin_configuration_space_3<K, P, R>::Spin_quadric_const_iterator Spin_configuration_space_3<K, P, R>::spin_quadrics_end() const
{
    return m_spin_quadrics.end();
}

template<class K, class P, class R>
typename Spin_configuration_space_3<K, P, R>::Spin_quadric_size_type Spin_configuration_space_3<K, P, R>::size_of_spin_quadrics() const
{
    return m_spin_quadrics.size();
}

template<class K, class P, class R>
typename Spin_configuration_space_3<K, P, R>::Predicate_const_iterator Spin_configuration_space_3<K, P, R>::predicates_begin() const
{
    return m_predicates.begin();
}

template<class K, class P, class R>
typename Spin_configuration_space_3<K, P, R>::Predicate_const_iterator Spin_configuration_space_3<K, P, R>::predicates_end() const
{
    return m_predicates.end();
}

template<class K, class P, class R>
typename Spin_configuration_space_3<K, P, R>::Predicate_size_type Spin_configuration_space_3<K, P, R>::size_of_predicates() const
{
    return m_predicates.size();
}

template<class K, class P, class R>
typename Spin_configuration_space_3<K, P, R>::Generic_predicate_const_iterator Spin_configuration_space_3<K, P, R>::generic_predicates_begin() const
{
    return m_generic_predicates.begin();
}

template<class K, class P, class R>
typename Spin_configuration_space_3<K, P, R>::Generic_predicate_const_iterator Spin_configuration_space_3<K, P, R>::generic_predicates_end() const
{
    return m_generic_predicates.end();
}

template<class K, class P, class R>
typename Spin_configuration_space_3<K, P, R>::Generic_predicate_size_type Spin_configuration_space_3<K, P, R>::size_of_generic_predicates() const
{
    return m_generic_predicates.size();
}
} // namespace CS