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
#ifndef LIBCS_SPIN_CONFIGURATION_SPACE_3_H
#define LIBCS_SPIN_CONFIGURATION_SPACE_3_H

#include <CGAL/Kernel/global_functions.h>
#include "Predicate_tt_3_list_generator.h"
#include "Predicate_bb_3_list_generator.h"
#include "Predicate_h_3_list_generator.h"
#include <log4cxx/logger.h>
#include <boost/scoped_ptr.hpp>
#include <boost/foreach.hpp>
#include <cassert>

namespace CS
{
// Spin_configuration_space_3:
//
// A general interface for representing SO(3) configuration space
//
// Kernel         is a spin kernel
// Predicate      is a base predicate: Predicate_tt_3 or Predicate_bb_3
// Representation is Spin_cell_graph_3 or Spin_exact_graph_3
//
template<class Kernel_,
         class Predicate_,
         class Representation_>
class Spin_configuration_space_3
{
    typedef typename Kernel_::RT                    RT;
    typedef typename Kernel_::FT                    FT;
    typedef typename Extended_generator<FT>::Type   ExtendedFT;

    typedef typename Kernel_::Point_3               Point_3;
    typedef typename Kernel_::Vector_3              Vector_3;

    typedef typename Kernel_::Predicate_g_3         Predicate_g_3;

    typedef Predicate_                              Predicate;
    typedef typename Predicate::Sub_predicate       Sub_predicate;
    static const size_t SUB_PREDICATE_COUNT       = Predicate::SUB_PREDICATE_COUNT;

    typedef typename Kernel_::Spin_quadric_3        Spin_quadric_3;

    typedef Predicate_list_generator<Kernel_, Predicate> Generator;

    typedef std::vector<Spin_quadric_3>             Spin_quadric_container;
    typedef std::vector<Predicate>                  Predicate_container;
    typedef std::vector<Predicate_g_3>              General_predicate_container;

public:
    typedef Representation_                         Representation;
    typedef typename Representation::Parameters     Parameters;

    // configuration space point / sample
    typedef typename Representation::Sample         Sample;

    // configuration space route
    typedef typename Representation::Route          Route;

    // quadrics
    typedef typename Spin_quadric_container::const_iterator Spin_quadric_const_iterator;
    typedef typename Spin_quadric_container::size_type      Spin_quadric_size_type;

    // predicates
    typedef typename Predicate_container::const_iterator Predicate_const_iterator;
    typedef typename Predicate_container::size_type      Predicate_size_type;

    // general predicates
    typedef typename General_predicate_container::const_iterator General_predicate_const_iterator;
    typedef typename General_predicate_container::size_type      General_predicate_size_type;

    // ctor
    Spin_configuration_space_3();

    // scene from primitives
    template<typename RobotInputIterator_, typename ObstacleInputIterator_>
    void                                    create_from_scene(RobotInputIterator_ robot_begin, RobotInputIterator_ robot_end,
                                                              ObstacleInputIterator_ obstacle_begin, ObstacleInputIterator_ obstacle_end,
                                                              const Parameters &parameters = Parameters());

    // route finder
    Route                                   find_route(const Sample &begin, const Sample &end);

    // representation
    Representation &                        rep();
    const Representation &                  rep() const;

    // predicates
    Predicate_const_iterator                predicates_begin() const;
    Predicate_const_iterator                predicates_end() const;

    Predicate_size_type                     size_of_predicates() const;

    // general predicates
    General_predicate_const_iterator        general_predicates_begin() const;
    General_predicate_const_iterator        general_predicates_end() const;

    General_predicate_size_type             size_of_general_predicates() const;

    // quadrics
    Spin_quadric_const_iterator             spin_quadrics_begin() const;
    Spin_quadric_const_iterator             spin_quadrics_end() const;

    Spin_quadric_size_type                  size_of_spin_quadrics() const;

private:
    // predicate level
    Predicate_container                     m_predicates;
    General_predicate_container             m_general_predicates;

    // quadrics level
    Spin_quadric_container                  m_spin_quadrics;

    // representation: cell graph or exact graph
    boost::scoped_ptr<Representation>       m_representation;

    // logger
    log4cxx::LoggerPtr                      m_logger;
};

// standard representations
template<class Kernel_,
         class Predicate_>
class Spin_null_configuration_space_3
    : public Spin_configuration_space_3<Kernel_,
                                        Predicate_,
                                        typename Kernel_::template Spin_null_graph_3_generator<Predicate_>::Type>
{
};

template<class Kernel_,
         class Predicate_>
class Spin_cell_configuration_space_3
    : public Spin_configuration_space_3<Kernel_,
                                        Predicate_,
                                        typename Kernel_::template Spin_cell_graph_3_generator<Predicate_>::Type>
{
};

template<class Kernel_,
         class Predicate_>
class Spin_raster_configuration_space_3
    : public Spin_configuration_space_3<Kernel_,
                                        Predicate_,
                                        typename Kernel_::template Spin_raster_graph_3_generator<Predicate_>::Type>
{
};

template<class Kernel_,
         class Predicate_>
class Spin_exact_configuration_space_3
    : public Spin_configuration_space_3<Kernel_,
                                        Predicate_,
                                        typename Kernel_::template Spin_exact_graph_3_generator<Predicate_>::Type>
{
};

} // namespace CS

#include "Spin_configuration_space_3.ipp"

#endif // LIBCS_SPIN_CONFIGURATION_SPACE_3_H
