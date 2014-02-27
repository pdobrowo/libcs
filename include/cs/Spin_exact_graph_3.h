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
#ifndef LIBCS_SPIN_EXACT_GRAPH_3_H
#define LIBCS_SPIN_EXACT_GRAPH_3_H

#include "Spin_3.h"             // FIXME: remove
#include "Benchmark.h"          // TODO: add Timing.h
#include <log4cxx/logger.h>
#include <boost/foreach.hpp>

namespace CS
{
// Spin_exact_graph_3:
//
// Exact arrangement of quadrics
//
// Predicate_type is Predicate_tt_3 or Predicate_bb_3
template<class K,
         class Predicate_>
class Spin_exact_graph_3
{
    typedef typename K::FT                          FT;
    typedef typename K::Spin_quadric_3              Spin_quadric_3;

    // sub-predicates
    typedef Predicate_                              Predicate;
    typedef typename Predicate::Sub_predicate       Sub_predicate;
    static const size_t SUB_PREDICATE_COUNT       = Predicate::SUB_PREDICATE_COUNT;

    // index for a predicate
    typedef bool Predicate_index[SUB_PREDICATE_COUNT];

public:
    typedef Spin_3<double>                          Sample;

    // primitives
    typedef typename K::Spin_qsic_3                 Spin_qsic_3;
    typedef typename K::Spin_qsip_3                 Spin_qsip_3;

    typedef Spin_qsic_3 *                           Qsic_handle;
    typedef Spin_qsip_3 *                           Qsip_handle;

private:
    // internal containers
    typedef std::vector<Spin_quadric_3>             Spin_quadric_container;
    typedef std::vector<Qsic_handle>                Qsic_container;
    typedef std::vector<Qsip_handle>                Qsip_container;    

public:
    // public iterators and handles
    typedef typename Spin_quadric_container::const_iterator Spin_quadric_const_iterator;
    typedef typename Qsic_container::const_iterator         Qsic_const_iterator;
    typedef typename Qsip_container::const_iterator         Qsip_const_iterator;

    typedef typename Spin_quadric_container::size_type      Spin_quadric_size_type;
    typedef typename Qsic_container::size_type              Qsic_size_type;
    typedef typename Qsip_container::size_type              Qsip_size_type;

    class Parameters
    {
    public:
        Parameters();
        Parameters(bool suppress_qsic_calculation, bool suppress_qsip_calculation);

        bool suppress_qsic_calculation() const;
        bool suppress_qsip_calculation() const;

    private:
        bool m_suppress_qsic_calculation;
        bool m_suppress_qsip_calculation;
    };

    // route container
    class Route
    {
    public:
        // invalid route
        Route();

        // valid route
        //Route(const std::vector<Voxel_link> &nodes);

        bool    is_valid() const;
        Sample  evaluate(double t) const;

    private:
        bool                    m_valid;
        //std::vector<Voxel_link> m_nodes;
    };

public:
    Spin_exact_graph_3(const std::vector<Predicate> &predicates,
                       const std::vector<Spin_quadric_3> &spin_quadrics,
                       const Parameters &parameters);

    ~Spin_exact_graph_3();

    // routing
    Route                       find_route(const Sample &begin, const Sample &end);

    // structure
    Spin_quadric_const_iterator spin_quadrics_begin() const;
    Spin_quadric_const_iterator spin_quadrics_end() const;

    Spin_quadric_size_type      size_of_spin_quadrics() const;

    Qsic_const_iterator         qsics_begin() const;
    Qsic_const_iterator         qsics_end() const;

    Qsic_size_type              size_of_qsics() const;

    Qsip_const_iterator         qsips_begin() const;
    Qsip_const_iterator         qsips_end() const;

    Qsip_size_type              size_of_qsips() const;

private:
    // internals
    Spin_quadric_container      m_spin_quadrics;
    Qsic_container              m_qsics;
    Qsip_container              m_qsips;

    log4cxx::LoggerPtr          m_logger;

    void release();
};
} // namespace CS

#include "Spin_exact_graph_3.ipp"

#endif // LIBCS_SPIN_EXACT_GRAPH_3_H
