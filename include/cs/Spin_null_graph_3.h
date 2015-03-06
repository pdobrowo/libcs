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
#ifndef LIBCS_SPIN_NULL_GRAPH_3_H
#define LIBCS_SPIN_NULL_GRAPH_3_H

#include "Coordinate.h"
#include "Spin_3.h"
#include "Voxel_3.h"
#include <log4cxx/logger.h>
#include <boost/scoped_array.hpp>
#include <boost/heap/fibonacci_heap.hpp>
#include <functional>
#include <limits>
#include <cmath>
#include <cassert>
#include <stack>
#include <vector>
#include <deque>

namespace CS
{
// Spin_null_graph_3:
//
// Null graph, used when representation is required
//
// Predicate_type is Predicate_tt_3 or Predicate_bb_3
template<class Kernel_,
         class Predicate_>
class Spin_null_graph_3
{
    typedef Kernel_ Kernel;

    typedef typename Kernel::FT                          FT;
    typedef typename Kernel::Spin_quadric_3              Spin_quadric_3;

    // sub-predicates
    typedef Predicate_                              Predicate;
    typedef typename Predicate::Sub_predicate       Sub_predicate;
    static const size_t SUB_PREDICATE_COUNT       = Predicate::SUB_PREDICATE_COUNT;

    // index for a predicate
    typedef bool Predicate_index[SUB_PREDICATE_COUNT];

public:
    typedef Spin_3<double>                          Sample;

    class Parameters
    {
    public:
        Parameters();
    };

    // route container
    class Route
    {
    public:
        // invalid route
        Route();

        // valid route
        //Route(const std::vector<Voxel_link> &nodes);

        bool                is_valid() const;
        Sample              evaluate(double t) const;

        std::vector<Sample> nodes() const;

    private:
        bool                    m_valid;
        //std::vector<Voxel_link> m_nodes;

        log4cxx::LoggerPtr      m_logger;
    };

public:
    Spin_null_graph_3(const std::vector<Predicate> &predicates,
                      const std::vector<Spin_quadric_3> &quadrics,
                      const Parameters &parameters);

    // routing
    Route                       find_route(const Sample &begin, const Sample &end);

private:
    log4cxx::LoggerPtr          m_logger;
};
} // namespace CS

#include "Spin_null_graph_3.ipp"

#endif // LIBCS_SPIN_NULL_GRAPH_3_H
