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
#include "Spin_null_graph_3.h"

namespace CS
{
template<class K, class P>
Spin_null_graph_3<K, P>::Parameters::Parameters()
{
}

template<class K, class P>
Spin_null_graph_3<K, P>::Spin_null_graph_3(const std::vector<Predicate> &predicates,
                                           const std::vector<Spin_quadric_3> &quadrics,
                                           const Parameters &parameters)
    : m_logger(log4cxx::Logger::getLogger("CS.Spin_null_graph_3"))
{
    (void)parameters;

    // check if everything is ok with predicates and quadrics
    assert(predicates.size() * SUB_PREDICATE_COUNT == quadrics.size());

    // ok to begin construction
    LOG4CXX_DEBUG(m_logger, "Creating spin null graph");
}

template<class K, class P>
typename Spin_null_graph_3<K, P>::Route Spin_null_graph_3<K, P>::find_route(const Sample &begin, const Sample &end)
{
    return Route();
}

template<class K, class P>
Spin_null_graph_3<K, P>::Route::Route()
    : m_valid(false),
      m_logger(log4cxx::Logger::getLogger("CS.Spin_null_graph_3::Route"))
{
}

#if 0
template<class K, class P>
Spin_null_graph_3<K, P>::Route::Route(const std::vector<Voxel_link> &nodes)
    : m_valid(true),
      m_nodes(nodes),
      m_logger(log4cxx::Logger::getLogger("CS.Spin_null_graph_3::Route"))
{
}
#endif

template<class K, class P>
bool Spin_null_graph_3<K, P>::Route::is_valid() const
{
    return m_valid;
}

template<class K, class P>
typename Spin_null_graph_3<K, P>::Sample Spin_null_graph_3<K, P>::Route::evaluate(double t) const
{
    return Sample();
}

template<class K, class P>
std::vector<typename Spin_null_graph_3<K, P>::Sample> Spin_null_graph_3<K, P>::Route::nodes() const
{
    std::vector<Sample> nodes;
    return nodes;
}

} // namespace CS
