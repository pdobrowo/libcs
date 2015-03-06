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
#ifndef LIBCS_SPIN_RASTER_GRAPH_3_H
#define LIBCS_SPIN_RASTER_GRAPH_3_H

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
// Spin_raster_graph_3:
//
// Plain rasterizer, useful for debugging, topology checking and volume computation
//
// TODO:
//
// currently sampling is performed in [0; 1]^3
//
// reimplement with a solution based on uniform SO(3) sampling, presented in (Hopf fibration):
//
// Generating Uniform Incremental Grids on SO(3) Using the Hopf Fibration (pdf)
// Anna Yershova, Swati Jain, Steven M. LaValle, and Julie C. Mitchell,
// International Journal of Robotics Research, IJRR 2009
//
// Predicate_type is Predicate_tt_3 or Predicate_bb_3
template<class Kernel_,
         class Predicate_>
class Spin_raster_graph_3
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
        Parameters(size_t resolution);

        size_t resolution() const;

    private:
        size_t m_resolution;
    };

    // voxel data
    struct Voxel_data
    {
        typedef std::pair<Voxel_3<Voxel_data> *, Cover> Voxel_link;

        struct Heap_element
        {
            Voxel_link  link;
            double      distance;

            bool operator >(const Heap_element &other) const
            {
                return distance > other.distance;
            }

            Heap_element(Voxel_link link_,
                         double distance_)
                : link(link_),
                  distance(distance_)
            {
            }
        };

        typedef boost::heap::fibonacci_heap<
                    Heap_element,
                    boost::heap::compare<std::greater<Heap_element> > > Heap;

        typedef typename Heap::handle_type                              Heap_element_handle;

        // helper for connected components
        int                 component[Cover_COUNT];

        // routing: BFS and Dijkstra processing
        int                 common_color[Cover_COUNT];
        double              common_distance[Cover_COUNT];   // natural distance in S^3
        Voxel_link          common_parent[Cover_COUNT];
        Heap_element_handle dijkstra_heap_node[Cover_COUNT];
    };

    // voxel colors
    static const int COLOR_WHITE    = 0;
    static const int COLOR_GRAY     = 1;
    static const int COLOR_BLACK    = 2;

    // typedef for three dimensional voxels
    typedef Voxel_3<Voxel_data>             Voxel;
    typedef Voxel *                         Voxel_iterator;
    typedef typename Voxel_data::Voxel_link Voxel_link;

    // Dijkstra heap
    typedef typename Voxel_data::Heap                   Heap;
    typedef typename Voxel_data::Heap_element           Heap_element;
    typedef typename Voxel_data::Heap_element_handle    Heap_element_handle;

    // route container
    class Route
    {
    public:
        // invalid route
        Route();

        // valid route
        Route(const std::vector<Voxel_link> &nodes);

        bool                is_valid() const;
        Sample              evaluate(double t) const;

        std::vector<Sample> nodes() const;

    private:
        bool                    m_valid;
        std::vector<Voxel_link> m_nodes;

        log4cxx::LoggerPtr      m_logger;
    };

public:
    Spin_raster_graph_3(const std::vector<Predicate> &predicates,
                        const std::vector<Spin_quadric_3> &quadrics,
                        const Parameters &parameters);

    // routing
    Route                       find_route(const Sample &begin, const Sample &end);

    // geometry
    size_t                      resolution() const;
    const Voxel &               voxel(int u, int v, int w) const;

private:
    // internals
    size_t                      m_resolution;   // target Cartesian density in a unit cube; only pi/6 = ~52% will hit Spin(3) ball

    boost::scoped_array<Voxel>  m_voxels;       // m_density ^ 3

    size_t                      m_number_of_empty_voxels;   // out of all real voxels
    size_t                      m_number_of_full_voxels;    // out of all real voxels
    size_t                      m_number_of_mixed_voxels;   // out of all real voxels
    size_t                      m_number_of_voxels;         // resolution ^ 3
    size_t                      m_number_of_real_voxels;    // pi/6 * voxels
    size_t                      m_number_of_border_voxels;  // 4/3 pi * voxels

    log4cxx::LoggerPtr          m_logger;

    Voxel &                     voxel(int u, int v, int w);
    Voxel &                     locate_voxel(const Sample &sample);

    // main subroutines
    void                        collect_samples(const std::vector<Predicate> &predicates,
                                                const std::vector<Spin_quadric_3> &quadrics);

    // helpers
    void                        calculate_raster_coordinate(const std::vector<Spin_quadric_3> &quadrics,
                                                            const Sample &sample,
                                                            Coordinate &out) const;

    bool                        evaluate_predicate_list_at_raster_coordinate(const std::vector<Predicate> &predicates,
                                                                             const Coordinate &coordinate) const;

    // dfs
    void                        mark_connected_components_dfs(const Voxel_link &link,
                                                              int component,
                                                              bool value);

    std::pair<size_t, size_t>   mark_connected_components(Voxel_iterator begin,
                                                          Voxel_iterator end);

#if 0
    void                        mark_routing_distances_bfs(Voxel_link source);
#endif

    void                        mark_routing_distances_dijkstra(Voxel_link source);

    template<typename OutputIterator>
    void                        collect_backtrace_route(Voxel_link target, OutputIterator out);

    double                      natural_voxel_distance(Voxel_link begin, Voxel_link end);
    double                      cartesian_voxel_distance(Voxel_link begin, Voxel_link end);
};
} // namespace CS

#include "Spin_raster_graph_3.ipp"

#endif // LIBCS_SPIN_RASTER_GRAPH_3_H
