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
#ifndef LIBCS_SPIN_CELL_GRAPH_3_H
#define LIBCS_SPIN_CELL_GRAPH_3_H

#include "Spin_straight_sample_generator_3.h"
#include "Spin_QC_sample_generator_3.h"
#include "Coordinate.h"
#include "Cell_n.h"
#include "Config.h"
#include "Benchmark.h"      // TODO: add Timing.h
#include "Index_tree_n.h"
#include <memory>
#include <log4cxx/logger.h>
#include <functional>
#include <cstddef>
#include <stack>

namespace CS
{
// Spin_cell_graph_3:
//
// QC sample generator is ~10-100 times better than straight one
// when comparing quanity of neighbour cells generated
//
// Assumptions:
// * the graph must not contain duplicated quadrics!
// * zero quadrics are forbidden
//
// Predicate_type is Predicate_tt_3 or Predicate_bb_3
template<class Kernel_,
         class Predicate_,
#if 1
         class Spin_sample_generator_ = Spin_QC_sample_generator_3<Kernel_>
#else
         class Spin_sample_generator_ = Spin_straight_sample_generator_3<Kernel_>
#endif
         >
class Spin_cell_graph_3
{
    typedef typename Kernel_::FT                    FT;
    typedef typename Kernel_::Spin_quadric_3        Spin_quadric_3;
    typedef Spin_sample_generator_                  Spin_sample_generator;

    // sub-predicates
    typedef Predicate_                              Predicate;
    typedef typename Predicate::Sub_predicate       Sub_predicate;
    static const size_t SUB_PREDICATE_COUNT       = Predicate::SUB_PREDICATE_COUNT;

    // index for a predicate
    typedef bool Predicate_index[SUB_PREDICATE_COUNT];

public:
    typedef typename Spin_sample_generator::Sample  Sample;

    class Parameters
    {
    public:
        Parameters();
        Parameters(size_t sample_point_count);

        size_t sample_point_count() const;

    private:
        size_t m_sample_point_count;
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

private:
    // cell data
    struct Cell_data
    {
        // helper for connected components
        int component;
    };

    template<class Cell, typename Cell_iterator>
    void mark_connected_components_dfs(Cell_iterator iterator, int component, bool value);

    template<typename Cell_iterator>
    std::pair<size_t, size_t> mark_connected_components(Cell_iterator begin, Cell_iterator end);

public:
    // cell type is based on sample type
    typedef Cell_n<Sample, Cell_data>           Cell;

    // cell list
    typedef std::vector<Cell>                   Cell_list;

    // size type
    typedef typename Cell_list::size_type       size_type;
    typedef typename Cell_list::const_iterator  Cell_const_iterator;

    Spin_cell_graph_3(const std::vector<Predicate> &predicates,
                      const std::vector<Spin_quadric_3> &quadrics,
                      const Parameters &parameters);

    size_type                   size_of_cells() const;

    Cell_const_iterator         cells_begin() const;
    Cell_const_iterator         cells_end() const;

    Cell_const_iterator         add_sample_point(const Sample &sample);

    // routing
    Route                       find_route(const Sample &begin, const Sample &end);

private:
    // internals
    std::vector<Predicate>      m_predicates;   // original predicate list. n elements
    std::vector<size_t>         m_links;        // links to quadrics in compressed list, SUB_PREDICATE_COUNT * n elements
    std::vector<bool>           m_inverse_flag; // inverse flags to quadrics in compressed list, SUB_PREDICATE_COUNT * n elements
    std::vector<Spin_quadric_3> m_quadrics;     // compressed quadric list for sub-predicates, <= SUB_PREDICATE_COUNT * n elements

    Cell_list                   m_cells;
    size_t                      m_empty_cell_count;
    size_t                      m_full_cell_count;
    size_t                      m_coordinate_size;

    log4cxx::LoggerPtr          m_logger;

    // main subroutines
    size_t compress_duplicated_quadrics(const std::vector<Spin_quadric_3> &quadrics);

    void collect_sample_cell_list(
            size_t sample_point_count,
            size_t &outRawEmptyCellCount,
            size_t &outRawFullCellCount,
            Cell_list &out_raw_cell_list);

    void collapse_cell_list(Cell_list &raw_cell_list);

    size_t collect_cell_neighbour_information();

    // helpers
    void calculate_cell_coordinate(const Sample &sample, Coordinate &out) const;

    bool evaluate_predicate_list_at_cell_coordinate(const Coordinate &coordinate) const;

    void display_coordinate_pops_histogram(bool verbose);
    void display_neighbour_histogram();
    void display_best_sampled_cells(size_t sample_point_count, size_t top_count);

    void sort_cell_list(Cell_list &cells);

    // find a cell-based route in cell graph
    //
    // CellOutputIterator is an iterator to an output container of type Cell_const_iterator which holds the result
    // the method returns indicator whether a route was found or not
    template<typename OutputCellIterator>
    bool find_cell_route(Cell_const_iterator begin, Cell_const_iterator end, OutputCellIterator outputCellIterator);

    Route find_route_same_cell(typename Spin_cell_graph_3::Cell_const_iterator cell,
                               const Sample &begin,
                               const Sample &end) const;

    Route find_route_same_component(typename Spin_cell_graph_3::Cell_const_iterator begin_cell,
                                    typename Spin_cell_graph_3::Cell_const_iterator end_cell,
                                    const Sample &begin, const Sample &end) const;

    Route find_route_general(typename Spin_cell_graph_3::Cell_const_iterator begin_cell,
                             typename Spin_cell_graph_3::Cell_const_iterator end_cell,
                             const Sample &begin, const Sample &end) const;
};
} // namespace CS

#include "Spin_cell_graph_3.ipp"

#endif // LIBCS_SPIN_CELL_GRAPH_3_H
