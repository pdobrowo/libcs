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
#include "Spin_cell_graph_3.h"

// use OMP
#define USE_PARALLEL    0

namespace CS
{
template<class K, class P, class G>
Spin_cell_graph_3<K, P, G>::Parameters::Parameters()
    : m_sample_point_count(10000)
{
}

template<class K, class P, class G>
Spin_cell_graph_3<K, P, G>::Parameters::Parameters(size_t sample_point_count)
    : m_sample_point_count(sample_point_count)
{
}

template<class K, class P, class G>
size_t Spin_cell_graph_3<K, P, G>::Parameters::sample_point_count() const
{
    return m_sample_point_count;
}

template<class K, class P, class G>
template<class Cell, typename Cell_iterator>
void Spin_cell_graph_3<K, P, G>::mark_connected_components_dfs(Cell_iterator iterator, int component, bool value)
{
    std::stack<Cell_iterator> stack;

    // initialize
    stack.push(iterator);

    // iterate
    while (!stack.empty())
    {
        // get next node
        Cell_iterator current = stack.top();
        stack.pop();

        // check if this cell has the proper value
        if (current->value() != value)
            continue;

        // check if this cell been already marked
        if (current->data().component != -1)
            continue;

        // mark this cell
        current->data().component = component;

        // propagate further
        typedef typename Cell::Edge_const_iterator Edge_const_iterator;

        for (Edge_const_iterator edgesIterator = current->edges_begin(); edgesIterator != current->edges_end(); ++edgesIterator)
            stack.push(*edgesIterator);
    }
}

template<class K, class P, class G>
template<typename Cell_iterator>
std::pair<size_t, size_t> Spin_cell_graph_3<K, P, G>::mark_connected_components(Cell_iterator begin, Cell_iterator end)
{
    typedef typename Cell_iterator::value_type Cell;

    // clear all components coordinates
    for (Cell_iterator iterator = begin; iterator != end; ++iterator)
        iterator->data().component = -1;

    // mark empty components
    std::pair<size_t, size_t> num_connected_components;
    int component = 0;

    for (Cell_iterator iterator = begin; iterator != end; ++iterator)
    {
        if (iterator->is_empty() && iterator->data().component == -1)
        {
            mark_connected_components_dfs<Cell>(iterator, component, false);
            ++component;
        }
    }

    // save empty component count
    num_connected_components.first = component;

    // mark full components
    for (Cell_iterator iterator = begin; iterator != end; ++iterator)
    {
        if (iterator->is_full() && iterator->data().component == -1)
        {
            mark_connected_components_dfs<Cell>(iterator, component, true);
            ++component;
        }
    }

    // save full component count
    num_connected_components.second = component - num_connected_components.first;

    return num_connected_components;
}

template<class K, class P, class G>
Spin_cell_graph_3<K, P, G>::Spin_cell_graph_3(const std::vector<Predicate> &predicates,
                                              const std::vector<Spin_quadric_3> &quadrics,
                                              const Parameters &parameters)
    : m_predicates(predicates),
      m_empty_cell_count(0),
      m_full_cell_count(0),
      m_logger(log4cxx::Logger::getLogger("CS.Spin_cell_graph_3"))
{
    // check if everything is ok with predicates and quadrics
    assert(predicates.size() * SUB_PREDICATE_COUNT == quadrics.size());

    // read parameters
    size_t sample_point_count = parameters.sample_point_count();

    // ok to begin construction
    unsigned long long timer_start, timer_end;

    LOG4CXX_DEBUG(m_logger, "Creating spin cell graph");

    // compress quadrics and create predicate-quadric link tab
    LOG4CXX_DEBUG(m_logger, "Compressing duplicated quadrics...");

    size_t duplicated_quadrics =
        compress_duplicated_quadrics(quadrics);

    LOG4CXX_DEBUG(m_logger, "Found " << duplicated_quadrics << " duplicated quadrics");

    // read coordinate size
    m_coordinate_size = m_quadrics.size();

    LOG4CXX_DEBUG(m_logger, "Compressed coordinate size is " << m_coordinate_size);

    LOG4CXX_DEBUG(m_logger, "Sampling " << sample_point_count << " spin points");

    timer_start = get_tick_count();

    size_t raw_empty_cell_count;
    size_t raw_full_cell_count;
    Cell_list raw_cell_list;

    collect_sample_cell_list(sample_point_count,
                             raw_empty_cell_count, raw_full_cell_count, raw_cell_list);

    timer_end = get_tick_count();

    LOG4CXX_INFO(m_logger, "Sampling result: " << raw_full_cell_count << " full cells, " << raw_empty_cell_count << " empty cells [" << (timer_end - timer_start) << " ms]");

    if (raw_empty_cell_count == 0)
        LOG4CXX_INFO(m_logger, "Configuration space is probably full");

    if (raw_full_cell_count == 0)
        LOG4CXX_INFO(m_logger, "Configuration space is probably empty");

#if 1
    // Note: sorting is not needed in index-tree version of cell collapsing algorithm

    // sort cells
    LOG4CXX_INFO(m_logger, "Sorting " << raw_cell_list.size() << " cells");

    sort_cell_list(raw_cell_list);

#endif

    // collapse cells
    LOG4CXX_INFO(m_logger, "Collapsing " << raw_cell_list.size() << " cells");

    collapse_cell_list(raw_cell_list);

    LOG4CXX_INFO(m_logger, "Collapsed to " << m_cells.size() << " different cells");
    LOG4CXX_INFO(m_logger, "Final structure: " << m_full_cell_count << " full cells, " << m_empty_cell_count << " empty cells ("
                 << std::fixed << std::setprecision(2) << (float(100 * m_empty_cell_count) / float(m_cells.size())) << "% empty)");

    // display best sampled cells
    display_best_sampled_cells(sample_point_count, 10);

    // collect neighbour information
    LOG4CXX_INFO(m_logger, "Collecting cell neighbour information");

    //display_coordinate_pops_histogram(false);
    //display_coordinate_pops_histogram(true);

    timer_start = get_tick_count();

    size_t num_neighbour_pairs =
        collect_cell_neighbour_information();

    timer_end = get_tick_count();

    LOG4CXX_INFO(m_logger, "Graph neighbour cell pairs: " << num_neighbour_pairs <<  " [" << (timer_end - timer_start) << " ms]");

    display_neighbour_histogram();

    // collect components
    LOG4CXX_INFO(m_logger, "Collecting connected components information");

    std::pair<int, int> num_connected_components =
        mark_connected_components(m_cells.begin(), m_cells.end());

    LOG4CXX_INFO(m_logger, "Number of empty connected components: " << num_connected_components.first);
    LOG4CXX_INFO(m_logger, "Number of full connected components: " << num_connected_components.second);
    LOG4CXX_INFO(m_logger, "Total number of connected components: " << (num_connected_components.first + num_connected_components.second));
}

template<class K, class P, class G>
typename Spin_cell_graph_3<K, P, G>::size_type Spin_cell_graph_3<K, P, G>::size_of_cells() const
{
    return m_cells.size();
}

template<class K, class P, class G>
typename Spin_cell_graph_3<K, P, G>::Cell_const_iterator Spin_cell_graph_3<K, P, G>::cells_begin() const
{
    return m_cells.begin();
}

template<class K, class P, class G>
typename Spin_cell_graph_3<K, P, G>::Cell_const_iterator Spin_cell_graph_3<K, P, G>::cells_end() const
{
    return m_cells.end();
}

template<class K, class P, class G>
typename Spin_cell_graph_3<K, P, G>::Cell_const_iterator Spin_cell_graph_3<K, P, G>::add_sample_point(const Sample &sample)
{
    LOG4CXX_DEBUG(m_logger, "Adding a custom sample point: " << sample);

    // lookup cell coordinate of new sample
    Coordinate coordinate;
    coordinate.resize(m_quadrics.size());

    calculate_cell_coordinate(sample, coordinate);

    // sample cell
    typename Cell_list::iterator iterator;

    // find corresponding cell
    for (iterator = m_cells.begin(); iterator != m_cells.end(); ++iterator)
    {
        if (coordinate == iterator->coordinate())
        {
            LOG4CXX_DEBUG(m_logger, "Sample point was added to an existing cell");

            // just add next sample point to the cell
            iterator->add_sample(sample);
            break;
        }
    }

    // still not found ?
    if (iterator == m_cells.end())
    {
        // create new cell
        LOG4CXX_DEBUG(m_logger, "Creating a new cell for the sample point");

        bool value = evaluate_predicate_list_at_cell_coordinate(coordinate);

        iterator = m_cells.insert(m_cells.end(), Cell(sample, coordinate, value));

        // update cell count information
        if (value)
            ++m_full_cell_count;
        else
            ++m_empty_cell_count;

        // update neighbour information
        // iterator is currently the last iterator, so omit it in enumeration
        for (typename Cell_list::iterator neighbourIterator = m_cells.begin(); neighbourIterator != iterator; ++neighbourIterator)
        {
            if (neighbour(*iterator, *neighbourIterator))
            {
                // add graph edges
                iterator->add_edge(neighbourIterator);
                neighbourIterator->add_edge(iterator);
            }
        }

        // update connectivity information
        //
        // FIXME: Reimplement with incremental version
        mark_connected_components(m_cells.begin(), m_cells.end());
    }

    // information
    LOG4CXX_DEBUG(m_logger, "Updated structure: " << m_full_cell_count << " full cells, " << m_empty_cell_count << " empty cells ("
                 << std::fixed << std::setprecision(2) << (float(100 * m_empty_cell_count) / float(m_cells.size())) << "% empty)");

    // return new cell iterator
    return iterator;
}

template<class K, class P, class G>
typename Spin_cell_graph_3<K, P, G>::Route Spin_cell_graph_3<K, P, G>::find_route(const Sample &begin, const Sample &end)
{
    // ensure that begin and end points are sampled in cell graph
    Cell_const_iterator begin_cell = add_sample_point(begin);
    Cell_const_iterator end_cell = add_sample_point(end);

    // check cells
    if (!begin_cell->is_empty())
    {
        LOG4CXX_INFO(m_logger, "Begin cell is forbidden: " << begin);
        return Route();
    }

    if (!end_cell->is_empty())
    {
        LOG4CXX_INFO(m_logger, "End cell is forbidden: " << end);
        return Route();
    }

    // the same cell
    if (begin_cell == end_cell)
        return find_route_same_cell(begin_cell, begin, end);

    // the same component
    if (begin_cell->data().component == end_cell->data().component)
        return find_route_same_component(begin_cell, end_cell, begin, end);

    // different components
    return find_route_general(begin_cell, end_cell, begin, end);
}

template<class K, class P, class G>
template<typename OutputCellIterator>
bool Spin_cell_graph_3<K, P, G>::find_cell_route(Cell_const_iterator begin, Cell_const_iterator end, OutputCellIterator outputCellIterator)
{
    // TODO: implement with Dijkstra
    return false;
}

template<class K, class P, class G>
typename Spin_cell_graph_3<K, P, G>::Route Spin_cell_graph_3<K, P, G>::find_route_same_cell(
        typename Spin_cell_graph_3::Cell_const_iterator cell,
        const Sample &begin, const Sample &end) const
{
    LOG4CXX_DEBUG(m_logger, "Find route: same cell");

    (void)cell;
    (void)begin;
    (void)end;

    return Route();
}

template<class K, class P, class G>
typename Spin_cell_graph_3<K, P, G>::Route Spin_cell_graph_3<K, P, G>::find_route_same_component(
        typename Spin_cell_graph_3::Cell_const_iterator begin_cell,
        typename Spin_cell_graph_3::Cell_const_iterator end_cell,
        const Sample &begin, const Sample &end) const
{
    LOG4CXX_DEBUG(m_logger, "Find route: same component");

    (void)begin_cell;
    (void)end_cell;
    (void)begin;
    (void)end;

    return Route();
}

template<class K, class P, class G>
typename Spin_cell_graph_3<K, P, G>::Route Spin_cell_graph_3<K, P, G>::find_route_general(
        typename Spin_cell_graph_3::Cell_const_iterator begin_cell, typename Spin_cell_graph_3::Cell_const_iterator end_cell,
        const Sample &begin, const Sample &end) const
{
    LOG4CXX_DEBUG(m_logger, "Find route: general");

    (void)begin_cell;
    (void)end_cell;
    (void)begin;
    (void)end;

    return Route();
}

template<class K, class P, class G>
size_t Spin_cell_graph_3<K, P, G>::compress_duplicated_quadrics(const std::vector<Spin_quadric_3> &quadrics)
{
    // take a copy of quadrics in order to normalize them
    std::vector<Spin_quadric_3> normalized_quadrics = quadrics;
    std::vector<Spin_quadric_3> normalized_inversed_quadrics = quadrics;

    // normalize quadrics
    for (size_t i = 0; i < normalized_quadrics.size(); ++i)
    {
        FT scale = normalized_quadrics[i].max_abs_coefficient();
        normalized_quadrics[i].scale(FT(1) / scale);
    }

    for (size_t i = 0; i < normalized_inversed_quadrics.size(); ++i)
    {
        normalized_inversed_quadrics[i].inverse();

        FT scale = normalized_inversed_quadrics[i].max_abs_coefficient();
        normalized_inversed_quadrics[i].scale(FT(1) / scale);
    }

    // prepare quadric links
    m_links.reserve(quadrics.size());

    // prepate quadric inverse flags
    m_inverse_flag.reserve(quadrics.size());

    // a list of compressed quadrics
    std::vector<Spin_quadric_3> compressed_quadrics;

    // for each quadric, find first quadric in list, that is the same as given (or its inverse)
    // in most cases the enumeration will finish when the given quadric is reached
    // store links and compressed quadrics during search
    //
    // result: neither the same quadric nor the inverse of quadric will occur more that once in final list
    size_t duplicated_quadrics = 0;

    const double MAX_ABS_DIFFERENCE = 10e-10;

    for (size_t current = 0; current < quadrics.size(); ++current)
    {
        bool found = false;

        // search in compressed quadrics
        for (size_t i = 0; i < compressed_quadrics.size(); ++i)
        {
            // check normal quadric
            if (max_abs_difference(normalized_quadrics[current], compressed_quadrics[i]) < MAX_ABS_DIFFERENCE)
            {
                m_links.push_back(i);
                m_inverse_flag.push_back(false);
                ++duplicated_quadrics;
                found = true;
                break;
            }

            // check inversed quadric
            if (max_abs_difference(normalized_inversed_quadrics[current], compressed_quadrics[i]) < MAX_ABS_DIFFERENCE)
            {
                m_links.push_back(i);
                m_inverse_flag.push_back(true);
                ++duplicated_quadrics;
                found = true;
                break;
            }
        }

        if (!found)
        {
            m_links.push_back(compressed_quadrics.size());
            m_inverse_flag.push_back(false);

            compressed_quadrics.push_back(normalized_quadrics[current]);
        }
    }

    // swap compressed quadrics with original list
    m_quadrics = compressed_quadrics;

    return duplicated_quadrics;
}

template<class K, class P, class G>
void Spin_cell_graph_3<K, P, G>::collect_sample_cell_list(
        size_t sample_point_count,
        size_t &out_raw_empty_cell_count,
        size_t &out_raw_full_cell_count,
        Cell_list &out_raw_cell_list)
{
    // reset counters
    out_raw_empty_cell_count = 0;
    out_raw_full_cell_count = 0;

    // prepare ouput cell list
    out_raw_cell_list.clear();
    out_raw_cell_list.resize(sample_point_count);

#if USE_PARALLEL

    // multi-threaded (OMP)
    int max_threads = omp_get_max_threads();

    LOG4CXX_DEBUG(m_logger, "Multi-threaded collect: " << max_threads << " threads");

    boost::scoped_array<Sample>                 samples(new Sample[max_threads]);
    boost::scoped_array<Spin_sample_generator>  generators(new Spin_sample_generator[max_threads]);
    boost::scoped_array<Coordinate>             coordinates(new Coordinate[max_threads]);

    for (int i = 0; i < max_threads; ++i)
        coordinates[i].resize(m_quadrics.size());

    boost::uint32_t empty = 0, full = 0;

    #pragma omp parallel for
    for (size_t z = 0; z < sample_point_count; ++z)
    {
        // choose next random sample
        generators[omp_get_thread_num()](m_quadrics.begin(), m_quadrics.end(), samples[omp_get_thread_num()]);

        // calculate sample coodrinate
        calculate_cell_coordinate(samples[omp_get_thread_num()], coordinates[omp_get_thread_num()]);

        // evaluate predicate value at index
        if (evaluate_predicate_list_at_cell_coordinate(coordinates[omp_get_thread_num()]))
        {
            // collect statistics
            boost::interprocess::detail::atomic_inc32(&full);
            out_raw_cell_list[z] = Cell(samples[omp_get_thread_num()], coordinates[omp_get_thread_num()], true);
        }
        else
        {
            // collect statistics
            boost::interprocess::detail::atomic_inc32(&empty);
            out_raw_cell_list[z] = Cell(samples[omp_get_thread_num()], coordinates[omp_get_thread_num()], false);
        }
    }

    // read atomic values
    out_raw_full_cell_count = static_cast<size_t>(boost::interprocess::detail::atomic_read32(&full));
    out_raw_empty_cell_count = static_cast<size_t>(boost::interprocess::detail::atomic_read32(&empty));

#else // USE_PARALLEL

    // single-threaded
    Sample sample;
    Spin_sample_generator generator;

    LOG4CXX_DEBUG(m_logger, "Single-threaded collect");

    Coordinate coordinate;
    coordinate.resize(m_quadrics.size());

    for (size_t z = 0; z < sample_point_count; z++)
    {
        // choose next random sample
        generator(m_quadrics.begin(), m_quadrics.end(), sample);

        // calculate sample coodrinate
        calculate_cell_coordinate(sample, coordinate);

        // evaluate predicate value at index
        if (evaluate_predicate_list_at_cell_coordinate(coordinate))
        {
            // collect statistics
            ++out_raw_full_cell_count;
            out_raw_cell_list[z] = Cell(sample, coordinate, true);
        }
        else
        {
            // collect statistics
            ++out_raw_empty_cell_count;
            out_raw_cell_list[z] = Cell(sample, coordinate, false);
        }
    }

#endif // USE_PARALLEL
}

template<class K, class P, class G>
size_t Spin_cell_graph_3<K, P, G>::collect_cell_neighbour_information()
{
    // Collect algorithm:
    //
    // 1: naive     O(n * n * l)
    // 2: optimized O(n * l * l)
    // 3: optimal   O(n * l)
    int neighbour_collect_algorithm = Config::neighbour_collect_algorithm();

    LOG4CXX_DEBUG(m_logger, "Config: neighbour_collect_algorithm = " << neighbour_collect_algorithm);

    if (neighbour_collect_algorithm == 3)
    {
        // O(l * n) optimal algorithm
        //
        // The algorithm is as follows:
        //
        // (preprocess) store word list in a tree for fast lookup
        //              maximum tree height is l
        //              O(n * l)
        // (iterate)    .
        // W grupach umieszczone sa te elementy ktore reprezentuja sciezki w drzewie bedace takimi samymi liczac od liscia.
        // Zadaniem grup i ich podzialow jest utrzymanie zbiorow sciezek bedacych tymi samymi koncowkami w indeksach.
        // Algorym iteruje po poziomach w drzewie, od lisci do korzenia. W pierwszym etapie mam jedna grupe koncowek, w drugim
        // etapie grupa ta dzieli sie na dwie - koncowki prawe i lewe. W elemencie grupy zawsze zachowany jest iterator na indeks (trzymany w lisciu)
        // ktory reprezentowany jest przez aktualna koncowke. Grupy jednoelelementowe mozna wyrzucac, gdyz nigdy sie nie sparuja.
        // Sprawdzanie par na danym poziomie jest liniowe do n. Najpierw oznaczamy rodzica dane rodzicow elmentow grupy na zero. Pozniej
        // ponownie przegladamy elementy i zaznczamy w danych rodzicach dane na aktualny elemnt lub parujemy, jesli juz takie dane byly zaznaczone.
        //
        size_t num_neighbour_pairs = 0;

        // preprocess - build up coordinate tree
        typedef Index_tree_node_visitor_n<Cell> Index_tree_node_visitor;
        typedef typename Index_tree_node_visitor::Index_tree_node Index_tree_node;
        //typedef Index_tree_node * Index_tree_node_ptr;
        typedef Multi_list_n<Index_tree_node_visitor> Multi_list;
        typedef Multi_list_rope_n<Index_tree_node_visitor> Multi_list_rope;
        typedef Multi_list_link_n<Index_tree_node_visitor> Multi_list_link;
        typedef typename Cell::Handle Handle;
        typedef typename Index_tree_node::Pool Pool;
        typedef typename Pool::Offset Offset;

        LOG4CXX_TRACE(m_logger, "Creating index tree pool");

        Pool pool(1 + m_cells.size() * m_coordinate_size);

        LOG4CXX_TRACE(m_logger, "Pool size is: " << std::fixed << std::setprecision(2) << (pool.memsize() / (1024 * 1024)) << " MB");

        LOG4CXX_TRACE(m_logger, "Building cell coordinate lookup tree");

        // index tree
        Offset cell_index_tree = Index_tree_node::create_root(&pool);

        // groups and initial full group
        Multi_list ropes(m_cells.size() * 2); // * 2 is for each one temporary left group during phase 2. TODO: Optimize

        LOG4CXX_TRACE(m_logger, "Ropes size is: " << std::fixed << std::setprecision(2) << (ropes.memsize() / (1024 * 1024)) << " MB");

        Multi_list_rope *initial_rope_iterator = ropes.alloc_rope();
        ropes.push_rope(initial_rope_iterator);

        Offset node;
        Multi_list_link *link;

        // create index tree and ropes at once
        for (typename Cell_list::iterator cell_iterator = m_cells.begin(); cell_iterator != m_cells.end(); ++cell_iterator)
        {
            // store node in tree
            node = Index_tree_node::store(cell_iterator->coordinate().begin(),
                                          cell_iterator->coordinate().end(),
                                          cell_iterator,
                                          &pool,
                                          cell_index_tree);

            // alloc link in rope
            link = ropes.alloc_link();

            // setup link
            link->set_value(Index_tree_node_visitor(node, cell_iterator));

            // add link to rope
            initial_rope_iterator->push_link(link);
        }

        // iterate
        LOG4CXX_TRACE(m_logger, "Searching for neighbour cell indexes");

        for (size_t l = 0; l < m_coordinate_size; ++l)
        {
            // first phase:
            //
            // for each group, search for siblings and output pairs
            // sibling search is realised in the following steps:
            // - all actions are taken only on current group
            // - for each element mark parent node data as null
            // - for each element
            // - - if parent node data is set, this is a sibling - output a pair where parent data points to a sibling element
            // - - otherwise set parent node data to this element visitor
            //
            // zeroing step is required because target node can be visible for more than one element

            // all actions are taken only on current group
            for (Multi_list_rope *rope_iterator = ropes.ropes_head(); rope_iterator != 0; rope_iterator = rope_iterator->next())
            {
                // for each element mark parent node data as null
                for (Multi_list_link *link_iterator = rope_iterator->links_head(); link_iterator != 0; link_iterator = link_iterator->next())
                    pool.ref(pool.ref(link_iterator->value().node)->parent())->set_data(0);

                // for each element
                for (Multi_list_link *link_iterator = rope_iterator->links_head(); link_iterator != 0; link_iterator = link_iterator->next())
                {
                    Index_tree_node_visitor *sibling_visitor = pool.ref(pool.ref(link_iterator->value().node)->parent())->data();

                    // if parent node data is set, this is a sibling - output a pair where parent data points to a sibling element
                    if (sibling_visitor)
                    {
                        // retrieve neighbour cell handles
                        const Handle &cell_handle = link_iterator->value().handle;
                        const Handle &other_cell_handle = sibling_visitor->handle;

                        assert(cell_handle != other_cell_handle);

                        ++num_neighbour_pairs;

                        // add graph edges
                        cell_handle->add_edge(other_cell_handle);
                        other_cell_handle->add_edge(cell_handle);
                    }
                    else // otherwise set parent node data to this element visitor
                    {
                        pool.ref(pool.ref(link_iterator->value().node)->parent())->set_data(&link_iterator->value());
                    }
                }
            }

            // second phase:
            //
            // split each group and advance level
            // a split of group is done with respect to node side (left or right)
            // after each split, current level nodes are set
            // all actions are taken only on current group
            for (Multi_list_rope *rope_iterator = ropes.ropes_head(); rope_iterator != 0; rope_iterator = rope_iterator->next())
            {
                // push empty group, so that copying is avoided
                Multi_list_rope *left_rope_iterator = ropes.alloc_rope();
                ropes.push_rope(left_rope_iterator);

                // splice left elements into new list
                Multi_list_link *link_iterator = rope_iterator->links_head();
                Multi_list_link *next_link_iterator;

                while (link_iterator != 0)
                {
                    // check whether element is left or right
                    bool is_left = Index_tree_node::is_left(&pool, link_iterator->value().node);

                    // update current level node
                    link_iterator->value().node = pool.ref(link_iterator->value().node)->parent();

                    // advance element prior to splice
                    next_link_iterator = link_iterator->next();

                    // splice if the element is left
                    if (is_left)
                    {
                        rope_iterator->pop_link(link_iterator);
                        left_rope_iterator->push_link(link_iterator);
                    }

                    link_iterator = next_link_iterator;
                }
            }

            // third phase:
            //
            // optimize traversal by removing empty or one-element groups (degenerated groups)
            Multi_list_rope *rope_iterator = ropes.ropes_head();
            Multi_list_rope *next_rope_iterator;

            while (rope_iterator != 0)
            {
                next_rope_iterator = rope_iterator->next();

                // empty or one-link ropes are degenerated
                if (!rope_iterator->links_head() || !rope_iterator->links_head()->next())
                {
                    ropes.pop_rope(rope_iterator);
                    ropes.release_rope(rope_iterator);
                }

                rope_iterator = next_rope_iterator;
            }
        }

        LOG4CXX_TRACE(m_logger, "Pool final usage: " << pool.usage() << " % out of total " << std::fixed << std::setprecision(2) << (pool.memsize() / (1024 * 1024)) << " MB");
        LOG4CXX_TRACE(m_logger, "Pool final allocted size is: " << std::fixed << std::setprecision(2) << (pool.allocated_memsize() / (1024 * 1024)) << " MB");

        return num_neighbour_pairs;
    }
    else if (neighbour_collect_algorithm == 2)
    {
        // O(l * l * n) optimized algorithm
        //
        // The algorithm is as follows:
        //
        // (preprocess) store word list in a tree for fast lookup
        //              maximum tree height is l
        //              O(n * l)
        // (iterate)    for each word, and each bit, flip the bit and check whether
        //              a new word exists in the tree
        //              O(n * l * l)
        size_t num_neighbour_pairs = 0;

        // preprocess - build up coordinate tree
        typedef Index_tree_node_n<Cell, Empty_node_data> Index_tree_node;
        typedef typename Index_tree_node::Pool Pool;
        typedef typename Pool::Offset Offset;

        LOG4CXX_TRACE(m_logger, "Creating index tree pool");

        Pool pool(1 + m_cells.size() * m_coordinate_size);

        LOG4CXX_TRACE(m_logger, "Pool size is: " << std::fixed << std::setprecision(2) << (pool.memsize() / (1024 * 1024)) << " MB");

        LOG4CXX_TRACE(m_logger, "Building cell coordinate lookup tree");

        Offset cell_index_tree = Index_tree_node::create_root(&pool);

        for (typename Cell_list::iterator cell_iterator = m_cells.begin(); cell_iterator != m_cells.end(); ++cell_iterator)
        {
            Index_tree_node::store(cell_iterator->coordinate().begin(),
                                   cell_iterator->coordinate().end(),
                                   cell_iterator,
                                   &pool,
                                   cell_index_tree);
        }

        // iterate
        LOG4CXX_TRACE(m_logger, "Searching for neighbour cell indexes");

        for (typename Cell_list::iterator cell_iterator = m_cells.begin(); cell_iterator != m_cells.end(); ++cell_iterator)
        {
            Coordinate neighbour_coordinate = cell_iterator->coordinate();

            for (size_t bitIndex = 0; bitIndex < m_coordinate_size; ++bitIndex)
            {
                // handle each pair only once - a switch from 0 to 1
                if (neighbour_coordinate[bitIndex])
                    continue;

                // flip coordinate bit
                //neighbour_coordinate[bitIndex].flip();
                neighbour_coordinate[bitIndex] = !neighbour_coordinate[bitIndex];

                // check neighbour index
                boost::optional<typename Cell::Handle> other_cell_handle =
                        Index_tree_node::contains(neighbour_coordinate.begin(),
                                                  neighbour_coordinate.end(),
                                                  &pool,
                                                  cell_index_tree);

                if (other_cell_handle)
                {
                    ++num_neighbour_pairs;

                    // add graph edges
                    cell_iterator->add_edge(*other_cell_handle);
                    (*other_cell_handle)->add_edge(cell_iterator);
                }

                // flip-back index
                //neighbour_coordinate[bitIndex].flip();
                neighbour_coordinate[bitIndex] = !neighbour_coordinate[bitIndex];
            }
        }

        LOG4CXX_TRACE(m_logger, "Pool final usage: " << pool.usage() << " % out of total " << std::fixed << std::setprecision(2) << (pool.memsize() / (1024 * 1024)) << " MB");
        LOG4CXX_TRACE(m_logger, "Pool final allocted size is: " << std::fixed << std::setprecision(2) << (pool.allocated_memsize() / (1024 * 1024)) << " MB");

        return num_neighbour_pairs;
    }
    else if (neighbour_collect_algorithm == 1)
    {
        // O(l * n * n) naive algorithm
        //
        // for each coodinate pair, check if these are neighbours
        size_t num_neighbour_pairs = 0;

        // this routine is collecting empty and full cell connections too
        for (typename Cell_list::iterator first_cell_iterator = m_cells.begin(); first_cell_iterator != m_cells.end(); ++first_cell_iterator)
        {
            // do not count pairs twice: enumerate only from next cell
            typename Cell_list::iterator next_cell_iterator = first_cell_iterator;
            ++next_cell_iterator;

            // check neighbours
            for (typename Cell_list::iterator second_cell_iterator = next_cell_iterator; second_cell_iterator != m_cells.end(); ++second_cell_iterator)
            {
                if (neighbour(*first_cell_iterator, *second_cell_iterator))
                {
                    ++num_neighbour_pairs;

                    // add graph edges
                    first_cell_iterator->add_edge(second_cell_iterator);
                    second_cell_iterator->add_edge(first_cell_iterator);
                }
            }
        }

        return num_neighbour_pairs;
    }
    else
    {
        assert(0 && "invalid algorithm");
    }

    return 0;
}

template<class K, class P, class G>
void Spin_cell_graph_3<K, P, G>::collapse_cell_list(
        Cell_list &raw_cell_list)
{
#if 1

    // Standard algorithm
    //
    // Note: the following algorithm assumes that cells are sorted

    for (typename Cell_list::iterator cell_iterator = raw_cell_list.begin(); cell_iterator != raw_cell_list.end(); ++cell_iterator)
    {
        Cell accumulated_cell = *cell_iterator;

        // start to accumulate cells
        typename Cell_list::iterator next_cell_iterator = cell_iterator;
        ++next_cell_iterator;

        // accumulate until next cell index is different from current
        while (next_cell_iterator != raw_cell_list.end() && *cell_iterator == *next_cell_iterator)
        {
            // this cell can be accumulated - only samples are different, values must be equal
            assert(cell_iterator->value() == next_cell_iterator->value());

            // there should one sample only for the next cell
            assert(next_cell_iterator->samples_size() == 1);

            // add next sample for current cell
            accumulated_cell.add_sample(*next_cell_iterator->samples_begin());

            // move on iterators
            ++cell_iterator;
            ++next_cell_iterator;
        }

        // store accumulated cell
        m_cells.push_back(accumulated_cell);

        // update statistics
        if (accumulated_cell.is_empty())
            ++m_empty_cell_count;
        else
            ++m_full_cell_count;
    }

    // clear raw cells
    raw_cell_list.clear();

#else

    // Optimized algorithm
    //
    // Use size_t instantiation of index tree to lookup repeating cell coordinates.
    // Cell iterator is unused.
    //
    // TODO: move cell iterator to templated cell data too.
    // TODO: memory usage is very high due to unknown final result size
    //
    // Important: Index tree node' cell data must no contain an iterator to vector element because it is MUTABLE here.
    //            This is because the final cell list is growing and vector resize may occur thus invalidating all iterators

    // preprocess - build up coordinate tree
    typedef Index_tree_node_n<Cell, size_t> Index_tree_node;
    typedef typename Index_tree_node::Pool Pool;

    LOG4CXX_TRACE(m_logger, "Creating reduced index tree pool");

    Pool pool(1 + raw_cell_list.size() * m_coordinate_size);

    LOG4CXX_TRACE(m_logger, "Reduced pool size is: " << std::fixed << std::setprecision(2) << (pool.memsize() / (1024 * 1024)) << " MB");

    LOG4CXX_TRACE(m_logger, "Building reduced cell coordinate lookup tree");

    // index tree
    Index_tree_node *cell_index_tree = Index_tree_node::create_root(&pool);
    std::pair<Index_tree_node *, bool> result;

    // create index tree and collapse cells
    for (typename Cell_list::iterator cell_iterator = raw_cell_list.begin(); cell_iterator != raw_cell_list.end(); ++cell_iterator)
    {
        // there should one sample per cell
        assert(cell_iterator->samples_size() == 1);

        // store node and interpret result
        result = cell_index_tree->check_and_store(cell_iterator->coordinate().begin(), cell_iterator->coordinate().end(), &pool);

        // new node
        if (result.second)
        {
            // add cell to collapsed list
            m_cells.push_back(*cell_iterator);

            // store new cell iterator
            result.first->set_data(m_cells.size() - 1); // size() must do it in O(1) !

            // update statistics
            if (cell_iterator->is_empty())
                ++m_empty_cell_count;
            else
                ++m_full_cell_count;
        }
        else
        {
            // node already exists
            // just add more samples to cell
            m_cells[result.first->data()].add_sample(*cell_iterator->samples_begin());
        }
    }

    // clear raw cells
    raw_cell_list.clear();

#endif
}

template<class K, class P, class G>
void Spin_cell_graph_3<K, P, G>::calculate_cell_coordinate(
        const Sample &sample, Coordinate &out) const
{
    assert(out.size() == m_quadrics.size());

    // FIXME: handle bad hits of quadrics
    CGAL::Sign sign;

    for (size_t i = 0; i < m_quadrics.size(); ++i)
    {
        sign = m_quadrics[i].evaluate_sign(sample);
        //double dsign = m_quadrics[i].evaluate(sample);

        // choosing a sample on a quadric is a bad idea
        if (sign == 0)
        //if (fabs(dsign) < 10e-6)
        {
//            LOG4CXX_WARN(m_logger, "Sample point " << sample << " lies on a quadric");

            // TODO: special code is needed to handle degenerate collision cases
        }

        //sign = CGAL::sign(dsign);

        out[i] = (sign >= 0);
    }
}

template<class K, class P, class G>
bool Spin_cell_graph_3<K, P, G>::evaluate_predicate_list_at_cell_coordinate(
        const Coordinate &coordinate) const
{
    // assert number of links to compressed quadric list
    assert(m_links.size() == SUB_PREDICATE_COUNT * m_predicates.size());

    // predicate collision sentence is predicates[0](spin) || predicates[1](spin) || ... || predicates[N](spin)
    // evaluate until result is known
    for (size_t i = 0; i < m_predicates.size(); ++i)
    {
        Predicate_index index;

        for (size_t c = 0; c < SUB_PREDICATE_COUNT; ++c)
        {
            index[c] = coordinate[m_links[SUB_PREDICATE_COUNT * i + c]];

            if (m_inverse_flag[SUB_PREDICATE_COUNT * i + c])
                index[c] = !index[c];
        }

        bool value = m_predicates[i].evaluate(index);

        if (value)
            return true; // collision
    }

    // no collision
    return false;
}

template<class K, class P, class G>
void Spin_cell_graph_3<K, P, G>::display_coordinate_pops_histogram(bool verbose)
{
    LOG4CXX_DEBUG(m_logger, "Pops histogram:");

    boost::scoped_array<size_t> counts(new size_t[m_coordinate_size]);

    for (size_t i = 0; i < m_coordinate_size; ++i)
        counts[i] = 0;

    size_t index = 0;

    for (typename Cell_list::iterator cell_iterator = m_cells.begin(); cell_iterator != m_cells.end(); ++cell_iterator)
    {
        Coordinate coordinate = cell_iterator->coordinate();
        size_t count = 0;

        std::string bits;

        for (size_t bit = 0; bit < m_coordinate_size; ++bit)
        {
            if (coordinate[bit])
                ++count;

            if (verbose)
                bits += (coordinate[bit] ? " 1" : " 0");
        }

        if (verbose)
            LOG4CXX_DEBUG(m_logger, "Cell(" << index << ") = " << bits);

        ++counts[count];
        ++index;
    }

    for (size_t i = 0; i < m_coordinate_size; ++i)
        if (counts[i] != 0)
            LOG4CXX_DEBUG(m_logger, "Pops(" << i << ") = " << counts[i]);

    // check cell coordinate compression
    size_t reduce = 0;

    for (size_t bit = 0; bit < m_coordinate_size; ++bit)
    {
        Coordinate::value_type expected_bit = m_cells.front().coordinate()[bit];
        bool success = true;

        for (typename Cell_list::iterator cell_iterator = m_cells.begin(); cell_iterator != m_cells.end(); ++cell_iterator)
        {
            const Coordinate &coordinate = cell_iterator->coordinate();

            if (coordinate[bit] != expected_bit)
            {
                success = false;
                break;
            }
        }

        if (success)
            ++reduce;
    }

    LOG4CXX_DEBUG(m_logger, "Coordinate description may be reduced by " << reduce << " bits");
}

template<class K, class P, class G>
void Spin_cell_graph_3<K, P, G>::display_neighbour_histogram()
{
#if 0

    // Memory consuming implementation
    size_t num_cells = m_cells.size();
    size_t max_edges_index = 1 + num_cells * (num_cells - 1) / 2; // all connected pair plus one

    boost::scoped_array<size_t> counts(new size_t[max_edges_index + 1]);

    for (size_t i = 0; i <= max_edges_index; ++i)
        counts[i] = 0;

    for (typename Cell_list::iterator cell_iterator = m_cells.begin(); cell_iterator != m_cells.end(); ++cell_iterator)
    {
        size_t index = cell_iterator->edges_size();
        assert(index <= max_edges_index);

        ++counts[index];
    }

    for (size_t i = 0; i <= max_edges_index; ++i)
        if (counts[i] != 0)
            LOG4CXX_DEBUG(m_logger, "Edges(" << i << ") = " << counts[i]);

#else

    // accumulator version
    typedef std::map<size_t, size_t> Counts;

    Counts general_counts;
    Counts empty_counts;
    Counts full_counts;

    for (typename Cell_list::iterator cell_iterator = m_cells.begin(); cell_iterator != m_cells.end(); ++cell_iterator)
    {
        // general
        size_t general_index = cell_iterator->edges_size();
        size_t empty_index = 0;
        size_t full_index = 0;

        typedef typename Cell::Edge_const_iterator Edge_const_iterator;

        for (Edge_const_iterator edges_iterator = cell_iterator->edges_begin(); edges_iterator != cell_iterator->edges_end(); ++edges_iterator)
        {
            // add empty
            if (cell_iterator->is_empty() && (*edges_iterator)->is_empty())
                ++empty_index;

            // add full
            if (cell_iterator->is_full() && (*edges_iterator)->is_full())
                ++full_index;
        }

        // general
        Counts::iterator iterator = general_counts.find(general_index);

        if (iterator != general_counts.end())
            ++iterator->second;
        else
            general_counts.insert(std::make_pair(general_index, 1));

        // empty
        if (cell_iterator->is_empty())
        {
            Counts::iterator empty_iterator = empty_counts.find(empty_index);

            if (empty_iterator != empty_counts.end())
                ++empty_iterator->second;
            else
                empty_counts.insert(std::make_pair(empty_index, 1));
        }

        // full
        if (cell_iterator->is_full())
        {
            Counts::iterator full_iterator = full_counts.find(full_index);

            if (full_iterator != full_counts.end())
                ++full_iterator->second;
            else
                full_counts.insert(std::make_pair(full_index, 1));
        }
    }

    LOG4CXX_DEBUG(m_logger, "All cells histogram:");

    for (Counts::const_iterator i = general_counts.begin(); i != general_counts.end(); ++i)
        LOG4CXX_DEBUG(m_logger, "#" << i->first << " degree: " << i->second << " cells");

    LOG4CXX_DEBUG(m_logger, "Empty cells histogram:");

    for (Counts::const_iterator i = empty_counts.begin(); i != empty_counts.end(); ++i)
        LOG4CXX_DEBUG(m_logger, "#" << i->first << " degree: " << i->second << " cells");

    LOG4CXX_DEBUG(m_logger, "Full cells histogram:");

    for (Counts::const_iterator i = full_counts.begin(); i != full_counts.end(); ++i)
        LOG4CXX_DEBUG(m_logger, "#" << i->first << " degree: " << i->second << " cells");

#endif
}

template<class K, class P, class G>
void Spin_cell_graph_3<K, P, G>::display_best_sampled_cells(size_t sample_point_count, size_t top_count)
{
    std::vector<typename Cell::Sample_size_type> sizes;
    sizes.reserve(m_cells.size());

    for (typename Cell_list::iterator cell_iterator = m_cells.begin(); cell_iterator != m_cells.end(); ++cell_iterator)
        sizes.push_back(cell_iterator->samples_size());

    std::sort(sizes.begin(), sizes.end(), std::greater<typename Cell::Sample_size_type>());

    LOG4CXX_DEBUG(m_logger, "Best sampled cells histogram:");

    if (top_count > m_cells.size())
        top_count = static_cast<size_t>(m_cells.size());

    double aggregated_fraction = 0.0;

    for (size_t i = 0; i < top_count; ++i)
    {
        double fraction = 100.0 * sizes[i] / sample_point_count;
        aggregated_fraction += fraction;

        LOG4CXX_DEBUG(m_logger, "#" << i + 1 << ": " << sizes[i] << " samples" <<
                      " (fraction: " << std::fixed << std::setprecision(2) << fraction << "%, aggregated fraction: " << aggregated_fraction << "%)");
    }

    aggregated_fraction = 0.0;

    for (size_t i = 0; i < sizes.size(); ++i)
    {
        double fraction = 100.0 * sizes[i] / sample_point_count;
        aggregated_fraction += fraction;

        if (aggregated_fraction >= 90.0)
        {
            LOG4CXX_DEBUG(m_logger, "" << i + 1 << " cells contain " << std::fixed << std::setprecision(2) << aggregated_fraction << "% of all samples");
            break;
        }
    }
}

template<class K, class P, class G>
void Spin_cell_graph_3<K, P, G>::sort_cell_list(Cell_list &cells)
{
    // Sorting algorithm:
    //
    // 1: radix-sort            O(n * l)
    // 2: merge-sort            O(n * log(n) * l)
    // 3: intrusive-sort        O(n * log(n) * l)
    // 3: multiway merge sort   O(n * log(n) * l) + parallel
#if USE_PARALLEL
    const int SORT_ALGORITHM = 4;
#else // USE_PARALLEL
    const int SORT_ALGORITHM = 3;
#endif // USE_PARALLEL

    if (SORT_ALGORITHM == 4)
    {
        // Requires: random-access iterator
        LOG4CXX_DEBUG(m_logger, "Multiway merge sort");

        // parallel sort
//        __gnu_parallel::sort(cells.begin(), cells.end());
        std::sort(cells.begin(), cells.end());
    }
    else if (SORT_ALGORITHM == 3)
    {
        // Requires: random-access iterator
        LOG4CXX_DEBUG(m_logger, "Intrusive sort");

        // intrusive sort
        std::sort(cells.begin(), cells.end());
    }
    else if (SORT_ALGORITHM == 2)
    {
#if 0
        LOG4CXX_DEBUG(m_logger, "Merge sort");

        // merge sort
        cells.sort();
#endif
    }
    else if (SORT_ALGORITHM == 1)
    {
#if 0
        LOG4CXX_DEBUG(m_logger, "Linear sort");

        typedef typename Cell_list::iterator Cell_list_iterator;

        const size_t number_of_cells = cells.size();

        // radix sort on coordinate
        boost::scoped_array<Cell_list_iterator> cell_index(new Cell_list_iterator[number_of_cells]);
        Cell_list_iterator iterator;
        size_t i, j;
        int bit;

        // initialize indexes to identity
        iterator = cells.begin();

        for (i = 0; i < number_of_cells; ++i)
            cell_index[i] = iterator++;

        // radix sort: for each bit, starting from the least significant
        boost::scoped_array<Cell_list_iterator> buffer_zero(new Cell_list_iterator[number_of_cells]);
        boost::scoped_array<Cell_list_iterator> buffer_one(new Cell_list_iterator[number_of_cells]);
        size_t buffer_zero_size, buffer_one_size;

        LOG4CXX_DEBUG(m_logger, "Iterating linear sort");

        for (bit = static_cast<int>(m_coordinate_size) - 1; bit >= 0; --bit)
        {
            // do a stable sort - counting sort
            buffer_zero_size = buffer_one_size = 0;

            // accumulate
            for (i = 0; i < number_of_cells; ++i)
            {
                if (cell_index[i]->coordinate()[bit])
                    buffer_one[buffer_one_size++] = cell_index[i];
                else
                    buffer_zero[buffer_zero_size++] = cell_index[i];
            }

            // propagate
            i = 0;

            for (j = 0; j < buffer_zero_size; ++j)
                cell_index[i++] = buffer_zero[j];

            for (j = 0; j < buffer_one_size; ++j)
                cell_index[i++] = buffer_one[j];
        }

        LOG4CXX_DEBUG(m_logger, "Back copy from linear sort");

        // copy back indexed cells
        for (i = 0; i < number_of_cells; ++i)
            cells.splice(cells.end(), cells, cell_index[i]);
#endif
    }
}

template<class K, class P, class G>
Spin_cell_graph_3<K, P, G>::Route::Route()
    : m_valid(false)
{
}

//template<class K, class P, class G>
//Spin_cell_graph_3<K, P, G>::Route::Route(const std::vector<Voxel_link> &nodes)
//    : m_valid(true),
//      m_nodes(nodes)
//{
//}

template<class K, class P, class G>
bool Spin_cell_graph_3<K, P, G>::Route::is_valid() const
{
    return m_valid;
}

template<class K, class P, class G>
typename Spin_cell_graph_3<K, P, G>::Sample Spin_cell_graph_3<K, P, G>::Route::evaluate(double t) const
{
    (void)t;

//    // for completeness
//    if (m_nodes.empty())
        return Sample();

//    // check corner cases
//    if (m_nodes.size() == 1)
//        return m_nodes.front().first->spinor(m_nodes.front().second);

//    // check bounds
//    if (t < 0.0) t = 0.0;
//    if (t > 1.0) t = 1.0;

//    // find a motion segment
//    int segment = static_cast<int>(t * (m_nodes.size() - 1));

//    if (segment == static_cast<int>(m_nodes.size() - 1))
//        --segment;

//    // use slerp
//    Sample segment_begin = m_nodes[segment].first->spinor(m_nodes[segment].second);
//    Sample segment_end = m_nodes[segment + 1].first->spinor(m_nodes[segment + 1].second);

//    double segment_time = 1.0 / (m_nodes.size() - 1);
//    double segment_start_time = segment * segment_time;

//    double segment_t = (t - segment_start_time) / segment_time;

//    return slerp(segment_begin, segment_end, segment_t);
}
} // namespace CS
