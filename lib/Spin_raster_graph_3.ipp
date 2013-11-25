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
#include "Spin_raster_graph_3.h"

namespace CS
{
const int DEFAULT_RASTER_DENSITY        = 22;   // 22 ^ 3 = 10 648, 16 degrees
const double MAXIMUM_DIJKSTRA_DISTANCE  = 9.9;  // !NOTE: this is a maximum distance in selected metric between two voxels (natural metric: pi)

// strong connectivity (26 sides)
static const int NUM_STRONG_NEIGHBOURS = 26;

static const int STRONG_NEIGHBOUR_U[26] = { -1,  0,  1, -1,  0,  1, -1,  0,  1,
                                            -1,  0,  1, -1,      1, -1,  0,  1,
                                            -1,  0,  1, -1,  0,  1, -1,  0,  1 };

static const int STRONG_NEIGHBOUR_V[26] = {  1,  1,  1,  0,  0,  0, -1, -1, -1,
                                             1,  1,  1,  0,      0, -1, -1, -1,
                                             1,  1,  1,  0,  0,  0, -1, -1, -1 };

static const int STRONG_NEIGHBOUR_W[26] = { -1, -1, -1, -1, -1, -1, -1, -1, -1,
                                             0,  0,  0,  0,      0,  0,  0,  0,
                                             1,  1,  1,  1,  1,  1,  1,  1,  1 };

// weak connectivity (6 sides)
static const int NUM_WEAK_NEIGHBOURS = 6;

static const int WEAK_NEIGHBOUR_U[6] = { -1,  1,  0,  0,  0,  0 };

static const int WEAK_NEIGHBOUR_V[6] = {  0,  0, -1,  1,  0,  0 };

static const int WEAK_NEIGHBOUR_W[6] = {  0,  0,  0,  0, -1,  1 };

template<class K, class P>
Spin_raster_graph_3<K, P>::Parameters::Parameters()
    : m_resolution(DEFAULT_RASTER_DENSITY)
{
}

template<class K, class P>
Spin_raster_graph_3<K, P>::Parameters::Parameters(size_t density)
    : m_resolution(density)
{
}

template<class K, class P>
size_t Spin_raster_graph_3<K, P>::Parameters::resolution() const
{
    return m_resolution;
}

template<class K, class P>
void Spin_raster_graph_3<K, P>::mark_connected_components_dfs(const Voxel_link &link, int component, bool value)
{
    std::stack<Voxel_link> stack;

    // initialize
    stack.push(link);

    // iterate
    while (!stack.empty())
    {
        // get next node
        Voxel_link current_link = stack.top();
        stack.pop();

        // unpack
        Voxel_iterator  current_iterator = current_link.first;
        Cover           current_cover = current_link.second;

        // check if this voxel has the proper value
        if (current_iterator->value(current_cover) != value)
            continue;

        // check if this voxel been already marked
        if (current_iterator->data().component[current_cover] != -1)
            continue;

        // mark this voxel
        current_iterator->data().component[current_cover] = component;

        // propagate further (raster graph is a 26-graph = cube without central point)
        for (int t = 0; t < NUM_STRONG_NEIGHBOURS; ++t)
        {
            int next_u = current_iterator->u() + STRONG_NEIGHBOUR_U[t];
            int next_v = current_iterator->v() + STRONG_NEIGHBOUR_V[t];
            int next_w = current_iterator->w() + STRONG_NEIGHBOUR_W[t];

            // add voxel if it is in the space
            if (next_u >= 0 && next_u < static_cast<int>(m_resolution) &&
                next_v >= 0 && next_v < static_cast<int>(m_resolution) &&
                next_w >= 0 && next_w < static_cast<int>(m_resolution))
            {
                Voxel_iterator next_iterator = &voxel(next_u, next_v, next_w);

                // add voxel if it is in the space
                if (next_iterator->is_real())
                    stack.push(Voxel_link(next_iterator, current_cover));
            }
        }

        // propagate to a different cover if applicable
        if (current_iterator->is_border())
        {
            // add voxel on opposite cover
            stack.push(Voxel_link(current_iterator, opposite_cover(current_cover)));
        }
    }
}

template<class K, class P>
std::pair<size_t, size_t> Spin_raster_graph_3<K, P>::mark_connected_components(Voxel_iterator begin, Voxel_iterator end)
{
    // clear all components coordinates
    for (Voxel_iterator iterator = begin; iterator != end; ++iterator)
    {
        iterator->data().component[Cover_Negative] = -1;
        iterator->data().component[Cover_Positive] = -1;
    }

    // mark empty components
    std::pair<size_t, size_t> num_connected_components;
    int component = 0;

    for (Voxel_iterator iterator = begin; iterator != end; ++iterator)
    {
        // mark empty negative components
        for (Cover cover = Cover_Negative; cover <= Cover_Positive; ++cover)
        {
            if (iterator->is_real() &&
                !iterator->value(cover) && // false
                iterator->data().component[cover] == -1)
            {
                mark_connected_components_dfs(Voxel_link(iterator, cover), component, false);
                ++component;
            }
        }
    }

    // save empty component count
    num_connected_components.first = component;

    // mark full components
    for (Voxel_iterator iterator = begin; iterator != end; ++iterator)
    {
        // mark empty negative components
        for (Cover cover = Cover_Negative; cover <= Cover_Positive; ++cover)
        {
            if (iterator->is_real() &&
                iterator->value(cover) && // true
                iterator->data().component[cover] == -1)
            {
                mark_connected_components_dfs(Voxel_link(iterator, cover), component, true);
                ++component;
            }
        }
    }

    // save full component count
    num_connected_components.second = component - num_connected_components.first;

    return num_connected_components;
}

template<class K, class P>
Spin_raster_graph_3<K, P>::Spin_raster_graph_3(const std::vector<Predicate> &predicates,
                                               const std::vector<Spin_quadric_3> &quadrics,
                                               const Parameters &parameters)
    : m_number_of_empty_voxels(0),
      m_number_of_full_voxels(0),
      m_number_of_mixed_voxels(0),
      m_number_of_voxels(0),
      m_number_of_real_voxels(0),
      m_number_of_border_voxels(0),
      m_logger(log4cxx::Logger::getLogger("CS.Spin_raster_graph_3"))
{
    // check if everything is ok with predicates and quadrics
    assert(predicates.size() * SUB_PREDICATE_COUNT == quadrics.size());

    // read parameters
    m_resolution = parameters.resolution();
    m_number_of_voxels = m_resolution * m_resolution * m_resolution;

    // ok to begin construction
    LOG4CXX_DEBUG(m_logger, "Creating spin raster graph");
    LOG4CXX_DEBUG(m_logger, "Sampling " << m_number_of_voxels << " voxels (density is: " << m_resolution << ")");

    collect_samples(predicates, quadrics);

    LOG4CXX_INFO(m_logger, "Final structure: ");
    LOG4CXX_INFO(m_logger, "    empty voxel(s): "   << m_number_of_empty_voxels);
    LOG4CXX_INFO(m_logger, "    full voxel(s): "    << m_number_of_full_voxels);
    LOG4CXX_INFO(m_logger, "    mixed voxel(s): "   << m_number_of_mixed_voxels);
    LOG4CXX_INFO(m_logger, "    voxel(s): "         << m_number_of_voxels);
    LOG4CXX_INFO(m_logger, "    real voxel(s): "    << m_number_of_real_voxels);
    LOG4CXX_INFO(m_logger, "    border voxel(s): "  << m_number_of_border_voxels);

    // collect components
    LOG4CXX_INFO(m_logger, "Collecting connected components information");

    std::pair<int, int> num_connected_components = mark_connected_components(m_voxels.get(), m_voxels.get() + m_number_of_voxels);

    LOG4CXX_INFO(m_logger, "Number of empty connected components: " << num_connected_components.first);
    LOG4CXX_INFO(m_logger, "Number of full connected components: " << num_connected_components.second);
    LOG4CXX_INFO(m_logger, "Total number of connected components: " << (num_connected_components.first + num_connected_components.second));
}

template<class K, class P>
typename Spin_raster_graph_3<K, P>::Route Spin_raster_graph_3<K, P>::find_route(const Sample &begin, const Sample &end)
{
    Voxel &begin_voxel = locate_voxel(begin);
    Voxel &end_voxel = locate_voxel(end);

    Cover begin_cover = which_cover(begin);
    Cover end_cover = which_cover(end);

    // initial checks
    if (begin_voxel.value(begin_cover)) // check for obstruction at begin or end rotation
    {
        LOG4CXX_WARN(m_logger, "Route not found: begin rotation is obstructed");
        return Route();
    }

    if (end_voxel.value(end_cover)) // check for obstruction at begin or end rotation
    {
        LOG4CXX_WARN(m_logger, "Route not found: end rotation is obstructed");
        return Route();
    }

    if (begin_voxel.data().component[begin_cover] != end_voxel.data().component[end_cover]) // both begin and end voxels must be in the same component
    {
        LOG4CXX_WARN(m_logger, "Route not found: different components for begin and end rotation");
        return Route();
    }

    std::vector<Voxel_link> nodes;

#if 1 // Dijkstra version

    // locate the shortest path in the component by using Dijkstra
    mark_routing_distances_dijkstra(Voxel_link(&begin_voxel, begin_cover));

#else // BFS version

    // locate the shortest path in the component by using bfs
    mark_routing_distances_bfs(Voxel_link(&begin_voxel, begin_cover));

#endif

    // backtrace a route
    collect_backtrace_route(Voxel_link(&end_voxel, end_cover), std::back_inserter(nodes));

    LOG4CXX_INFO(m_logger, "Found a route with " << nodes.size() << " node(s)");

    return Route(nodes);
}

#if 0
template<class K, class P>
void Spin_raster_graph_3<K, P>::mark_routing_distances_bfs(Voxel_link source)
{
    LOG4CXX_INFO(m_logger, "BFS routing");

    // unpack link
    Voxel_iterator source_iterator = source.first;
    Cover source_cover = source.second;

    // clear all components coordinates
    for (Voxel_iterator it = m_voxels.get(); it != m_voxels.get() + m_number_of_voxels; ++it)
    {
        for (Cover cover = Cover_Negative; cover <= Cover_Positive; ++cover)
        {
            it->data().common_color[cover] = COLOR_WHITE;
            it->data().common_distance[cover] = 9.9;                           // !NOTE: this is a maximum distance in selected metric between two voxels (natural metric: pi)
            it->data().common_parent[cover] = Voxel_link(0, Cover_COUNT);
        }
    }

    std::deque<Voxel_link> queue;

    // initialize
    source_iterator->data().common_color[source_cover] = COLOR_GRAY;
    source_iterator->data().common_distance[source_cover] = 0;
    source_iterator->data().common_parent[source_cover] = Voxel_link(0, Cover_COUNT);
    queue.push_back(source);

    // BFS
    while (!queue.empty())
    {
        // get next node
        Voxel_link current_link = queue.front();
        queue.pop_front();

        // unpack
        Voxel_iterator  current_iterator = current_link.first;
        Cover           current_cover = current_link.second;

        // propagate further (raster graph is a 26-graph = cube without central point)
        for (int t = 0; t < NUM_STRONG_NEIGHBOURS; ++t)
        {
            int next_u = current_iterator->u() + STRONG_NEIGHBOUR_U[t];
            int next_v = current_iterator->v() + STRONG_NEIGHBOUR_V[t];
            int next_w = current_iterator->w() + STRONG_NEIGHBOUR_W[t];

            // check bounds
            if (next_u < 0 || next_u >= static_cast<int>(m_resolution) ||
                next_v < 0 || next_v >= static_cast<int>(m_resolution) ||
                next_w < 0 || next_w >= static_cast<int>(m_resolution))
            {
                continue;
            }

            Voxel_iterator next_iterator = &voxel(next_u, next_v, next_w);

            // check component
            if (next_iterator->data().component[current_cover] != source_iterator->data().component[source_cover])
                continue;

            // check bfs condition
            if (next_iterator->data().common_color[current_cover] == COLOR_WHITE)
            {
                next_iterator->data().common_color[current_cover] = COLOR_GRAY;
                next_iterator->data().common_distance[current_cover] = current_iterator->data().common_distance[current_cover] + 1;
                next_iterator->data().common_parent[current_cover] = current_link;

                queue.push_back(Voxel_link(next_iterator, current_cover));
            }
        }

        // propagate to a different cover if applicable
        if (current_iterator->is_border())
        {
            Cover opposite_current_cover = opposite_cover(current_cover);

            // check component
            if (current_iterator->data().component[opposite_current_cover] == source_iterator->data().component[source_cover])
            {
                // check bfs condition
                if (current_iterator->data().common_color[opposite_current_cover] == COLOR_WHITE)
                {
                    current_iterator->data().common_color[opposite_current_cover] = COLOR_GRAY;
                    current_iterator->data().common_distance[opposite_current_cover] = current_iterator->data().common_distance[current_cover] + 1;
                    current_iterator->data().common_parent[opposite_current_cover] = current_link;

                    queue.push_back(Voxel_link(current_iterator, opposite_current_cover));
                }
            }
        }

        // mark current node as black
        current_iterator->data().common_color[current_cover] = COLOR_BLACK;
    }
}
#endif

template<class K, class P>
template<typename OutputIterator>
void Spin_raster_graph_3<K, P>::collect_backtrace_route(Voxel_link target, OutputIterator out)
{
    std::stack<Voxel_link> stack;

    Voxel_link current = target;

    // bactrace
    do
    {
        stack.push(current);
        current = current.first->data().common_parent[current.second];
    }
    while (current.first && current.first->data().common_distance[current.second] > 0);

    // reverse order
    while (!stack.empty())
    {
        *out++ = stack.top();
        stack.pop();
    }
}

template<class K, class P>
size_t Spin_raster_graph_3<K, P>::resolution() const
{
    return m_resolution;
}

template<class K, class P>
const typename Spin_raster_graph_3<K, P>::Voxel &Spin_raster_graph_3<K, P>::voxel(int u, int v, int w) const
{
    return m_voxels[m_resolution * (m_resolution * w + v) + u];
}

template<class K, class P>
typename Spin_raster_graph_3<K, P>::Voxel &Spin_raster_graph_3<K, P>::voxel(int u, int v, int w)
{
    return m_voxels[m_resolution * (m_resolution * w + v) + u];
}

template<class K, class P>
typename Spin_raster_graph_3<K, P>::Voxel &Spin_raster_graph_3<K, P>::locate_voxel(const Sample &sample)
{
    // a sample is always in bounds; no need to check
    // the formula is: s12 = 2.0 * double(u) / double(m_resolution - 1) - 1.0;
    int u = static_cast<int>(round((sample.s12() + 1.0) * double(m_resolution - 1) / 2.0));
    int v = static_cast<int>(round((sample.s23() + 1.0) * double(m_resolution - 1) / 2.0));
    int w = static_cast<int>(round((sample.s31() + 1.0) * double(m_resolution - 1) / 2.0));

    assert(u >= 0 && u < static_cast<int>(m_resolution));
    assert(v >= 0 && v < static_cast<int>(m_resolution));
    assert(w >= 0 && w < static_cast<int>(m_resolution));

    Voxel &result = voxel(u, v, w);

    assert(result.is_real());

    return result;
}

template<class K, class P>
void Spin_raster_graph_3<K, P>::collect_samples(const std::vector<Predicate> &predicates,
                                                const std::vector<Spin_quadric_3> &quadrics)
{
    LOG4CXX_DEBUG(m_logger, "Single-threaded collect");

    // initialize buffer
    m_voxels.reset(new Voxel[m_number_of_voxels]);

    // collect samples
    Sample positive_sample;
    Sample negative_sample;
    double s12, s23, s31, abs_s0, partial_norm;

    Coordinate positive_coordinate;
    positive_coordinate.resize(quadrics.size());

    Coordinate negative_coordinate;
    negative_coordinate.resize(quadrics.size());

    for (size_t u = 0; u < m_resolution; ++u)
    {
        for (size_t v = 0; v < m_resolution; ++v)
        {
            for (size_t w = 0; w < m_resolution; ++w)
            {
                // calculate random sample
                s12 = 2.0 * double(u) / double(m_resolution - 1) - 1.0;
                s23 = 2.0 * double(v) / double(m_resolution - 1) - 1.0;
                s31 = 2.0 * double(w) / double(m_resolution - 1) - 1.0;

                partial_norm = s12 * s12 + s23 * s23 + s31 * s31;

                if (partial_norm > 1.0)
                {
                    // this voxel is out of bounds; not real
                    voxel(u, v, w) = Voxel(u, v, w, s12, s23, s31, false, false, false);
                    continue;
                }

                // this voxel is in Spin(3)
                ++m_number_of_real_voxels;

                // create sample
                abs_s0 = std::sqrt(1.0 - partial_norm);
                positive_sample = Sample(s12, s23, s31, abs_s0);
                negative_sample = Sample(s12, s23, s31, -abs_s0);

                // calculate sample coodrinates
                calculate_raster_coordinate(quadrics, positive_sample, positive_coordinate);
                calculate_raster_coordinate(quadrics, negative_sample, negative_coordinate);

                // evaluate predicate value at index
                bool positive_value = evaluate_predicate_list_at_raster_coordinate(predicates, positive_coordinate);
                bool negative_value = evaluate_predicate_list_at_raster_coordinate(predicates, negative_coordinate);

                // setup a real voxel
                voxel(u, v, w) = Voxel(u, v, w, s12, s23, s31, true, positive_value, negative_value);

                // collect statistics
                if (positive_value && negative_value)
                    ++m_number_of_full_voxels;
                else if (!positive_value && !negative_value)
                    ++m_number_of_empty_voxels;
                else
                    ++m_number_of_mixed_voxels;
            }
        }
    }

    LOG4CXX_DEBUG(m_logger, "Collected " << m_number_of_real_voxels << " real voxel(s)");

    // find border (not optimized)
    // border are those voxels that have at least one empty neighbour in standard direction
    for (size_t u = 0; u < m_resolution; ++u)
    {
        for (size_t v = 0; v < m_resolution; ++v)
        {
            for (size_t w = 0; w < m_resolution; ++w)
            {
                // check only real voxels
                Voxel *current = &voxel(u, v, w);

                if (current->is_real())
                {
                    // check all directions
                    for (int t = 0; t < NUM_WEAK_NEIGHBOURS; ++t)
                    {
                        int nu = u + WEAK_NEIGHBOUR_U[t];
                        int nv = v + WEAK_NEIGHBOUR_V[t];
                        int nw = w + WEAK_NEIGHBOUR_W[t];

                        // is next voxel imaginary
                        if (nu >= 0 && nu < static_cast<int>(m_resolution) &&
                            nv >= 0 && nv < static_cast<int>(m_resolution) &&
                            nw >= 0 && nw < static_cast<int>(m_resolution))
                        {
                            Voxel *next = &voxel(nu, nv, nw);

                            if (next->is_imaginary())
                            {
                                current->set_is_border(true);
                                ++m_number_of_border_voxels;
                                break;
                            }
                        }
                    }
                }
            }
        }
    }

    LOG4CXX_DEBUG(m_logger, "Collected " << m_number_of_border_voxels << " border voxel(s)");
}

template<class K, class P>
void Spin_raster_graph_3<K, P>::calculate_raster_coordinate(
        const std::vector<Spin_quadric_3> &quadrics,
        const Sample &sample,
        Coordinate &out) const
{
    assert(out.size() == quadrics.size());

    // FIXME: handle bad hits of quadrics
    CGAL::Sign sign;

    for (size_t i = 0; i < quadrics.size(); ++i)
    {
        sign = quadrics[i].evaluate_sign(sample);
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

template<class K, class P>
bool Spin_raster_graph_3<K, P>::evaluate_predicate_list_at_raster_coordinate(
        const std::vector<Predicate> &predicates,
        const Coordinate &coordinate) const
{
    // predicate collision sentence is predicates[0](spin) || predicates[1](spin) || ... || predicates[N](spin)
    // evaluate until result is known
    for (size_t i = 0; i < predicates.size(); ++i)
    {
        Predicate_index index;

        for (size_t c = 0; c < SUB_PREDICATE_COUNT; ++c)
            index[c] = coordinate[SUB_PREDICATE_COUNT * i + c];

        bool value = predicates[i].evaluate(index);

        if (value)
            return true; // collision
    }

    // no collision
    return false;
}

template<class K, class P>
Spin_raster_graph_3<K, P>::Route::Route()
    : m_valid(false),
      m_logger(log4cxx::Logger::getLogger("CS.Spin_raster_graph_3::Route"))
{
}

template<class K, class P>
Spin_raster_graph_3<K, P>::Route::Route(const std::vector<Voxel_link> &nodes)
    : m_valid(true),
      m_nodes(nodes),
      m_logger(log4cxx::Logger::getLogger("CS.Spin_raster_graph_3::Route"))
{
#if 0
    LOG4CXX_DEBUG(m_logger, "Route nodes (" << m_nodes.size() << "):");

    std::vector<Sample> sample_nodes = this->nodes();
    size_t index = 0;

    for (std::vector<Sample>::const_iterator it = sample_nodes.begin(); it != sample_nodes.end(); ++it)
        LOG4CXX_DEBUG(m_logger, "    Node[" << index++ << "] = " << *it);

#endif
}

template<class K, class P>
bool Spin_raster_graph_3<K, P>::Route::is_valid() const
{
    return m_valid;
}

template<class K, class P>
typename Spin_raster_graph_3<K, P>::Sample Spin_raster_graph_3<K, P>::Route::evaluate(double t) const
{
    // for completeness
    if (m_nodes.empty())
        return Sample();

    // check corner cases
    if (m_nodes.size() == 1)
        return m_nodes.front().first->spinor(m_nodes.front().second);

    // check bounds
    if (t < 0.0) t = 0.0;
    if (t > 1.0) t = 1.0;

    // find a motion segment
    int segment = static_cast<int>(t * (m_nodes.size() - 1));

    if (segment == static_cast<int>(m_nodes.size() - 1))
        --segment;

    // use slerp
    Sample segment_begin = m_nodes[segment].first->spinor(m_nodes[segment].second);
    Sample segment_end = m_nodes[segment + 1].first->spinor(m_nodes[segment + 1].second);

    double segment_time = 1.0 / (m_nodes.size() - 1);
    double segment_start_time = segment * segment_time;

    double segment_t = (t - segment_start_time) / segment_time;

    return slerp(segment_begin, segment_end, segment_t);
}

template<class K, class P>
std::vector<typename Spin_raster_graph_3<K, P>::Sample> Spin_raster_graph_3<K, P>::Route::nodes() const
{
    std::vector<Sample> nodes;

    for (typename std::vector<Voxel_link>::const_iterator it = m_nodes.begin(); it != m_nodes.end(); ++it)
        nodes.push_back(it->first->spinor(it->second));

    return nodes;
}

template<class K, class P>
void Spin_raster_graph_3<K, P>::mark_routing_distances_dijkstra(Voxel_link source)
{
    LOG4CXX_INFO(m_logger, "Dijkstra routing");

    // unpack link
    Voxel_iterator source_iterator = source.first;
    Cover source_cover = source.second;

    // clear all components coordinates
    for (Voxel_iterator it = m_voxels.get(); it != m_voxels.get() + m_number_of_voxels; ++it)
    {
        for (Cover cover = Cover_Negative; cover <= Cover_Positive; ++cover)
        {
            it->data().common_color[cover] = COLOR_WHITE;
            it->data().common_distance[cover] = MAXIMUM_DIJKSTRA_DISTANCE;
            it->data().dijkstra_heap_node[cover] = Heap_element_handle();
            it->data().common_parent[cover] = Voxel_link(0, Cover_COUNT);
        }
    }

    // Dijkstra on a Fibonacci heap
    Heap heap;

    // initialize
    source_iterator->data().common_color[source_cover] = COLOR_GRAY;
    source_iterator->data().common_distance[source_cover] = 0.0;
    source_iterator->data().common_parent[source_cover] = Voxel_link(0, Cover_COUNT);
    Heap_element_handle heap_node = heap.push(Heap_element(source, 0.0));
    source_iterator->data().dijkstra_heap_node[source_cover] = heap_node;

    // Dijkstra
    while (!heap.empty())
    {
        Heap_element top = heap.top();
        heap.pop();

        Voxel_link current_link = top.link;

        // unpack
        Voxel_iterator  current_iterator = current_link.first;
        Cover           current_cover = current_link.second;

        // propagate further (raster graph is a 26-graph = cube without central point)
        for (int t = 0; t < NUM_STRONG_NEIGHBOURS; ++t)
        {
            int next_u = current_iterator->u() + STRONG_NEIGHBOUR_U[t];
            int next_v = current_iterator->v() + STRONG_NEIGHBOUR_V[t];
            int next_w = current_iterator->w() + STRONG_NEIGHBOUR_W[t];

            // check bounds
            if (next_u < 0 || next_u >= static_cast<int>(m_resolution) ||
                next_v < 0 || next_v >= static_cast<int>(m_resolution) ||
                next_w < 0 || next_w >= static_cast<int>(m_resolution))
            {
                continue;
            }

            Voxel_iterator next_iterator = &voxel(next_u, next_v, next_w);

            // check component
            if (next_iterator->data().component[current_cover] != source_iterator->data().component[source_cover])
                continue;

            Voxel_link next_link = Voxel_link(next_iterator, current_cover);

            // calculate the distance to this neighbour
            double distance = natural_voxel_distance(current_link, next_link);
            double updated_distance = current_iterator->data().common_distance[current_cover] + distance;

            assert(updated_distance < MAXIMUM_DIJKSTRA_DISTANCE);

            // check Dijkstra condition
            if (next_iterator->data().common_distance[current_cover] > updated_distance)
            {
                next_iterator->data().common_distance[current_cover] = updated_distance;
                next_iterator->data().common_parent[current_cover] = current_link;

                if (next_iterator->data().common_color[current_cover] == COLOR_WHITE)
                {
                    next_iterator->data().common_color[current_cover] = COLOR_GRAY;
                    heap_node = heap.push(Heap_element(next_link, updated_distance));
                    next_iterator->data().dijkstra_heap_node[current_cover] = heap_node;
                }
                else if (next_iterator->data().common_color[current_cover] == COLOR_GRAY)
                {
                    heap_node = next_iterator->data().dijkstra_heap_node[current_cover];
                    heap.decrease(heap_node, Heap_element(next_link, updated_distance));
                }
            }
        }

        // propagate to a different cover if applicable
        if (current_iterator->is_border())
        {
            Cover opposite_current_cover = opposite_cover(current_cover);

            // check component
            if (current_iterator->data().component[opposite_current_cover] == source_iterator->data().component[source_cover])
            {
                Voxel_link next_link = Voxel_link(current_iterator, opposite_current_cover);

                // calculate the distance to this neighbour
                double distance = natural_voxel_distance(current_link, next_link);
                double updated_distance = current_iterator->data().common_distance[current_cover] + distance;

                // check bfs condition
                if (current_iterator->data().common_distance[opposite_current_cover] > updated_distance)
                {
                    current_iterator->data().common_distance[opposite_current_cover] = updated_distance;
                    current_iterator->data().common_parent[opposite_current_cover] = current_link;

                    if (current_iterator->data().common_color[opposite_current_cover] == COLOR_WHITE)
                    {
                        current_iterator->data().common_color[opposite_current_cover] = COLOR_GRAY;
                        heap_node = heap.push(Heap_element(next_link, updated_distance));
                        current_iterator->data().dijkstra_heap_node[opposite_current_cover] = heap_node;
                    }
                    else if (current_iterator->data().common_color[opposite_current_cover] == COLOR_GRAY)
                    {
                        heap_node = current_iterator->data().dijkstra_heap_node[opposite_current_cover];
                        heap.decrease(heap_node, Heap_element(next_link, updated_distance));
                    }
                }
            }
        }

        // mark current node as black
        current_iterator->data().common_color[current_cover] = COLOR_BLACK;
    }
}

template<class K, class P>
double Spin_raster_graph_3<K, P>::natural_voxel_distance(Voxel_link begin, Voxel_link end)
{
    // natural metric on S^3
    Spin_3<double> begin_spin = begin.first->spinor(begin.second);
    Spin_3<double> end_spin = end.first->spinor(end.second);

    return 2.0 * std::acos(begin_spin.s12() * end_spin.s12() + begin_spin.s23() * end_spin.s23() +
                           begin_spin.s31() * end_spin.s31() + begin_spin.s0() * end_spin.s0());

}

template<class K, class P>
double Spin_raster_graph_3<K, P>::cartesian_voxel_distance(Voxel_link begin, Voxel_link end)
{
    // ambient Cartesian metric R^4
    Spin_3<double> begin_spin = begin.first->spinor(begin.second);
    Spin_3<double> end_spin = end.first->spinor(end.second);

    return std::sqrt((begin_spin.s12() - end_spin.s12()) * (begin_spin.s12() - end_spin.s12()) +
                     (begin_spin.s23() - end_spin.s23()) * (begin_spin.s23() - end_spin.s23()) +
                     (begin_spin.s31() - end_spin.s31()) * (begin_spin.s31() - end_spin.s31()) +
                     (begin_spin.s0()  - end_spin.s0())  * (begin_spin.s0()  - end_spin.s0()));
}
} // namespace CS
