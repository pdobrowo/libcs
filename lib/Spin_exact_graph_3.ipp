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
#include "Spin_exact_graph_3.h"

namespace CS
{
template<class K, class P>
Spin_exact_graph_3<K, P>::Parameters::Parameters()
    : m_suppress_qsic_calculation(false),
      m_suppress_qsip_calculation(false)
{
}

template<class K, class P>
Spin_exact_graph_3<K, P>::Parameters::Parameters(bool suppress_qsic_calculation, bool suppress_qsip_calculation)
    : m_suppress_qsic_calculation(suppress_qsic_calculation),
      m_suppress_qsip_calculation(suppress_qsip_calculation)
{
}

template<class K, class P>
bool Spin_exact_graph_3<K, P>::Parameters::suppress_qsic_calculation() const
{
    return m_suppress_qsic_calculation;
}

template<class K, class P>
bool Spin_exact_graph_3<K, P>::Parameters::suppress_qsip_calculation() const
{
    return m_suppress_qsip_calculation;
}

template<class K, class P>
Spin_exact_graph_3<K, P>::Spin_exact_graph_3(const std::vector<Predicate> &predicates,
                                             const std::vector<Spin_quadric_3> &spin_quadrics,
                                             const Parameters &parameters)
    : m_spin_quadrics(spin_quadrics),
      m_logger(log4cxx::Logger::getLogger("CS.Spin_exact_graph_3"))
{
    (void)parameters;

    // check if everything is ok with predicates and quadrics
    assert(predicates.size() * SUB_PREDICATE_COUNT == spin_quadrics.size());

    // timing
    unsigned long long now;

    // create qsics
    unsigned long long checkQsicsBegin = get_tick_count();

    if (!parameters.suppress_qsic_calculation())
    {
        m_qsics.reserve(spin_quadrics.size() * spin_quadrics.size());

        for (size_t i = 0; i < spin_quadrics.size(); ++i)
            for (size_t j = i + 1; j < spin_quadrics.size(); ++j)
                m_qsics.push_back(new Spin_qsic_3(spin_quadrics[i], spin_quadrics[j]));
    }

    now = get_tick_count();

    if (!m_qsics.empty())
    {
        LOG4CXX_INFO(m_logger, "Computed " << m_qsics.size() << " spin-QSICs [" << (now - checkQsicsBegin) << " ms] [avg: " << (now - checkQsicsBegin) / m_qsics.size() << " ms]");

        // create qsips
        unsigned long long checkQsipsBegin = get_tick_count();

        if (!parameters.suppress_qsip_calculation())
        {
            m_qsips.reserve(spin_quadrics.size() * m_qsics.size());

            for (size_t i = 0; i < spin_quadrics.size(); ++i)       // quadric i
                for (size_t j = 0; j < m_qsics.size(); ++j)         // qsic j
                    m_qsips.push_back(new Spin_qsip_3(spin_quadrics[i], *m_qsics[j]));
        }

        now = get_tick_count();

        if (!m_qsips.empty())
        {
            LOG4CXX_INFO(m_logger, "Computed " << m_qsips.size() << " spin-QSIPs [" << (now - checkQsipsBegin) << " ms] [avg: " << (now - checkQsipsBegin) / m_qsips.size() << " ms]");
        }
        else
        {
            LOG4CXX_INFO(m_logger, "no qsips created");
        }
    }
    else
    {
        LOG4CXX_INFO(m_logger, "no qsics created");
    }

    LOG4CXX_INFO(m_logger, "Exact representation: " << m_qsics.size() << " spin-QSICs, " << m_qsips.size() << " spin-QSIPs");
}

template<class K, class P>
Spin_exact_graph_3<K, P>::~Spin_exact_graph_3()
{
    release();
}

template<class K, class P>
typename Spin_exact_graph_3<K, P>::Route Spin_exact_graph_3<K, P>::find_route(const Sample &begin, const Sample &end)
{
    (void)begin;
    (void)end;
    return Route();
}

template<class K, class P>
void Spin_exact_graph_3<K, P>::release()
{
    BOOST_FOREACH(Qsic_handle handle, m_qsics)
        delete handle;

    BOOST_FOREACH(Qsip_handle handle, m_qsips)
        delete handle;

    m_qsics.clear();
    m_qsips.clear();
}

template<class K, class P>
typename Spin_exact_graph_3<K, P>::Spin_quadric_const_iterator Spin_exact_graph_3<K, P>::spin_quadrics_begin() const
{
    return m_spin_quadrics.begin();
}

template<class K, class P>
typename Spin_exact_graph_3<K, P>::Spin_quadric_const_iterator Spin_exact_graph_3<K, P>::spin_quadrics_end() const
{
    return m_spin_quadrics.end();
}

template<class K, class P>
typename Spin_exact_graph_3<K, P>::Spin_quadric_size_type Spin_exact_graph_3<K, P>::size_of_spin_quadrics() const
{
    return m_spin_quadrics.size();
}

template<class K, class P>
typename Spin_exact_graph_3<K, P>::Qsic_const_iterator Spin_exact_graph_3<K, P>::qsics_begin() const
{
    return m_qsics.begin();
}

template<class K, class P>
typename Spin_exact_graph_3<K, P>::Qsic_const_iterator Spin_exact_graph_3<K, P>::qsics_end() const
{
    return m_qsics.end();
}

template<class K, class P>
typename Spin_exact_graph_3<K, P>::Qsic_size_type Spin_exact_graph_3<K, P>::size_of_qsics() const
{
    return m_qsics.size();
}

template<class K, class P>
typename Spin_exact_graph_3<K, P>::Qsip_const_iterator Spin_exact_graph_3<K, P>::qsips_begin() const
{
    return m_qsips.begin();
}

template<class K, class P>
typename Spin_exact_graph_3<K, P>::Qsip_const_iterator Spin_exact_graph_3<K, P>::qsips_end() const
{
    return m_qsips.end();
}

template<class K, class P>
typename Spin_exact_graph_3<K, P>::Qsip_size_type Spin_exact_graph_3<K, P>::size_of_qsips() const
{
    return m_qsips.size();
}

template<class K, class P>
Spin_exact_graph_3<K, P>::Route::Route()
    : m_valid(false)
{
}

//template<class K, class P>
//Spin_exact_graph_3<K, P>::Route::Route(const std::vector<Voxel_link> &nodes)
//    : m_valid(true),
//      m_nodes(nodes)
//{
//}

template<class K, class P>
bool Spin_exact_graph_3<K, P>::Route::is_valid() const
{
    return m_valid;
}

template<class K, class P>
typename Spin_exact_graph_3<K, P>::Sample Spin_exact_graph_3<K, P>::Route::evaluate(double t) const
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
