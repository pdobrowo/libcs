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
#include "Cell_n.h"

namespace CS
{
template<class S, class D>
Cell_n<S, D>::Cell_n()
{
}

template<class S, class D>
Cell_n<S, D>::Cell_n(const Sample &sample, const Coordinate &coordinate, bool value)
    : m_coordinate(coordinate),
      m_value(value)
{
    m_samples.push_back(sample);

    // calculate pops
    m_pops = std::count(m_coordinate.begin(), m_coordinate.end(), true);
}

template<class S, class D>
const Coordinate &Cell_n<S, D>::coordinate() const
{
    return m_coordinate;
}

template<class S, class D>
size_t Cell_n<S, D>::pops_coordinate() const
{
    return m_pops;
}

template<class S, class D>
bool Cell_n<S, D>::value() const
{
    return m_value;
}

template<class S, class D>
bool Cell_n<S, D>::is_empty() const
{
    return m_value == false;
}

template<class S, class D>
bool Cell_n<S, D>::is_full() const
{
    return m_value == true;
}

template<class S, class D>
const typename Cell_n<S, D>::Data &Cell_n<S, D>::data() const
{
    return m_data;
}

template<class S, class D>
typename Cell_n<S, D>::Data &Cell_n<S, D>::data()
{
    return m_data;
}

template<class S, class D>
void Cell_n<S, D>::add_sample(const Sample &sample)
{
    m_samples.push_back(sample);
}

template<class S, class D>
typename Cell_n<S, D>::Sample_const_iterator Cell_n<S, D>::samples_begin() const
{
    return m_samples.begin();
}

template<class S, class D>
typename Cell_n<S, D>::Sample_const_iterator Cell_n<S, D>::samples_end() const
{
    return m_samples.end();
}

template<class S, class D>
typename Cell_n<S, D>::Sample_size_type Cell_n<S, D>::samples_size() const
{
    return m_samples.size();
}

template<class S, class D>
void Cell_n<S, D>::add_edge(Handle other)
{
    m_edges.push_back(other);
}

template<class S, class D>
typename Cell_n<S, D>::Edge_const_iterator Cell_n<S, D>::edges_begin() const
{
    return m_edges.begin();
}

template<class S, class D>
typename Cell_n<S, D>::Edge_const_iterator Cell_n<S, D>::edges_end() const
{
    return m_edges.end();
}

template<class S, class D>
typename Cell_n<S, D>::Edge_size_type Cell_n<S, D>::edges_size() const
{
    return m_edges.size();
}

template<class S, class D>
bool operator <(const Cell_n<S, D> &lhs, const Cell_n<S, D> &rhs)
{
    return lhs.coordinate() < rhs.coordinate();
}

template<class S, class D>
bool operator >(const Cell_n<S, D> &lhs, const Cell_n<S, D> &rhs)
{
    return lhs.coordinate() > rhs.coordinate();
}

template<class S, class D>
bool operator <=(const Cell_n<S, D> &lhs, const Cell_n<S, D> &rhs)
{
    return lhs.coordinate() <= rhs.coordinate();
}

template<class S, class D>
bool operator >=(const Cell_n<S, D> &lhs, const Cell_n<S, D> &rhs)
{
    return lhs.coordinate() >= rhs.coordinate();
}

template<class S, class D>
bool operator ==(const Cell_n<S, D> &lhs, const Cell_n<S, D> &rhs)
{
    // heuristic based on pops
    if (lhs.pops_coordinate() != rhs.pops_coordinate())
        return false;

    // straight algorithm
    return lhs.coordinate() == rhs.coordinate();
}

template<class S, class D>
bool operator !=(const Cell_n<S, D> &lhs, const Cell_n<S, D> &rhs)
{
    // use equality implementation
    return !(lhs.coordinate() == rhs.coordinate());
}

template<class S, class D>
bool neighbour(const Cell_n<S, D> &lhs, const Cell_n<S, D> &rhs)
{
    assert(lhs.coordinate().size() == rhs.coordinate().size());

    // heuristic: check pops values
    // note that pops values are unsigned
    if (lhs.pops_coordinate() + 1 != rhs.pops_coordinate() &&
        rhs.pops_coordinate() + 1 != lhs.pops_coordinate())
    {
        // these cannot be neighbours for sure
        return false;
    }

    // slow version
    size_t differences = 0;
    Coordinate::size_type coordinateSize = lhs.coordinate().size();

    for (Coordinate::size_type i = 0; i < coordinateSize; ++i)
        if (lhs.coordinate()[i] != rhs.coordinate()[i])
            if (++differences > 1)
                return false;

    return differences == 1;
}
} // namespace CS
