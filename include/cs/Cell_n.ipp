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
#include "Cell_n.h"

namespace CS
{
template<class Sample_, class Data_>
Cell_n<Sample_, Data_>::Cell_n()
{
}

template<class Sample_, class Data_>
Cell_n<Sample_, Data_>::Cell_n(const Sample &sample, const Coordinate &coordinate, bool value)
    : m_coordinate(coordinate),
      m_value(value)
{
    m_samples.push_back(sample);

    // calculate pops
    m_pops = std::count(m_coordinate.begin(), m_coordinate.end(), true);
}

template<class Sample_, class Data_>
const Coordinate &Cell_n<Sample_, Data_>::coordinate() const
{
    return m_coordinate;
}

template<class Sample_, class Data_>
size_t Cell_n<Sample_, Data_>::pops_coordinate() const
{
    return m_pops;
}

template<class Sample_, class Data_>
bool Cell_n<Sample_, Data_>::value() const
{
    return m_value;
}

template<class Sample_, class Data_>
bool Cell_n<Sample_, Data_>::is_empty() const
{
    return m_value == false;
}

template<class Sample_, class Data_>
bool Cell_n<Sample_, Data_>::is_full() const
{
    return m_value == true;
}

template<class Sample_, class Data_>
const typename Cell_n<Sample_, Data_>::Data &Cell_n<Sample_, Data_>::data() const
{
    return m_data;
}

template<class Sample_, class Data_>
typename Cell_n<Sample_, Data_>::Data &Cell_n<Sample_, Data_>::data()
{
    return m_data;
}

template<class Sample_, class Data_>
void Cell_n<Sample_, Data_>::add_sample(const Sample &sample)
{
    m_samples.push_back(sample);
}

template<class Sample_, class Data_>
typename Cell_n<Sample_, Data_>::Sample_const_iterator Cell_n<Sample_, Data_>::samples_begin() const
{
    return m_samples.begin();
}

template<class Sample_, class Data_>
typename Cell_n<Sample_, Data_>::Sample_const_iterator Cell_n<Sample_, Data_>::samples_end() const
{
    return m_samples.end();
}

template<class Sample_, class Data_>
typename Cell_n<Sample_, Data_>::Sample_size_type Cell_n<Sample_, Data_>::samples_size() const
{
    return m_samples.size();
}

template<class Sample_, class Data_>
void Cell_n<Sample_, Data_>::add_edge(Handle other)
{
    m_edges.push_back(other);
}

template<class Sample_, class Data_>
typename Cell_n<Sample_, Data_>::Edge_const_iterator Cell_n<Sample_, Data_>::edges_begin() const
{
    return m_edges.begin();
}

template<class Sample_, class Data_>
typename Cell_n<Sample_, Data_>::Edge_const_iterator Cell_n<Sample_, Data_>::edges_end() const
{
    return m_edges.end();
}

template<class Sample_, class Data_>
typename Cell_n<Sample_, Data_>::Edge_size_type Cell_n<Sample_, Data_>::edges_size() const
{
    return m_edges.size();
}

template<class Sample_, class Data_>
bool operator <(const Cell_n<Sample_, Data_> &lhs, const Cell_n<Sample_, Data_> &rhs)
{
    return lhs.coordinate() < rhs.coordinate();
}

template<class Sample_, class Data_>
bool operator >(const Cell_n<Sample_, Data_> &lhs, const Cell_n<Sample_, Data_> &rhs)
{
    return lhs.coordinate() > rhs.coordinate();
}

template<class Sample_, class Data_>
bool operator <=(const Cell_n<Sample_, Data_> &lhs, const Cell_n<Sample_, Data_> &rhs)
{
    return lhs.coordinate() <= rhs.coordinate();
}

template<class Sample_, class Data_>
bool operator >=(const Cell_n<Sample_, Data_> &lhs, const Cell_n<Sample_, Data_> &rhs)
{
    return lhs.coordinate() >= rhs.coordinate();
}

template<class Sample_, class Data_>
bool operator ==(const Cell_n<Sample_, Data_> &lhs, const Cell_n<Sample_, Data_> &rhs)
{
    // heuristic based on pops
    if (lhs.pops_coordinate() != rhs.pops_coordinate())
        return false;

    // straight algorithm
    return lhs.coordinate() == rhs.coordinate();
}

template<class Sample_, class Data_>
bool operator !=(const Cell_n<Sample_, Data_> &lhs, const Cell_n<Sample_, Data_> &rhs)
{
    // use equality implementation
    return !(lhs.coordinate() == rhs.coordinate());
}

template<class Sample_, class Data_>
bool neighbour(const Cell_n<Sample_, Data_> &lhs, const Cell_n<Sample_, Data_> &rhs)
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
