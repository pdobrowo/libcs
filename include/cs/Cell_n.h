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
#ifndef LIBCS_CELL_N_H
#define LIBCS_CELL_N_H

#include "Coordinate.h"
#include <algorithm>
#include <vector>
#include <cassert>

namespace CS
{
struct Empty_cell_data {};

// generic n-dimensional cell
template<class Sample_,
         class Data_ = Empty_cell_data>
class Cell_n
{
public:
    typedef Sample_ Sample;
    typedef Data_   Data;

    // copying Cell_n is forbidden
    // a cell is always pointed by a valid handle, which is
    // an iterator to a NON-MUTABLE list with cells
    typedef typename std::vector<Cell_n>::iterator  Handle;
    typedef std::vector<Handle>                     Handle_list;

private:
    typedef std::vector<Sample>                     Sample_list;

public:
    Cell_n();

    Cell_n(const Sample &sample, const Coordinate &coordinate, bool value);

    const Coordinate &coordinate() const;
    size_t pops_coordinate() const;

    bool value() const;
    bool is_empty() const;
    bool is_full() const;

    const Data &data() const;
    Data &data();

    void add_sample(const Sample &sample);

    typedef typename Sample_list::const_iterator  Sample_const_iterator;
    typedef typename Sample_list::size_type       Sample_size_type;

    Sample_const_iterator samples_begin() const;
    Sample_const_iterator samples_end() const;
    Sample_size_type samples_size() const;

    void add_edge(Handle other);

    typedef typename Handle_list::const_iterator  Edge_const_iterator;
    typedef typename Handle_list::size_type       Edge_size_type;

    Edge_const_iterator edges_begin() const;
    Edge_const_iterator edges_end() const;
    Edge_size_type edges_size() const;

private:
    Coordinate          m_coordinate;
    bool                m_value;
    Handle_list         m_edges;
    Sample_list         m_samples;
    Data                m_data;
    size_t              m_pops;
};

// cell operators are base on coordinate comparisions
template<class Sample, class Data>
bool operator <(const Cell_n<Sample, Data> &lhs, const Cell_n<Sample, Data> &rhs);

template<class Sample, class Data>
bool operator >(const Cell_n<Sample, Data> &lhs, const Cell_n<Sample, Data> &rhs);

template<class Sample, class Data>
bool operator <=(const Cell_n<Sample, Data> &lhs, const Cell_n<Sample, Data> &rhs);

template<class Sample, class Data>
bool operator >=(const Cell_n<Sample, Data> &lhs, const Cell_n<Sample, Data> &rhs);

template<class Sample, class Data>
bool operator ==(const Cell_n<Sample, Data> &lhs, const Cell_n<Sample, Data> &rhs);

template<class Sample, class Data>
bool operator !=(const Cell_n<Sample, Data> &lhs, const Cell_n<Sample, Data> &rhs);

template<class Sample, class Data>
bool neighbour(const Cell_n<Sample, Data> &lhs, const Cell_n<Sample, Data> &rhs);
} // namespace CS

#include "Cell_n.ipp"

#endif // LIBCS_CELL_N_H
