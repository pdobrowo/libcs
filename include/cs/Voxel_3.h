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
#ifndef LIBCS_VOXEL_3_H
#define LIBCS_VOXEL_3_H

#include "Spin_3.h"
#include <cmath>
#include <cassert>

namespace CS
{
enum Cover
{
    Cover_Negative = 0,
    Cover_Positive = 1,

    Cover_COUNT = 2
};

inline Cover &operator++(Cover &cover)
{
    cover = static_cast<Cover>(static_cast<int>(cover) + 1);
    return cover;
}

inline Cover opposite_cover(Cover cover)
{
    switch (cover)
    {
    case Cover_Negative: return Cover_Positive;
    case Cover_Positive: return Cover_Negative;
    case Cover_COUNT:
    default:
        assert(0);
        return Cover_COUNT;
    }
}

inline Cover which_cover(const Spin_3<double> &spin)
{
    if (spin.s0() >= 0.0)
        return Cover_Positive;
    else
        return Cover_Negative;
}

// 3-dimensional voxel (double covered)
template<class Data_>
class Voxel_3
{
public:
    typedef Data_   Data;

    Voxel_3();

    Voxel_3(int u, int v, int w,
            double px, double py, double pz,
            bool real,
            bool negative_value,
            bool positive_value);

    // raster index
    int         u() const;
    int         v() const;
    int         w() const;

    // absolute position
    double      px() const;
    double      py() const;
    double      pz() const;

    Spin_3<double>  spinor(Cover cover) const;

    // is the voxel in bounds of Spin(3)
    bool        is_real() const;
    bool        is_imaginary() const;

    // border
    bool        is_border() const;
    void        set_is_border(bool boolean);

    // values on the positive and the negative leaf s0 = \pm \sqrt(1 - s_12^2 - s_23^2 - s_31^2)
    bool        value(Cover cover) const;

    const Data &data() const;
    Data &      data();

private:
    // coordinates in raster
    int     m_u, m_v, m_w;

    // absolute position
    double  m_px, m_py, m_pz;

    // is in space
    bool    m_real;
    bool    m_border;

    // voxel value
    bool    m_value[Cover_COUNT];

    // external data
    Data    m_data;
};
} // namespace CS

#include "Voxel_3.ipp"

#endif // LIBCS_VOXEL_3_H
