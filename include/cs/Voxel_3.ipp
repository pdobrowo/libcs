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
#include "Voxel_3.h"

namespace CS
{
template<class Data_>
Voxel_3<Data_>::Voxel_3()
{
    // dummy constructor
}

template<class Data_>
Voxel_3<Data_>::Voxel_3(int u, int v, int w,
                    double px, double py, double pz,
                    bool real,
                    bool negative_value,
                    bool positive_value)
    : m_u(u),
      m_v(v),
      m_w(w),
      m_px(px),
      m_py(py),
      m_pz(pz),
      m_real(real),
      m_border(false)
{
    m_value[Cover_Negative] = negative_value;
    m_value[Cover_Positive] = positive_value;
}

template<class Data_>
int Voxel_3<Data_>::u() const
{
    return m_u;
}

template<class Data_>
int Voxel_3<Data_>::v() const
{
    return m_v;
}

template<class Data_>
int Voxel_3<Data_>::w() const
{
    return m_w;
}

template<class Data_>
double Voxel_3<Data_>::px() const
{
    return m_px;
}

template<class Data_>
double Voxel_3<Data_>::py() const
{
    return m_py;
}

template<class Data_>
double Voxel_3<Data_>::pz() const
{
    return m_pz;
}

template<class Data_>
Spin_3<double> Voxel_3<Data_>::spinor(Cover cover) const
{
    switch (cover)
    {
    case Cover_Negative: return Spin_3<double>(m_px, m_py, m_pz, -std::sqrt(1.0 - m_px * m_px - m_py * m_py - m_pz * m_pz));
    case Cover_Positive: return Spin_3<double>(m_px, m_py, m_pz, std::sqrt(1.0 - m_px * m_px - m_py * m_py - m_pz * m_pz));
    case Cover_COUNT:
    default:
        assert(0);
        return Spin_3<double>();
    }
}

template<class Data_>
bool Voxel_3<Data_>::is_real() const
{
    return m_real;
}

template<class Data_>
bool Voxel_3<Data_>::is_imaginary() const
{
    return !is_real();
}

template<class Data_>
bool Voxel_3<Data_>::is_border() const
{
    return m_border;
}

template<class Data_>
void Voxel_3<Data_>::set_is_border(bool boolean)
{
    m_border = boolean;
}

template<class Data_>
bool Voxel_3<Data_>::value(Cover cover) const
{
    return m_value[cover];
}

template<class Data_>
const typename Voxel_3<Data_>::Data &Voxel_3<Data_>::data() const
{
    return m_data;
}

template<class Data_>
typename Voxel_3<Data_>::Data &Voxel_3<Data_>::data()
{
    return m_data;
}
} // namespace CS
