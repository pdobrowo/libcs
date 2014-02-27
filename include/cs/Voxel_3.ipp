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
template<class D>
Voxel_3<D>::Voxel_3()
{
    // dummy constructor
}

template<class D>
Voxel_3<D>::Voxel_3(int u, int v, int w,
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

template<class D>
int Voxel_3<D>::u() const
{
    return m_u;
}

template<class D>
int Voxel_3<D>::v() const
{
    return m_v;
}

template<class D>
int Voxel_3<D>::w() const
{
    return m_w;
}

template<class D>
double Voxel_3<D>::px() const
{
    return m_px;
}

template<class D>
double Voxel_3<D>::py() const
{
    return m_py;
}

template<class D>
double Voxel_3<D>::pz() const
{
    return m_pz;
}

template<class D>
Spin_3<double> Voxel_3<D>::spinor(Cover cover) const
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

template<class D>
bool Voxel_3<D>::is_real() const
{
    return m_real;
}

template<class D>
bool Voxel_3<D>::is_imaginary() const
{
    return !is_real();
}

template<class D>
bool Voxel_3<D>::is_border() const
{
    return m_border;
}

template<class D>
void Voxel_3<D>::set_is_border(bool boolean)
{
    m_border = boolean;
}

template<class D>
bool Voxel_3<D>::value(Cover cover) const
{
    return m_value[cover];
}

template<class D>
const typename Voxel_3<D>::Data &Voxel_3<D>::data() const
{
    return m_data;
}

template<class D>
typename Voxel_3<D>::Data &Voxel_3<D>::data()
{
    return m_data;
}
} // namespace CS
