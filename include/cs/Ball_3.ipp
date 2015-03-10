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
#include "Ball_3.h"

namespace CS
{
template<class Kernel_>
Ball_3<Kernel_>::Ball_3()
    : m_center(Vector_3(RT(0), RT(0), RT(0))),
      m_radius(RT(0))
{
}

template<class Kernel_>
Ball_3<Kernel_>::Ball_3(const RT &x_,
       const RT &y_,
       const RT &z_,
       const RT &radius_)
    : m_center(Vector_3(x_, y_, z_)),
      m_radius(radius_)
{
}

template<class Kernel_>
Ball_3<Kernel_>::Ball_3(const Vector_3 &center_,
       const RT &radius_)
    : m_center(center_),
      m_radius(radius_)
{
}

template<class Kernel_>
const typename Ball_3<Kernel_>::Vector_3 & Ball_3<Kernel_>::center() const
{
    return m_center;
}

template<class Kernel_>
const typename Ball_3<Kernel_>::RT & Ball_3<Kernel_>::radius() const
{
    return m_radius;
}

template<class Kernel_>
void Ball_3<Kernel_>::scale(const RT &scale)
{
    m_center = Vector_3(m_center.x() * scale,
                        m_center.y() * scale,
                        m_center.z() * scale);

    m_radius *= scale;
}
} // namespace CS
