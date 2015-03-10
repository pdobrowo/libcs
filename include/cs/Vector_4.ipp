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
#include "Vector_4.h"

namespace CS
{
template<class RT_>
Vector_4<RT_>::Vector_4()
{
    set_zero();
}

template<class RT_>
Vector_4<RT_>::Vector_4(const RT &x, const RT &y, const RT &z, const RT &w)
{
    m_e[0] = x;
    m_e[1] = y;
    m_e[2] = z;
    m_e[3] = w;
}

template<class RT_>
void Vector_4<RT_>::set_zero()
{
    for (int i = 0; i < 4; ++i)
        m_e[i] = 0;
}

template<class RT_>
const RT_ &Vector_4<RT_>::get(int i) const
{
    return m_e[i];
}

template<class RT_>
void Vector_4<RT_>::set(int i, const RT &v)
{
    m_e[i] = v;
}
} // namespace CS
