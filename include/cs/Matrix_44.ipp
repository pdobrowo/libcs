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
#include "Matrix_44.h"

namespace CS
{
template<class RT_>
Matrix_44<RT_>::Matrix_44()
{
    set_zero();
}

template<class RT_>
void Matrix_44<RT_>::set_zero()
{
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            m_e[i][j] = 0;
}

template<class RT_>
void Matrix_44<RT_>::set_identity()
{
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            m_e[i][j] = (i == j);
}

template<class RT_>
const RT_ &Matrix_44<RT_>::get(int i, int j) const
{
    return m_e[i][j];
}

template<class RT_>
void Matrix_44<RT_>::set(int i, int j, const RT &v)
{
    m_e[i][j] = v;
}
} // namespace CS
