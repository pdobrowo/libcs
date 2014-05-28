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
const typename Matrix_44<RT_>::RT &Matrix_44<RT_>::get(int i, int j) const
{
    return m_e[i][j];
}

template<class RT_>
void Matrix_44<RT_>::set(int i, int j, const RT &v)
{
    m_e[i][j] = v;
}

template<class RT_>
void Matrix_44<RT_>::row_swap(int ra, int rb)
{
    for (int i = 0; i < 4; ++i)
    {
        RT t = m_e[ra][i];
        m_e[ra][i] = m_e[rb][i];
        m_e[rb][i] = t;
    }
}

template<class RT_>
void Matrix_44<RT_>::row_mul(int r, RT s)
{
    for (int i = 0; i < 4; ++i)
        m_e[r][i] *= s;
}

template<class RT_>
void Matrix_44<RT_>::row_mad(int rt, int rs, RT s)
{
    for (int i = 0; i < 4; ++i)
        m_e[rt][i] += m_e[rs][i] * s;
}

template<class RT_>
int Matrix_44<RT_>::rank() const
{
    int r = 0;

    for (int i = 0; i < 4; ++i)
    {
        int z = 0;

        for (int j = 0; j < 4; ++j)
        {
            if (m_e[i][j] == RT(0))
                ++z;
            else
                break;
        }

        if (z < 4)
            ++r;
    }

    return r;
}

template<class RT_>
bool Matrix_44<RT_>::row_pivot(int rc)
{
    if (m_e[rc][rc] == RT(0))
    {
        for (int r = rc + 1; r < 4; ++r)
        {
            if (m_e[r][rc] != RT(0))
            {
                row_swap(rc, r);
                return true;
            }
        }

        return false;
    }

    return true;
}

template<class RT_>
int Matrix_44<RT_>::kernel(Vector_4<RT> out_base[4]) const
{
    Matrix_44<RT> m = *this;

    for (int c = 0; c < 3; ++c)
    {
        // pivot
        if (!m.row_pivot(c))
            continue;

        // mul
        m.row_mul(c, RT(1) / m.get(c, c));

        // mad
        for (int r = c + 1; r < 4; ++r)
            m.row_mad(r, c, -m.get(r, c));
    }

    int r = m.rank();
    int k = 4 - r;

    switch (k)
    {
    case 0:
        break;

    case 1:
        {
            RT ev0[4];
            ev0[3] = RT(1);

            for (int i = 2; i >= 0; --i)
            {
                ev0[i] = RT(0);

                for (int j = 3; j >= i + 1; --j)
                    ev0[i] -= m.get(i, j) * ev0[j];
            }

            out_base[0] = Vector_4<RT>(ev0[0], ev0[1], ev0[2], ev0[3]);
        }
        break;

    case 2:
        {
        }
        break;

    case 3:
        out_base[0] = Vector_4<RT>(RT(1), RT(1), RT(0), RT(0));
        out_base[1] = Vector_4<RT>(RT(1), RT(0), RT(1), RT(0));
        out_base[2] = Vector_4<RT>(RT(1), RT(0), RT(0), RT(1));
        break;

    case 4:
        out_base[0] = Vector_4<RT>(RT(1), RT(0), RT(0), RT(0));
        out_base[1] = Vector_4<RT>(RT(0), RT(1), RT(0), RT(0));
        out_base[2] = Vector_4<RT>(RT(0), RT(0), RT(1), RT(0));
        out_base[3] = Vector_4<RT>(RT(0), RT(0), RT(0), RT(1));
        break;
    }

    return k;
}

} // namespace CS
