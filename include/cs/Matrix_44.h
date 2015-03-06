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
#ifndef LIBCS_MATRIX_44_H
#define LIBCS_MATRIX_44_H

#include "Vector_4.h"
#include <ostream>

namespace CS
{
// Matrix_44:
//
// | a11 a12 a13 a14 |
// | a12 a22 a23 a24 |
// | a13 a23 a33 a34 |
// | a14 a24 a34 a44 |
//
template<class RT_>
class Matrix_44
{
public:    
    typedef RT_ RT;

    Matrix_44();

    template<class RT2_>
    Matrix_44(const Matrix_44<RT2_> &other);

    static Matrix_44    diagonal(const RT &s);

    void        set_zero();
    void        set_identity();

    const RT &  get(int i, int j) const;
    void        set(int i, int j, const RT &v);

    int         rank() const;

    void        row_swap(int ra, int rb);
    void        row_mul(int r, const RT &s);
    void        row_mad(int rt, int rs, const RT &s);
    bool        row_pivot(int rc);

    int         kernel(Vector_4<RT> out_base[4]) const;

private:
    RT m_e[4][4];
};

template<class RT_>
std::ostream &operator <<(std::ostream &os, const Matrix_44<RT_> &matrix);
} // namespace CS

#include "Matrix_44.ipp"

#endif // LIBCS_MATRIX_44_H
