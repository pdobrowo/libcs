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
#ifndef LIBCS_MATRIX_44_H
#define LIBCS_MATRIX_44_H

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

    void        set_zero();
    void        set_identity();

    const RT &  get(int i, int j) const;
    void        set(int i, int j, const RT &v);

private:
    RT m_e[4][4];
};
} // namespace CS

#include "Matrix_44.ipp"

#endif // LIBCS_MATRIX_44_H
