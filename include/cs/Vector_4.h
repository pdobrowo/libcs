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
#ifndef LIBCS_VECTOR_4_H
#define LIBCS_VECTOR_4_H

namespace CS
{
// Vector_4:
//
// | a1 |
// | a2 |
// | a3 |
// | a4 |
//
template<class RT_>
class Vector_4
{
public:    
    typedef RT_ RT;

    Vector_4();
    Vector_4(const RT &x, const RT &y, const RT &z, const RT &w);

    void        set_zero();

    const RT &  get(int i) const;
    void        set(int i, const RT &v);

private:
    RT m_e[4];
};
} // namespace CS

#include "Vector_4.ipp"

#endif // LIBCS_VECTOR_4_H
