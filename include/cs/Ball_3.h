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
#ifndef LIBCS_BALL_3_H
#define LIBCS_BALL_3_H

namespace CS
{
template<class Kernel_>
class Ball_3
{
    typedef typename Kernel_::RT        RT;
    typedef typename Kernel_::Vector_3  Vector_3;

public:
    typedef Kernel_                     Kernel;

    Ball_3();

    Ball_3(const RT &x_,
           const RT &y_,
           const RT &z_,
           const RT &radius_);

    Ball_3(const Vector_3 &center_,
           const RT &radius_);

    const Vector_3 &    center() const;
    const RT &          radius() const;

    void                scale(const RT &scale);

private:
    Vector_3    m_center;
    RT          m_radius;
};
} // namespace CS

#include "Ball_3.ipp"

#endif // LIBCS_BALL_3_H
