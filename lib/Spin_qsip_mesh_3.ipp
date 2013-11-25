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
#include "Spin_qsip_mesh_3.h"

namespace CS
{
template<class R>
Spin_qsip_mesh_3<R>::Spin_qsip_mesh_3(const Spin_qsip_3 &spin_qsip)
    : m_spin_qsip(spin_qsip)
{
}

template<class R>
size_t Spin_qsip_mesh_3<R>::size_of_points() const
{
    // forward
    return m_spin_qsip.size_of_points();
}

template<class R>
void Spin_qsip_mesh_3<R>::mesh_point(Spin_3 &outSpin, size_t index) const
{
    // refine point up to maximum in Epick
    Spin_qsip_point point = m_spin_qsip.point_at(index);
    Hom_root root = point.root();

    // parametrization
    RT a;
    RT b;

    // check point type
    if (root.is_infinite())
    {
        // parameter value is infinity
        a = 1;
        b = 0;
    }
    else
    {
        static const double QSIC_PARAMETER_REFINE_SCALE = 1000000;

        // parameter value is finite
        Algebraic_real_1 parameter = root.root();
        double parameter_real = CGAL::to_double(parameter);

        a = QSIC_PARAMETER_REFINE_SCALE * parameter_real;
        b = QSIC_PARAMETER_REFINE_SCALE;
    }

    // evaluate component
    Spin_qsic_mesh_3::evaluate_component(outSpin, point.component(), a, b);
}
} // namespace CS
