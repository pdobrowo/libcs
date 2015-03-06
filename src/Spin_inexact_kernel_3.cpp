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
#include <cs/Spin_inexact_kernel_3.h>
#include <CGAL/Cartesian.h>

// force instantiation
template<typename RT_>
struct Force_spin_inexact_kernel_3_instantiation
{
    typedef RT_                                     RT;
    typedef CGAL::Cartesian<RT>                     Base_kernel;
    typedef CS::Spin_inexact_kernel_3<Base_kernel>  Kernel;

    typedef typename Kernel::Predicate_tt_3 Predicate_tt_3;
    typedef typename Kernel::Predicate_bb_3 Predicate_bb_3;

    typedef typename Kernel::template Spin_cell_configuration_space_3<Predicate_tt_3>::Type CS_TT_C;
    typedef typename Kernel::template Spin_cell_configuration_space_3<Predicate_bb_3>::Type CS_BB_C;

    typedef typename Kernel::template Spin_raster_configuration_space_3<Predicate_tt_3>::Type CS_TT_R;
    typedef typename Kernel::template Spin_raster_configuration_space_3<Predicate_bb_3>::Type CS_BB_R;
};

typedef Force_spin_inexact_kernel_3_instantiation<float> Force_spin_inexact_kernel_3_instantiation_float;
typedef Force_spin_inexact_kernel_3_instantiation<double> Force_spin_inexact_kernel_3_instantiation_double;
