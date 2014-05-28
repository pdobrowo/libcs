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
#ifndef LIBCS_SPIN_QSIP_MESH_3
#define LIBCS_SPIN_QSIP_MESH_3

#include <CGAL/Kernel/global_functions.h>
#include <cstddef>

namespace CS
{
// Spin_qsip_mesh_3:
//     mesher
//
// This is a decorator. It does not take ownership of the object.
template<class Kernel_>
class Spin_qsip_mesh_3
{
    typedef typename Kernel_::RT                                         RT;

public:
    typedef typename Kernel_::Spin_qsip_3                                Spin_qsip_3;
    typedef typename Kernel_::Spin_qsic_mesh_3                           Spin_qsic_mesh_3;
    typedef typename Kernel_::Algebraic_kernel_with_sqrt                 Algebraic_kernel_with_sqrt;
    typedef typename Algebraic_kernel_with_sqrt::Algebraic_real_1       Algebraic_real_1;

    typedef typename Spin_qsic_mesh_3::Spin_3                           Spin_3;

    typedef typename Kernel_::Spin_qsip_point                            Spin_qsip_point;
    typedef typename Kernel_::Hom_root                                   Hom_root;

    typedef Kernel_                                                      Kernel;

public:
    Spin_qsip_mesh_3(const Spin_qsip_3 &spin_qsip);

    size_t  size_of_points() const;
    void    mesh_point(Spin_3 &outPoint, size_t index) const;

private:
    const Spin_qsip_3       &m_spin_qsip;
};
} // namespace CS

#include "Spin_qsip_mesh_3.ipp"

#endif // LIBCS_SPIN_QSIP_MESH_3
