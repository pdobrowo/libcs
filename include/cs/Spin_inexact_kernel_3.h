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
#ifndef LIBCS_SPIN_INEXACT_KERNEL_H
#define LIBCS_SPIN_INEXACT_KERNEL_H

#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Polyhedron_3.h>

// math utils
#define LIBCS_SPIN_INEXACT_KERNEL_INCLUDE_FILE
#include "Math_utils.h"
#undef LIBCS_SPIN_INEXACT_KERNEL_INCLUDE_FILE

// predicates
#include "Predicate_h_3.h"
#include "Predicate_s_3.h"
#include "Predicate_g_3.h"
#include "Predicate_tt_3.h"
#include "Predicate_bb_3.h"
#include "Ball_3.h"

// spin
#include "Spin_quadric_3.h"
#include "Spin_reduced_quadric_3.h"

// meshers
#include "Spin_quadric_mesh_3.h"

// sentence
#include "Spin_quadric_tree_3.h"

// configuration space
#include "Spin_null_graph_3.h"
#include "Spin_cell_graph_3.h"
#include "Spin_raster_graph_3.h"
#include "Spin_configuration_space_3.h"

// matrix
#include "Matrix_44.h"

namespace CS
{
// Spin_inexact_kernel_3 is an adapter for a standard CGAL geometric kernel
// this is an inexact and reduced version
// several things are not available in this version - QSICs, QSIPs and similar
template<class K_>
struct Spin_inexact_kernel_3
    : public K_
{
    // derived types
    typedef typename K_::RT                         RT;

    // spin kernel
    typedef Spin_inexact_kernel_3                   Kernel;

    // spin types
    typedef CS::Spin_quadric_3<Kernel>              Spin_quadric_3;
    typedef CS::Spin_reduced_quadric_3<Kernel>      Spin_reduced_quadric_3;

    // predicates
    typedef CS::Predicate_h_3<Kernel>               Predicate_h_3;
    typedef CS::Predicate_s_3<Kernel>               Predicate_s_3;
    typedef CS::Predicate_g_3<Kernel>               Predicate_g_3;
    typedef CS::Predicate_tt_3<Kernel>              Predicate_tt_3;
    typedef CS::Predicate_bb_3<Kernel>              Predicate_bb_3;

    typedef CS::Ball_3<Kernel>                      Ball_3;
    typedef CGAL::Polyhedron_3<Kernel>              Polyhedron_3;

    // meshers
    typedef CS::Spin_quadric_mesh_3<Kernel>         Spin_quadric_mesh_3;

    // inexact spin
    typedef CS::Spin_3<RT>                          Spin_3;
    typedef CS::Diff_spin_3<RT>                     Diff_spin_3;

    // sentence
    typedef CS::Spin_quadric_tree_3<Kernel>         Spin_quadric_tree_3;

    // spin null graph and configuration space generators based on base predicate type
    template<class Predicate>
    struct Spin_null_graph_3_generator
    {
        typedef CS::Spin_null_graph_3<Kernel, Predicate> Type;
    };

    // spin cell graph and configuration space generators based on base predicate type
    template<class Predicate>
    struct Spin_cell_graph_3_generator
    {
        typedef CS::Spin_cell_graph_3<Kernel, Predicate> Type;
    };

    // spin raster graph and configuration space generators based on base predicate type
    template<class Predicate>
    struct Spin_raster_graph_3_generator
    {
        typedef CS::Spin_raster_graph_3<Kernel, Predicate> Type;
    };

    // configuration space
    template<class Predicate, class Representation>
    struct Spin_configuration_space_3
    {
        typedef CS::Spin_configuration_space_3<Kernel, Predicate, Representation> Type;
    };

    // shortcuts
    template<class Predicate>
    struct Spin_null_configuration_space_3
    {
        typedef CS::Spin_null_configuration_space_3<Kernel, Predicate> Type;
    };

    template<class Predicate>
    struct Spin_cell_configuration_space_3
    {
        typedef CS::Spin_cell_configuration_space_3<Kernel, Predicate> Type;
    };

    template<class Predicate>
    struct Spin_raster_configuration_space_3
    {
        typedef CS::Spin_raster_configuration_space_3<Kernel, Predicate> Type;
    };

    // matrix
    typedef CS::Matrix_44<RT>   Matrix;
};

// basic conversions
inline double to_double(const float &x)
{
    return static_cast<double>(x);
}

inline double to_double(const double &x)
{
    return x;
}

// extern instantiations
extern template struct Spin_inexact_kernel_3<CGAL::Cartesian<double> >;

// default kernel
typedef Spin_inexact_kernel_3<CGAL::Cartesian<double> > Default_inexact_kernel;

} // namespace CS

#include "Spin_inexact_kernel_3.ipp"

#endif // LIBCS_SPIN_INEXACT_KERNEL_H
