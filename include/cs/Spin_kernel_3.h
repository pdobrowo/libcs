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
#ifndef LIBCS_SPIN_KERNEL_H
#define LIBCS_SPIN_KERNEL_H

// define this to use external RS library for root isolation
//#define CS_USE_RS_LIBRARY

// algebraic kernel
#if 1
#   define CGAL_USE_MPFI 1
#endif

#include <CGAL/basic.h>
#include <CGAL/Gmpz.h>

// math utils
#define LIBCS_SPIN_KERNEL_INCLUDE_FILE
#include "Math_utils_exact.h"
#undef LIBCS_SPIN_KERNEL_INCLUDE_FILE

#ifdef CS_USE_RS_LIBRARY
#    define CGAL_USE_RS3 1
//   RS kernel for integer polynomials
#    include <CGAL/RS/Algebraic_kernel_rs_1.h>
#    include <CGAL/Algebraic_kernel_rs_gmpz_d_1.h>
//   standard kernel for sqrt-extended polynomials
#    include <CGAL/Algebraic_kernel_d_1.h>
#else
#    include <CGAL/Algebraic_kernel_d_1.h>
#endif

#include <CGAL/Quotient.h>
#include <CGAL/Sqrt_extension.h>
#include <CGAL/Filtered_kernel.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include "Bigint.h"

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
#include "Spin_qsic_3.h"
#include "Spin_qsip_3.h"
#include "Spin_qsip_point.h"
#include "Hom_root.h"

// meshers
#include "Spin_quadric_mesh_3.h"
#include "Spin_qsic_mesh_3.h"
#include "Spin_qsip_mesh_3.h"
#include "Smooth_triangle_3.h"

// sentence
#include "Spin_quadric_tree_3.h"

// configuration space
#include "Spin_null_graph_3.h"
#include "Spin_cell_graph_3.h"
#include "Spin_raster_graph_3.h"
#include "Spin_exact_graph_3.h"
#include "Spin_configuration_space_3.h"

// matrix
#include "Matrix_44.h"
#include <lidia/matrix.h>

// qsic
#include <libqi.h>

namespace CS
{
// Spin_kernel_3 is an adapter for a standard CGAL geometric kernel
// this is an exact version
template<class K_>
struct Spin_kernel_3
    : public K_
{
    // derived types
    typedef typename K_::RT                         RT;
    typedef CGAL::Quotient<RT>                      FT;

    // spin kernel
    typedef Spin_kernel_3                           Kernel;

    // spin types
    typedef CS::Spin_quadric_3<Kernel>              Spin_quadric_3;
    typedef CS::Spin_reduced_quadric_3<Kernel>      Spin_reduced_quadric_3;
    typedef CS::Spin_qsic_3<Kernel>                 Spin_qsic_3;
    typedef CS::Spin_qsip_3<Kernel>                 Spin_qsip_3;

    // QSIC
    typedef QI::quad_inter<RT>                      Qsic;
    typedef QI::component<RT>                       Qsic_component;
    typedef QI::curve_param<RT>                     Qsic_curve;
    typedef QI::surface_param<RT>                   Qsic_surface;

    // polynomial
    typedef QI::hom_polynomial<RT>                  Hom_polynomial;
    typedef QI::hom_hom_polynomial<RT>              Hom_hom_polynomial;

    typedef CGAL::Sqrt_extension<RT, RT>            Sqrt_extension;

    typedef QI::hom_polynomial<Sqrt_extension>      Hom_polynomial_with_sqrt;

    typedef CS::Hom_root<Kernel>                    Hom_root;
    typedef CS::Spin_qsip_point<Kernel>             Spin_qsip_point;

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
    typedef CS::Spin_qsic_mesh_3<Kernel>            Spin_qsic_mesh_3;
    typedef CS::Spin_qsip_mesh_3<Kernel>            Spin_qsip_mesh_3;
    typedef CS::Smooth_triangle_3<Kernel>           Smooth_triangle_3;

    // exact spin
    typedef CS::Spin_3<FT>                          Spin_3;
    typedef CS::Diff_spin_3<FT>                     Diff_spin_3;

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

    // spin exact graph and configuration space generators based on base predicate type
    template<class Predicate>
    struct Spin_exact_graph_3_generator
    {
        typedef CS::Spin_exact_graph_3<Kernel, Predicate> Type;
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

    template<class Predicate>
    struct Spin_exact_configuration_space_3
    {
        typedef CS::Spin_exact_configuration_space_3<Kernel, Predicate> Type;
    };

    // algebraic kernel
#ifdef CS_USE_RS_LIBRARY
    typedef Algebraic_kernel_rs_1<CGAL::Gmpz>                                           Algebraic_kernel;
//  typedef Algebraic_kernel_rs_1<CGAL::Sqrt_extension<CGAL::Gmpz, CGAL::Gmpz> >        Algebraic_kernel_with_sqrt;
    typedef CGAL::Algebraic_kernel_d_1<CGAL::Sqrt_extension<CGAL::Gmpz, CGAL::Gmpz> >   Algebraic_kernel_with_sqrt;
#else
    typedef CGAL::Algebraic_kernel_d_1<CGAL::Gmpz>                                      Algebraic_kernel;
    typedef CGAL::Algebraic_kernel_d_1<CGAL::Sqrt_extension<CGAL::Gmpz, CGAL::Gmpz> >   Algebraic_kernel_with_sqrt;
#endif

    // matrix
    typedef CS::Matrix_44<RT>   Matrix;
    typedef LiDIA::matrix<RT>   LidiaMatrix;

    static LidiaMatrix to_lidia_matrix(const Matrix &matrix)
    {
        LidiaMatrix result(4, 4);

        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 4; ++j)
                result.sto(i, j, matrix.get(i, j));

        return result;
    }
};

// default kernel
typedef Spin_kernel_3<
            CGAL::Filtered_kernel<
                CGAL::Cartesian<
                    CGAL::Bigint>
                >
            >
        Default_kernel;
} // namespace CS

#include "Spin_kernel_3.ipp"

#endif // LIBCS_SPIN_KERNEL_H
