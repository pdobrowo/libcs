###
 # Copyright (C) 2009-2015  Przemysław Dobrowolski
 #
 # This file is part of the Configuration Space Library (libcs), a library
 # for creating configuration spaces of various motion planning problems.
 #
 # This program is free software: you can redistribute it and/or modify
 # it under the terms of the GNU General Public License as published by
 # the Free Software Foundation, either version 3 of the License, or
 # (at your option) any later version.
 #
 # This program is distributed in the hope that it will be useful,
 # but WITHOUT ANY WARRANTY; without even the implied warranty of
 # MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 # GNU General Public License for more details.
 #
 # You should have received a copy of the GNU General Public License
 # along with this program.  If not, see <http://www.gnu.org/licenses/>.
 ##
PROJECT(cs)
CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

# cmake find
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# cgal
FIND_PACKAGE(CGAL REQUIRED)
INCLUDE(${CGAL_USE_FILE})

# boost
FIND_PACKAGE(Boost COMPONENTS thread REQUIRED)
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})

# log4cxx
find_package(LOG4CXX REQUIRED)
INCLUDE_DIRECTORIES(${LOG4CXX_INCLUDE_DIRS})

# compiler settings
SET(LIBCS_DEFAULT_CLANG_AND_GNU_FLAGS "-std=c++14 -Wall -Wextra")

IF(${CMAKE_CXX_COMPILER_ID} MATCHES Clang)
    ADD_DEFINITIONS(${LIBCS_DEFAULT_CLANG_AND_GNU_FLAGS})
ENDIF(${CMAKE_CXX_COMPILER_ID} MATCHES Clang)

IF(${CMAKE_CXX_COMPILER_ID} MATCHES GNU)
    ADD_DEFINITIONS(${LIBCS_DEFAULT_CLANG_AND_GNU_FLAGS})
ENDIF(${CMAKE_CXX_COMPILER_ID} MATCHES GNU)

# output directory
SET(LIBRARY_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/lib)

# base include directory
include_directories(
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/3rd
)

# configuration: with exact kernel
SET(LIBCS_WITH_EXACT_KERNEL $ENV{EXACT})

IF(LIBCS_WITH_EXACT_KERNEL)
    MESSAGE(STATUS "LIBCS: environment variable EXACT is set")
ELSE(LIBCS_WITH_EXACT_KERNEL)
    SET(LIBCS_WITH_EXACT_KERNEL 0)
ENDIF(LIBCS_WITH_EXACT_KERNEL)

MESSAGE(STATUS "LIBCS: exact kernel: " ${LIBCS_WITH_EXACT_KERNEL})

# non-packaged 3rd libs
SET(ARCH ${CMAKE_SOURCE_DIR}/3rd)

# architecture
IF(CMAKE_SIZEOF_VOID_P MATCHES 8)
    SET(ARCH ${ARCH}/amd64)
ENDIF(CMAKE_SIZEOF_VOID_P MATCHES 8)

IF(CMAKE_SIZEOF_VOID_P MATCHES 4)
    SET(ARCH ${ARCH}/i386)
ENDIF(CMAKE_SIZEOF_VOID_P MATCHES 4)

# system
IF(WIN32)
    SET(ARCH ${ARCH}/windows)
ENDIF(WIN32)

IF(UNIX)
    SET(ARCH ${ARCH}/unix)
ENDIF(UNIX)

message(STATUS "LIBCS: 3rd libs: " ${ARCH})

IF(LIBCS_WITH_EXACT_KERNEL)
    # lidia
    SET(LIDIA ${ARCH}/lidia-2.3.0)
    include_directories(${LIDIA}/include)

    # libqi
    SET(LIBQI ${ARCH}/libqi-0.9.33)
    include_directories(${LIBQI}/include)

    # realroot
    SET(REALROOT ${ARCH}/realroot)
    include_directories(${REALROOT}/include)
ENDIF(LIBCS_WITH_EXACT_KERNEL)

# optimization
#add_definitions(-DCGAL_NDEBUG)

IF(WIN32)
    # boost needs that
#    add_definitions(-DBOOST_THREAD_USE_LIB)

    # boost
#    SET(BOOST ${ARCH}/boost-1.49)
#    include_directories(${BOOST}/include)

    # GMP
    SET(GMP ${ARCH}/gmp-5.0.5)
    include_directories(${GMP}/include)

    # MPFR
    SET(MPFR ${ARCH}/mpfr-3.1.0)
    include_directories(${MPFR}/include)

    # MPFI
    SET(MPFI ${ARCH}/mpfi-1.5.1)
    include_directories(${MPFI}/include)

    # CGAL
#    SET(CGAL ${ARCH}/CGAL-4.0)
#    include_directories(${CGAL}/include)

    # APR
    SET(APR ${ARCH}/apr-1.4.6)
    include_directories(${APR}/include)

    # APRUTIL
    SET(APRUTIL ${ARCH}/apr-util-1.4.1)
    include_directories(${APRUTIL}/include)

    # LOG4CXX
    SET(LOG4CXX ${ARCH}/log4cxx-0.10.0)
    include_directories(${LOG4CXX}/include)
ENDIF(WIN32)

SET(cs_HEADERS
    3rd/libqi.h
    include/cs/Ball_3.h
    include/cs/Ball_3.ipp
    include/cs/Benchmark.h
    include/cs/Benchmark.ipp
    include/cs/Bigint.h
    include/cs/Bigint.ipp
    include/cs/Build.h
    include/cs/Cell_n.h
    include/cs/Cell_n.ipp
    include/cs/Config.h
    include/cs/Contfrac.h
    include/cs/Convert_polynomial.h
    include/cs/Convert_polynomial.ipp
    include/cs/Coordinate.h
    include/cs/Cubic.h
    include/cs/Cubic.ipp
    include/cs/Extended.h
    include/cs/Hom_root.h
    include/cs/Hom_root.ipp
    include/cs/Index_tree_n.h
    include/cs/Index_tree_n.ipp
    include/cs/Isolate_polynomial_roots.h
    include/cs/Isolate_polynomial_roots.ipp
    include/cs/Linear_system.h
    include/cs/Linear_system.ipp
    include/cs/Loader_scene.h
    include/cs/Loader_scene.ipp
    include/cs/Loader_sphere_tree.h
    include/cs/Loader_sphere_tree.ipp
    include/cs/Math_utils.h
    include/cs/Math_utils_exact.h
    include/cs/Math_utils.ipp
    include/cs/Math_utils_exact.ipp
    include/cs/Matrix_44.h
    include/cs/Matrix_44.ipp
    include/cs/Polynomial_interval.h
    include/cs/Polynomial_interval.ipp
    include/cs/Polynomial_sign_at.h
    include/cs/Polynomial_sign_at.ipp
    include/cs/Predicate_bb_3.h
    include/cs/Predicate_bb_3.ipp
    include/cs/Predicate_bb_3_list_generator.h
    include/cs/Predicate_bb_3_list_generator.ipp
    include/cs/Predicate_h_3_list_generator.h
    include/cs/Predicate_h_3_list_generator.ipp
    include/cs/Predicate_g_3.h
    include/cs/Predicate_g_3.ipp
    include/cs/Predicate_g_parametrization_3.h
    include/cs/Predicate_g_parametrization_3.ipp
    include/cs/Predicate_h_3.h
    include/cs/Predicate_h_3.ipp
    include/cs/Predicate_list_generator.h
    include/cs/Predicate_s_3.h
    include/cs/Predicate_s_3.ipp
    include/cs/Predicate_tt_3.h
    include/cs/Predicate_tt_3.ipp
    include/cs/Predicate_tt_3_list_generator.h
    include/cs/Predicate_tt_3_list_generator.ipp
    include/cs/Quadratic.h
    include/cs/Quadratic.ipp
    include/cs/Quartic.h
    include/cs/Quartic.ipp
    include/cs/Random.h
    include/cs/Random_spin_circle_3.h
    include/cs/Random_spin_circle_3.ipp
    include/cs/Realroot.h
    include/cs/Smooth_triangle_3.h
    include/cs/Smooth_triangle_3.ipp
    include/cs/Sphere_triangle_intersection_3.h
    include/cs/Sphere_triangle_intersection_3.ipp
    include/cs/Spin_3.h
    include/cs/Spin_3.ipp
    include/cs/Spin_cell_graph_3.h
    include/cs/Spin_cell_graph_3.ipp
    include/cs/Spin_configuration_space_3.h
    include/cs/Spin_configuration_space_3.ipp
    include/cs/Spin_exact_graph_3.h
    include/cs/Spin_exact_graph_3.ipp
    include/cs/Spin_inexact_kernel_3.h
    include/cs/Spin_inexact_kernel_3.ipp
    include/cs/Spin_kernel_3.h
    include/cs/Spin_kernel_3.ipp
    include/cs/Spin_QC_sample_generator_3.h
    include/cs/Spin_QC_sample_generator_3.ipp
    include/cs/Spin_qsic_3.h
    include/cs/Spin_qsic_3.ipp
    include/cs/Spin_qsic_mesh_3.h
    include/cs/Spin_qsic_mesh_3.ipp
    include/cs/Spin_qsic_trim_3.h
    include/cs/Spin_qsip_3.h
    include/cs/Spin_qsip_3.ipp
    include/cs/Spin_qsip_mesh_3.h
    include/cs/Spin_qsip_mesh_3.ipp
    include/cs/Spin_qsip_point.h
    include/cs/Spin_qsip_point.ipp
    include/cs/Spin_quadric_3.h
    include/cs/Spin_quadric_3.ipp
    include/cs/Spin_quadric_mesh_3.h
    include/cs/Spin_quadric_mesh_3.ipp
    include/cs/Spin_quadric_tree_3.h
    include/cs/Spin_quadric_tree_3.ipp
    include/cs/Spin_quadric_tree_generator.h
    include/cs/Spin_quadric_tree_generator.ipp
    include/cs/Spin_raster_graph_3.h
    include/cs/Spin_raster_graph_3.ipp
    include/cs/Spin_null_graph_3.h
    include/cs/Spin_null_graph_3.ipp
    include/cs/Spin_reduced_quadric_3.h
    include/cs/Spin_reduced_quadric_3.ipp
    include/cs/Spin_straight_sample_generator_3.h
    include/cs/Spin_straight_sample_generator_3.ipp
    include/cs/Uniform_random_spin_3.h
    include/cs/Uniform_random_spin_3.ipp
    include/cs/Uspensky.h
    include/cs/Uspensky_Hanrot.h
    include/cs/Vector_4.h
    include/cs/Vector_4.ipp
    include/cs/Voxel_3.h
    include/cs/Voxel_3.ipp
)

SET(cs_BASE_SOURCES
    ${cs_HEADERS}
    src/Spin_inexact_kernel_3.cpp
    src/Spin_qsic_trim_3.cpp
    src/Spin_3.cpp
    src/Scene_loader.cpp
    src/Random.cpp
    src/Math_utils.cpp
    src/Coordinate.cpp
    src/Benchmark.cpp
    src/Build.cpp
    src/Config.cpp
)

IF(LIBCS_WITH_EXACT_KERNEL)
    SET(cs_SOURCES
        ${cs_BASE_SOURCES}
        src/Spin_kernel_3.cpp
        src/Uspensky_Hanrot.c
        src/Hom_root.cpp
        src/Contfrac.cpp
        src/Math_utils_exact.cpp
    )
ELSE(LIBCS_WITH_EXACT_KERNEL)
    SET(cs_SOURCES
        ${cs_BASE_SOURCES}
    )
ENDIF(LIBCS_WITH_EXACT_KERNEL)

add_library(cs STATIC ${cs_SOURCES})
