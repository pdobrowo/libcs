###
 # Copyright (C) 2009-2013  Przemysław Dobrowolski
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
project(cs)
cmake_minimum_required(VERSION 2.6)

# configuration: with exact kernel
set(LIBCS_WITH_EXACT_KERNEL 0)
message("with exact kernel: " ${LIBCS_WITH_EXACT_KERNEL})

# CGAL
find_package(CGAL REQUIRED)
include(${CGAL_USE_FILE})

# log4cxx
#find_package(log4cxx REQUIRED)
#include_directories(${log4cxx_INCLUDE_DIRS})

# boost
find_package(Boost COMPONENTS thread REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

# compiler settings
if(UNIX)
    add_definitions(
    #   -Werror
        -Wall
        -Wextra
    #   -pedantic
    #   -Weffc++
        -Wwrite-strings
        -Winit-self
        -Wcast-align
        -Wcast-qual
    #   -Wold-style-cast
        -Wpointer-arith
        -Wstrict-aliasing
        -Wuninitialized
    #   -Wmissing-declarations
    #   -Woverloaded-virtual
    #   -Wnon-virtual-dtor
    #   -Wctor-dtor-privacy
    #   -Wlong-long
        -Wunreachable-code
    #   -Wfloat-equal
        -Wabi
    #   -Wswitch-enum
        -Wformat=2
    #   -Wundef
    #   -Wshadow
    #   -Wconversion
    #   -Wlogical-op
    #   -Waggregate-return
        -Wmissing-field-initializers
        -Wredundant-decls
        -Wmissing-include-dirs
    #   -Wswitch-default
    #   -Wsync-nand
        -Wunused
    #   -Wstrict-overflow=5
    #   -Wunsafe-loop-optimizations
        -Wmissing-format-attribute
    #   -Wpadded
    #   -Winline
        -Wdisabled-optimization
    )
endif(UNIX)

# output directory
set(LIBRARY_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/bin)

# base include directory
include_directories(
    ${CMAKE_SOURCE_DIR}/include
)

# non-packaged 3rdparty external libraries
set(ARCH ${CMAKE_SOURCE_DIR}/3rdparty)

# architecture
if(CMAKE_SIZEOF_VOID_P MATCHES "8")
    set(ARCH ${ARCH}/amd64)
endif(CMAKE_SIZEOF_VOID_P MATCHES "8")

if(CMAKE_SIZEOF_VOID_P MATCHES "4")
    set(ARCH ${ARCH}/i386)
endif(CMAKE_SIZEOF_VOID_P MATCHES "4")

# system
if(WIN32)
    set(ARCH ${ARCH}/windows)
endif()

if(UNIX)
    set(ARCH ${ARCH}/linux)
endif()

message("3rdparty external libraries: " ${ARCH})

if(LIBCS_WITH_EXACT_KERNEL)
    # lidia
    set(LIDIA ${ARCH}/lidia-2.3.0)
    include_directories(${LIDIA}/include)
    link_directories(${LIDIA}/lib)

    # libqi
    set(LIBQI ${ARCH}/libqi-0.9.33)
    include_directories(${LIBQI}/include)
    link_directories(${LIBQI}/lib)

    # realroot
    set(REALROOT ${ARCH}/realroot)
    include_directories(${REALROOT}/include)
    link_directories(${REALROOT}/lib)
endif(LIBCS_WITH_EXACT_KERNEL)

# optimization
#add_definitions(-DCGAL_NDEBUG)

if(WIN32)
    # boost needs that
#    add_definitions(-DBOOST_THREAD_USE_LIB)

    # boost
#    set(BOOST ${ARCH}/boost-1.49)
#    include_directories(${BOOST}/include)
#    link_directories(${BOOST}/lib)

    # GMP
    set(GMP ${ARCH}/gmp-5.0.5)
    include_directories(${GMP}/include)
    link_directories(${GMP}/lib)

    # MPFR
    set(MPFR ${ARCH}/mpfr-3.1.0)
    include_directories(${MPFR}/include)
    link_directories(${MPFR}/lib)

    # MPFI
    set(MPFI ${ARCH}/mpfi-1.5.1)
    include_directories(${MPFI}/include)
    link_directories(${MPFI}/lib)

    # CGAL
#    set(CGAL ${ARCH}/CGAL-4.0)
#    include_directories(${CGAL}/include)
#    link_directories(${CGAL}/lib)

    # APR
    set(APR ${ARCH}/apr-1.4.6)
    include_directories(${APR}/include)
    link_directories(${APR}/lib)

    # APRUTIL
    set(APRUTIL ${ARCH}/apr-util-1.4.1)
    include_directories(${APRUTIL}/include)
    link_directories(${APRUTIL}/lib)

    # LOG4CXX
    set(LOG4CXX ${ARCH}/log4cxx-0.10.0)
    include_directories(${LOG4CXX}/include)
    link_directories(${LOG4CXX}/lib)
endif()

set(cs_HEADERS
    3rdparty/libqi.h
    include/cs/Uspensky_Hanrot.h
    include/cs/Uniform_random_spin_3.ipp
    include/cs/Spin_straight_sample_generator_3.ipp
    include/cs/Spin_quadric_tree_generator.ipp
    include/cs/Spin_quadric_tree_3.ipp
    include/cs/Spin_quadric_mesh_3.ipp
    include/cs/Spin_quadric_3.ipp
    include/cs/Spin_qsip_point.ipp
    include/cs/Spin_qsip_mesh_3.ipp
    include/cs/Spin_qsip_3.ipp
    include/cs/Spin_qsic_3.ipp
    include/cs/Spin_QC_sample_generator_3.ipp
    include/cs/Spin_kernel_3.ipp
    include/cs/Spin_inexact_kernel_3.ipp
    include/cs/Spin_configuration_space_3.ipp
    include/cs/Spin_cell_graph_3.ipp
    include/cs/Spin_3.ipp
    include/cs/Smooth_triangle_3.ipp
    include/cs/Random_spin_circle_3.ipp
    include/cs/Quartic.ipp
    include/cs/Quadratic.ipp
    include/cs/Predicate_tt_3.ipp
    include/cs/Predicate_s_3.ipp
    include/cs/Predicate_h_3.ipp
    include/cs/Predicate_g_3.ipp
    include/cs/Predicate_bb_3.ipp
    include/cs/Polynomial_sign_at.ipp
    include/cs/Polynomial_interval.ipp
    include/cs/Math_utils.ipp
    include/cs/Isolate_polynomial_roots.ipp
    include/cs/Hom_root.ipp
    include/cs/Cubic.ipp
    include/cs/Convert_polynomial.ipp
    include/cs/Cell_n.ipp
    include/cs/Bigint.ipp
    include/cs/Benchmark.ipp
    include/cs/Uniform_random_spin_3.h
    include/cs/Spin_straight_sample_generator_3.h
    include/cs/Spin_quadric_tree_generator.h
    include/cs/Spin_quadric_tree_3.h
    include/cs/Spin_quadric_mesh_3.h
    include/cs/Spin_quadric_3.h
    include/cs/Spin_qsip_point.h
    include/cs/Spin_qsip_mesh_3.h
    include/cs/Spin_qsip_3.h
    include/cs/Spin_qsic_mesh_3.h
    include/cs/Spin_qsic_3.h
    include/cs/Spin_QC_sample_generator_3.h
    include/cs/Spin_kernel_3.h
    include/cs/Spin_inexact_kernel_3.h
    include/cs/Spin_configuration_space_3.h
    include/cs/Spin_cell_graph_3.h
    include/cs/Spin_3.h
    include/cs/Sphere_triangle_intersection_3.h
    include/cs/Smooth_triangle_3.h
    include/cs/Random_spin_circle_3.h
    include/cs/Random.h
    include/cs/Quartic.h
    include/cs/Quadratic.h
    include/cs/Predicate_tt_3.h
    include/cs/Predicate_s_3.h
    include/cs/Predicate_h_3.h
    include/cs/Predicate_g_3.h
    include/cs/Predicate_bb_3.h
    include/cs/Polynomial_sign_at.h
    include/cs/Polynomial_interval.h
    include/cs/Math_utils.h
    include/cs/Isolate_polynomial_roots.h
    include/cs/Hom_root.h
    include/cs/Extended.h
    include/cs/Cubic.h
    include/cs/Coordinate.h
    include/cs/Convert_polynomial.h
    include/cs/Contfrac.h
    include/cs/Cell_n.h
    include/cs/Bigint.h
    include/cs/Benchmark.h
    include/cs/Sphere_triangle_intersection_3.ipp
    include/cs/Spin_qsic_mesh_3.ipp
    include/cs/Index_tree_n.ipp
    include/cs/Index_tree_n.h
    include/cs/Predicate_list_generator.h
    include/cs/Predicate_tt_3_list_generator.ipp
    include/cs/Predicate_tt_3_list_generator.h
    include/cs/Predicate_bb_3_list_generator.ipp
    include/cs/Predicate_bb_3_list_generator.h
    include/cs/Ball_3.ipp
    include/cs/Ball_3.h
    include/cs/Loader_sphere_tree.ipp
    include/cs/Loader_sphere_tree.h
    include/cs/Loader_scene.h
    include/cs/Loader_scene.ipp
    include/cs/Build.h
    include/cs/Spin_raster_graph_3.ipp
    include/cs/Spin_raster_graph_3.h
    include/cs/Voxel_3.h
    include/cs/Voxel_3.ipp
    include/cs/Spin_exact_graph_3.ipp
    include/cs/Spin_exact_graph_3.h
    include/cs/Config.h
)

set(cs_BASE_SOURCES
    ${cs_HEADERS}
    src/Spin_inexact_kernel_3.cpp
    src/Spin_3.cpp
    src/Scene_loader.cpp
    src/Random.cpp
    src/Math_utils.cpp
    src/Coordinate.cpp
    src/Benchmark.cpp
    src/Build.cpp
    src/Config.cpp
)

if(LIBCS_WITH_EXACT_KERNEL)
    set(cs_SOURCES
        ${cs_BASE_SOURCES}
        src/Spin_kernel_3.cpp
        src/Uspensky_Hanrot.c
        src/Hom_root.cpp
        src/Contfrac.cpp
    )
else(LIBCS_WITH_EXACT_KERNEL)
    set(cs_SOURCES
        ${cs_BASE_SOURCES}
    )
endif(LIBCS_WITH_EXACT_KERNEL)

add_library(cs SHARED ${cs_SOURCES})

target_link_libraries(cs 
    ${CGAL_LIBRARIES}
    ${Boost_LIBRARIES}
    ${log4cxx_LIBRARIES}
)

if(LIBCS_WITH_EXACT_KERNEL)
    target_link_libraries(cs
        LiDIA
        qi
        mpfi
        mpfr
        gmp
        realroot
    )
endif(LIBCS_WITH_EXACT_KERNEL)

if(WIN32)
    target_link_libraries(cs
        apr-1
        aprutil-1
        ws2_32
        expat
        iconv
        mswsock
    )
endif()

