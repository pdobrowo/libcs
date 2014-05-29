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
project(test)
cmake_minimum_required(VERSION 2.6)

# configuration: with exact kernel
set(LIBCS_WITH_EXACT_KERNEL $ENV{EXACT})

if(LIBCS_WITH_EXACT_KERNEL)
    message("ok")
else(LIBCS_WITH_EXACT_KERNEL)
    set(LIBCS_WITH_EXACT_KERNEL 0)
endif(LIBCS_WITH_EXACT_KERNEL)

message("with exact kernel: " ${LIBCS_WITH_EXACT_KERNEL})

# CGAL
find_package(CGAL REQUIRED)
include(${CGAL_USE_FILE})

# boost
find_package(Boost COMPONENTS thread REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

# cs
set(cs_INCLUDE_DIRS ${CMAKE_SOURCE_DIR}/../include)
set(cs_LINK_DIRS ${CMAKE_SOURCE_DIR}/../bin)
set(cs_LIBRARIES cs)

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
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/../bin)

# base include directory
include_directories(
    ${cs_INCLUDE_DIRS}
)

link_directories(
    ${cs_LINK_DIRS}
)

# non-packaged 3rdparty external libraries
set(ARCH ${CMAKE_SOURCE_DIR}/../3rdparty)

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

set(test_HEADERS
    test_inexact_perf.h
)

set(test_BASE_SOURCES
    ${test_HEADERS}
    main.cpp
    test_inexact_perf.cpp
    test_various.cpp
)

if(LIBCS_WITH_EXACT_KERNEL)
    set(test_SOURCES
        ${test_BASE_SOURCES}
        
    )
else(LIBCS_WITH_EXACT_KERNEL)
    set(test_SOURCES
        ${test_BASE_SOURCES}
    )
endif(LIBCS_WITH_EXACT_KERNEL)

add_executable(test ${test_SOURCES})

target_link_libraries(test
    ${cs_LIBRARIES}
    ${CGAL_LIBRARIES}
    ${Boost_LIBRARIES}
    log4cxx
)

if(LIBCS_WITH_EXACT_KERNEL)
    target_link_libraries(test
        LiDIA
        qi
        mpfi
        mpfr
        gmp
        realroot
    )
endif(LIBCS_WITH_EXACT_KERNEL)

if(WIN32)
    target_link_libraries(test
        apr-1
        aprutil-1
        ws2_32
        expat
        iconv
        mswsock
    )
endif()

