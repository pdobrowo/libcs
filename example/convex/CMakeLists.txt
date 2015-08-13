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
PROJECT(convex)
CMAKE_MINIMUM_REQUIRED(VERSION 2.6)

# cmake find
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/../../cmake)

# cs
FIND_PACKAGE(CS REQUIRED)
INCLUDE_DIRECTORIES(${CS_INCLUDE_DIR})

# compiler settings
SET(DEFAULT_CLANG_AND_GNU_FLAGS "-std=c++11 -Wall -Wextra")

IF(${CMAKE_CXX_COMPILER_ID} MATCHES Clang)
ADD_DEFINITIONS(${DEFAULT_CLANG_AND_GNU_FLAGS})
ENDIF(${CMAKE_CXX_COMPILER_ID} MATCHES Clang)

IF(${CMAKE_CXX_COMPILER_ID} MATCHES GNU)
ADD_DEFINITIONS(${DEFAULT_CLANG_AND_GNU_FLAGS})
ENDIF(${CMAKE_CXX_COMPILER_ID} MATCHES GNU)

# header
SET(convex_HEADERS
    convex.h
)

# source
SET(convex_SOURCES
    ${convex_HEADERS}
    main.cpp
    convex.cpp
)

# bin
ADD_EXECUTABLE(convex ${convex_SOURCES})
TARGET_LINK_LIBRARIES(convex ${CS_LIBRARIES})
