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
project(find_path)
cmake_minimum_required(VERSION 2.6)

# cmake find
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/../../cmake)

# cs
find_package(CS REQUIRED)
include_directories(${CS_INCLUDE_DIR})

# compiler settings
SET(DEFAULT_CLANG_AND_GNU_FLAGS "-std=c++14 -Wall -Wextra")

IF(${CMAKE_CXX_COMPILER_ID} MATCHES Clang)
ADD_DEFINITIONS(${DEFAULT_CLANG_AND_GNU_FLAGS})
ENDIF(${CMAKE_CXX_COMPILER_ID} MATCHES Clang)

IF(${CMAKE_CXX_COMPILER_ID} MATCHES GNU)
ADD_DEFINITIONS(${DEFAULT_CLANG_AND_GNU_FLAGS})
ENDIF(${CMAKE_CXX_COMPILER_ID} MATCHES GNU)

# header
set(find_path_HEADERS
    find_path.h
)

# source
set(find_path_SOURCES
    ${find_path_HEADERS}
    main.cpp
    find_path.cpp
)

# bin
add_executable(find_path ${find_path_SOURCES})
target_link_libraries(find_path ${CS_LIBRARIES})
