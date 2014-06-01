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
#ifndef LIBCS_MATH_UTILS_H
#define LIBCS_MATH_UTILS_H

#if !defined(LIBCS_MATH_UTILS_INTERNAL_FILE)
#if !defined(LIBCS_SPIN_KERNEL_INCLUDE_FILE) && !defined(LIBCS_SPIN_INEXACT_KERNEL_INCLUDE_FILE)
#error Do not include this file directly, include Spin_kernel_3.h or Spin_inexact_kernel_3.h instead
#endif // !defined(LIBCS_SPIN_KERNEL_INCLUDE_FILE) && !defined(LIBCS_SPIN_INEXACT_KERNEL_INCLUDE_FILE)
#endif // !defined(LIBCS_MATH_UTILS_INTERNAL_FILE)

namespace CS
{
namespace Math
{
// add missing definitions in some toolchains
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif // M_PI

float sin(const float &x);
float cos(const float &x);
float acos(const float &x);
float atan(const float &x);
float pow(const float &x, const float &y);
float atan2(const float &y, const float &x);
float sqrt(const float &x);
float fabs(const float &x);

double sin(const double &x);
double cos(const double &x);
double acos(const double &x);
double atan(const double &x);
double pow(const double &x, const double &y);
double atan2(const double &y, const double &x);
double sqrt(const double &x);
double fabs(const double &x);

template<class FT_>
FT_ pi();

template<>
float pi<float>();

template<>
double pi<double>();
} // namespace Math
} // namespace CS

#include "Math_utils.ipp"

#endif // LIBCS_MATH_UTILS_H
