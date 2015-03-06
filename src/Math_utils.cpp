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
#define LIBCS_MATH_UTILS_INTERNAL_FILE
#include <cs/Math_utils.h>
#undef LIBCS_MATH_UTILS_INTERNAL_FILE

#include <cmath>

namespace CS
{
namespace Math
{
template<>
float pi<float>()
{
    return static_cast<float>(M_PI);
}

template<>
double pi<double>()
{
    return M_PI;
}

float sin(const float &x)
{
    return std::sin(x);
}

float cos(const float &x)
{
    return std::cos(x);
}

float acos(const float &x)
{
    return std::acos(x);
}

float atan(const float &x)
{
    return std::atan(x);
}

float pow(const float &x, const float &y)
{
    return std::pow(x, y);
}

float atan2(const float &y, const float &x)
{
    return std::atan2(y, x);
}

float sqrt(const float &x)
{
    return std::sqrt(x);
}

float fabs(const float &x)
{
    return std::fabs(x);
}

double sin(const double &x)
{
    return std::sin(x);
}

double cos(const double &x)
{
    return std::cos(x);
}

double acos(const double &x)
{
    return std::acos(x);
}

double atan(const double &x)
{
    return std::atan(x);
}

double pow(const double &x, const double &y)
{
    return std::pow(x, y);
}

double atan2(const double &y, const double &x)
{
    return std::atan2(y, x);
}

double sqrt(const double &x)
{
    return std::sqrt(x);
}

double fabs(const double &x)
{
    return std::fabs(x);
}

} // namespace Math
} // namespace CS
