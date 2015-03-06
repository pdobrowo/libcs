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
#ifndef LIBCS_QUADRATIC_H
#define LIBCS_QUADRATIC_H

#include <CGAL/Gmpfr.h>
#include <cassert>
#include <cmath>

namespace CS
{
namespace Math
{
/*
 * Original version generalized to templated one
 *
 * The code comes from GSL package
 */
using Math::sqrt;
using Math::fabs;

/* poly/solve_quadratic.c
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2007 Brian Gough
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/* solve_quadratic.c - finds the real roots of a x^2 + b x + c = 0 */

template<class FT_>
int solve_quadratic(const FT_ &a, const FT_ &b, const FT_ &c,
                    FT_ *x0, FT_ *x1);

} // namespace Math
} // namespace CS

#include "Quadratic.ipp"

#endif // LIBCS_QUADRATIC_H
