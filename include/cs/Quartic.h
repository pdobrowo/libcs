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
#ifndef LIBCS_QUARTIC_H
#define LIBCS_QUARTIC_H

#include "Quadratic.h"
#include <CGAL/Gmpfr.h>
#include <algorithm>
#include <cassert>
#include <cmath>

namespace CS
{
namespace Math
{
/*
 * Original version generalized to templated one
 *
 * The code comes from GSL package extension, which in fact comes from CERN project
 */
using Math::sin;
using Math::cos;
using Math::sqrt;
using Math::acos;
using Math::pow;
using Math::atan2;

using Math::pi;

/* poly/solve_quartic.c
 *
 * Copyright (C) 2003 CERN and K.S. K\"{o}lbig
 *
 * Converted to C and implemented into the GSL Library
 * by Andrew W. Steiner and Andy Buckley
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

/* solve_quartic.c - finds the real roots of
 *  x^4 + a x^3 + b x^2 + c x + d = 0
 */
template<class FT_>
int solve_quartic(const FT_ &a, const FT_ &b, const FT_ &c, const FT_ &d,
                  FT_ *x0, FT_ *x1, FT_ *x2, FT_ *x3);

/*  DYNAMO:- Event driven molecular dynamics simulator
    http://www.marcusbannerman.co.uk/dynamo
    Copyright (C) 2010  Marcus N Campbell Bannerman <m.bannerman@gmail.com>

    This program is free software: you can redistribute it and/or
    modify it under the terms of the GNU General Public License
    version 3 as published by the Free Software Foundation.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#if 0
#include <quartic.hpp>

// BUG: this implementation gives ocassionaly errors
// BUG: maybe the problematic is a my bugfix in ferrari code
template<>
int solve_quartic<double>(const double &a, const double &b, const double &c, const double &d,
                          double *x0, double *x1, double *x2, double *x3)
{
    return static_cast<int>(magnet::math::quarticSolve(a, b, c, d, *x0, *x1, *x2, *x3));
}
#endif
} // namespace Math
} // namespace CS

#include "Quartic.ipp"

#endif // LIBCS_QUARTIC_H
