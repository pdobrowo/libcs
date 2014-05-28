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
#ifndef LIBCS_UNIFORM_RANDOM_SPIN_3_H
#define LIBCS_UNIFORM_RANDOM_SPIN_3_H

#include <CGAL/Gmpfr.h>
#include "Spin_3.h"
#include "Math_utils.h"
#include <cstdlib>
#include <cassert>

namespace CS
{
template<typename FT_>
FT_ uniform_rand();

template<class FT_>
void uniform_random_spin_3(Spin_3<FT_> &out);
} // namespace CS

#include "Uniform_random_spin_3.ipp"

#endif // LIBCS_UNIFORM_RANDOM_SPIN_3_H
