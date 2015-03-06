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
#ifndef LIBCS_CONTFRAC_H
#define LIBCS_CONTFRAC_H

#include <CGAL/Gmpq.h>
#include <gmp.h>
#include <vector>

struct contfrac_interval
{
    CGAL::Gmpq low;
    CGAL::Gmpq high;

    contfrac_interval(const CGAL::Gmpq &low_, const CGAL::Gmpq &high_);
};

// real root isolation based on continued fraction expansion
void contfrac(mpz_t *coeffs, unsigned int degree, std::vector<contfrac_interval> &out);

#endif // LIBCS_CONTFRAC_H
