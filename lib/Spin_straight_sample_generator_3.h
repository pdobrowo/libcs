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
#ifndef LIBCS_SPIN_STRAIGHT_SAMPLE_GENERATOR_3_H
#define LIBCS_SPIN_STRAIGHT_SAMPLE_GENERATOR_3_H

#include "Extended.h"
#include "Spin_3.h"
#include <CGAL/Gmpfr.h>
#include "Uniform_random_spin_3.h"

namespace CS
{
// straight sample generator: choose spin sample uniformly in SO(3)
template<class Kernel, class FT = typename Kernel::FT>
struct Spin_straight_sample_generator_3
{
    typedef typename Extended_generator<FT>::Type   ExtendedFT;
    typedef Spin_3<ExtendedFT>                      Sample;

    template<typename Spin_quadric_iterator>
    void operator()(Spin_quadric_iterator quadrics_begin, Spin_quadric_iterator quadrics_end, Sample &sample);
};

// TODO: the same implementation is for all inexact types
//       create inexact template implementation, and derive this implementation in all inexact types
//       see QC sample generator for reference
template<class Kernel>
struct Spin_straight_sample_generator_3<Kernel, double>
{
    typedef Spin_3<double>          Sample;

    template<typename Spin_quadric_iterator>
    void operator()(Spin_quadric_iterator quadrics_begin, Spin_quadric_iterator quadrics_end, Sample &sample);
};

template<class Kernel>
struct Spin_straight_sample_generator_3<Kernel, CGAL::Gmpfr>
{
    typedef Spin_3<CGAL::Gmpfr>     Sample;

    template<typename Spin_quadric_iterator>
    void operator()(Spin_quadric_iterator quadrics_begin, Spin_quadric_iterator quadrics_end, Sample &sample);
};
} // namespace CS

#include "Spin_straight_sample_generator_3.ipp"

#endif // LIBCS_SPIN_STRAIGHT_SAMPLE_GENERATOR_3_H
