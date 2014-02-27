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
#ifndef LIBCS_SPIN_QC_SAMPLE_GENERATOR_3_H
#define LIBCS_SPIN_QC_SAMPLE_GENERATOR_3_H

#include "Spin_straight_sample_generator_3.h"
#include "Spin_3.h"
#include "Uniform_random_spin_3.h"
#include "Random_spin_circle_3.h"
#include <log4cxx/logger.h>
#include <deque>

namespace CS
{
// QC (quadric-circle) sample generator: generate samples in between of quadrics in SO(3)
template<class Kernel, class FT = typename Kernel::FT>
struct Spin_QC_sample_generator_3
{
    // no general implementation
};

// QC sample generator designed for inexact types
template<class Kernel, class FT = typename Kernel::FT>
struct Spin_QC_inexact_sample_generator_3
{
    typedef Spin_3<FT>              Sample;
    typedef std::deque<Sample>      Sample_list;

    // cached samples
    Sample_list                     m_cached_samples;

    log4cxx::LoggerPtr              m_logger;

    Spin_QC_inexact_sample_generator_3();

    template<typename Spin_quadric_iterator>
    void fetch(Spin_quadric_iterator quadrics_begin, Spin_quadric_iterator quadrics_end);

    template<typename Spin_quadric_iterator>
    void operator()(Spin_quadric_iterator quadrics_begin, Spin_quadric_iterator quadrics_end, Sample &sample);
};

// QC sample generator
//
// specialize depending on type exactness
// for now, only specialize QC for inexact types
// the specialization is provided by Spin_QC_inexact_sample_generator_3
template<class Kernel>
struct Spin_QC_sample_generator_3<Kernel, float>
    : Spin_QC_inexact_sample_generator_3<Kernel, float>
{
};

template<class Kernel>
struct Spin_QC_sample_generator_3<Kernel, double>
    : Spin_QC_inexact_sample_generator_3<Kernel, double>
{
};

template<class Kernel>
struct Spin_QC_sample_generator_3<Kernel, CGAL::Gmpfr>
    : Spin_QC_inexact_sample_generator_3<Kernel, CGAL::Gmpfr>
{
};
} // namespace CS

#include "Spin_QC_sample_generator_3.ipp"

#endif // LIBCS_SPIN_QC_SAMPLE_GENERATOR_3_H
