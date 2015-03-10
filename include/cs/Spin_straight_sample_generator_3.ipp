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
#include "Spin_straight_sample_generator_3.h"

namespace CS
{
template<class Kernel_, class FT_>
template<typename Spin_quadric_iterator_>
void Spin_straight_sample_generator_3<Kernel_, FT_>::operator()(Spin_quadric_iterator_ quadrics_begin, Spin_quadric_iterator_ quadrics_end, Sample &sample)
{
    (void)quadrics_begin;
    (void)quadrics_end;

    const double SAMPLE_SPIN_DIGITS = 1000000;
    Spin_3<double> sampleSpin;

    // random spin components
    FT s12, s23, s31, sqr;

    // randomize
    do
    {
        // choose uniform random sample
        uniform_random_spin_3(sampleSpin);

        s12 = FT(SAMPLE_SPIN_DIGITS * sampleSpin.s12()) / FT(SAMPLE_SPIN_DIGITS);
        s23 = FT(SAMPLE_SPIN_DIGITS * sampleSpin.s23()) / FT(SAMPLE_SPIN_DIGITS);
        s31 = FT(SAMPLE_SPIN_DIGITS * sampleSpin.s31()) / FT(SAMPLE_SPIN_DIGITS);

        sqr = s12 * s12 + s23 * s23 + s31 * s31;

        // it is possible that the square root is larger that 1 because of floating
        // point rounding errors in that case choose random sample once again
    } while (sqr > FT(1));

    // choose exact spin sample
    sample = Sample(s12, s23, s31, ExtendedFT(FT(0), (rand() % 2) ? FT(1) : FT(-1), FT(1) - sqr));
}

template<class Kernel_>
template<typename Spin_quadric_iterator>
void Spin_straight_sample_generator_3<Kernel_, double>::operator()(Spin_quadric_iterator quadrics_begin, Spin_quadric_iterator quadrics_end, Sample &sample)
{
    (void)quadrics_begin;
    (void)quadrics_end;

    // choose inexact spin sample
    uniform_random_spin_3(sample);
}

template<class Kernel_>
template<typename Spin_quadric_iterator>
void Spin_straight_sample_generator_3<Kernel_, CGAL::Gmpfr>::operator()(Spin_quadric_iterator quadrics_begin, Spin_quadric_iterator quadrics_end, Sample &sample)
{
    (void)quadrics_begin;
    (void)quadrics_end;

    // choose inexact spin sample
    uniform_random_spin_3(sample);
}
} // namespace CS
