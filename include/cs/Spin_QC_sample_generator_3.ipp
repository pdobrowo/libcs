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
#include "Spin_QC_sample_generator_3.h"

namespace CS
{
template<class Kernel_, class FT_>
Spin_QC_inexact_sample_generator_3<Kernel_, FT_>::Spin_QC_inexact_sample_generator_3()
{
}

template<class Kernel_, class FT_>
template<typename Spin_quadric_iterator_>
void Spin_QC_inexact_sample_generator_3<Kernel_, FT_>::fetch(Spin_quadric_iterator_ quadrics_begin, Spin_quadric_iterator_ quadrics_end)
{
    typedef Kernel_ Kernel;
    typedef FT_ FT;

    // random rotation spin
    Spin_3<FT> random_rotation;
    uniform_random_spin_3(random_rotation);

    // random spin circle
    Random_spin_circle_3<Kernel> random_circle(random_rotation);

    // circle domain cuts
    std::vector<FT> circle_cuts;

    FT solutions[4];
    size_t number_of_solutions, index;

    // generate QC intersections
    for (Spin_quadric_iterator_ iterator = quadrics_begin; iterator != quadrics_end; ++iterator)
    {
        number_of_solutions = random_circle.intersect_quadric(*iterator, solutions);

        for (size_t j = 0; j < number_of_solutions; ++j)
        {
            // check that solution lies on a quadric
            const double MAX_SOLUTION_ERROR = 10-8;
            (void)MAX_SOLUTION_ERROR;

            using Math::fabs;

            assert(fabs(iterator->evaluate(random_circle.evaluate(solutions[j]))) < MAX_SOLUTION_ERROR);
        }

        for (index = 0; index < number_of_solutions; ++index)
            circle_cuts.push_back(solutions[index]);
    }

    // sort solutions
    std::sort(circle_cuts.begin(), circle_cuts.end());

    // consider only unique ones
    typename std::vector<FT>::iterator new_end = std::unique(circle_cuts.begin(), circle_cuts.end());
    circle_cuts.erase(new_end, circle_cuts.end());

    // gather midpoints
    if (circle_cuts.size() < 2)
    {
        CS_logger_debug(MODULE, "Cannot fetch more samples!");
        return;
    }

    FT midarg;

    // add sample which is between the first and last one
    midarg = FT(0.5) * (circle_cuts.front() + circle_cuts.back());

    m_cached_samples.push_back(random_circle.evaluate(midarg));

    // add other samples
    typename std::vector<FT>::const_iterator previous_iterator = circle_cuts.begin();
    typename std::vector<FT>::const_iterator current_iterator = ++circle_cuts.begin();

    while (current_iterator != circle_cuts.end())
    {
        midarg = FT(0.5) * (*previous_iterator + *current_iterator);

        m_cached_samples.push_back(random_circle.evaluate(midarg));

        ++previous_iterator;
        ++current_iterator;
    }

//  CS_logger_debug(MODULE, "Fetched " << m_cached_samples.size() << " samples");
}

template<class Kernel_, class FT_>
template<typename Spin_quadric_iterator>
void Spin_QC_inexact_sample_generator_3<Kernel_, FT_>::operator()(Spin_quadric_iterator quadrics_begin, Spin_quadric_iterator quadrics_end, Sample &sample)
{
    Spin_quadric_iterator next;

    // avoid using QC sample generator for zero or one quadric
    if (quadrics_begin == quadrics_end)
        goto qc_failed;

    next = quadrics_begin;
    ++next;

    if (next == quadrics_end)
        goto qc_failed;

    if (m_cached_samples.empty())
        fetch(quadrics_begin, quadrics_end);

    if (m_cached_samples.empty())
        goto qc_failed;

    sample = m_cached_samples.front();
    m_cached_samples.pop_front();

    return;

qc_failed:
    // return whatever - better luck next time
    Spin_straight_sample_generator_3<Kernel_, FT_> straight_generator;
    straight_generator(quadrics_begin, quadrics_end, sample);
}
} // namespace CS
