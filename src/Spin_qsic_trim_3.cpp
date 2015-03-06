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
#include <cs/Spin_qsic_trim_3.h>
#include <cs/Random.h>
#include <cs/Linear_system.h>

namespace CS
{
Spin_qsic_trim_3::Spin_qsic_trim_3(const Spin_reduced_quadric_3 &spin_reduced_quadric_a, const Spin_reduced_quadric_3 &spin_reduced_quadric_b)
    : m_spin_reduced_quadric_a(spin_reduced_quadric_a),
      m_spin_reduced_quadric_b(spin_reduced_quadric_b)
{
}

bool Spin_qsic_trim_3::begin(double &trim_s12, double &trim_s23, double &trim_s31) const
{
    const int MAXIMUM_NUMBER_OF_INITIAL_GUESSES = 250;
    int number_of_initial_guesses = 0;

    while (++number_of_initial_guesses < MAXIMUM_NUMBER_OF_INITIAL_GUESSES)
    {
        double guess_s12 = random_double(-1, 1);
        double guess_s23 = random_double(-1, 1);
        double guess_s31 = random_double(-1, 1);
        double guess_radius = random_double(0, 1);

        if (!try_refine_begin(guess_s12, guess_s23, guess_s31, guess_radius, trim_s12, trim_s23, trim_s31))
        {
            //std::cerr << "begin FAILED after " << MAXIMUM_NUMBER_OF_INITIAL_GUESSES << std::endl;
        }
        else
        {
            //std::cerr << "begin succeeded after " << number_of_initial_guesses << " initial guesess: " << trim_s12 << ", " << trim_s23 << ", " << trim_s31 << std::endl;
            return true;
        }
    }

    return false;
}

bool Spin_qsic_trim_3::next(double &trim_s12, double &trim_s23, double &trim_s31) const
{
    const double INITIAL_DELTA_STEP = 0.001;
    const double NEXT_DELTA_STEP_FACTOR = 0.5;
    const int MAXIMUM_NUMBER_OF_DELTA_STEP_GUESSES = 4;

    int number_of_delta_guesses = 0;
    double current_delta_step = INITIAL_DELTA_STEP;

    double grad_a_s12, grad_a_s23, grad_a_s31;
    double grad_b_s12, grad_b_s23, grad_b_s31;

    m_spin_reduced_quadric_a.gradient(trim_s12, trim_s23, trim_s31, grad_a_s12, grad_a_s23, grad_a_s31);
    m_spin_reduced_quadric_b.gradient(trim_s12, trim_s23, trim_s31, grad_b_s12, grad_b_s23, grad_b_s31);

    double dir_s12 = grad_a_s23 * grad_b_s31 - grad_a_s31 * grad_b_s23;
    double dir_s23 = grad_a_s31 * grad_b_s12 - grad_a_s12 * grad_b_s31;
    double dir_s31 = grad_a_s12 * grad_b_s23 - grad_a_s23 * grad_b_s12;

    normalize(dir_s12, dir_s23, dir_s31);

    double next_s12, next_s23, next_s31;

    while (++number_of_delta_guesses < MAXIMUM_NUMBER_OF_DELTA_STEP_GUESSES)
    {
        if (!try_refine_next(trim_s12, trim_s23, trim_s31,
                             dir_s12, dir_s23, dir_s31, current_delta_step,
                             next_s12, next_s23, next_s31))
        {
            //std::cerr << "next FAILED for " << current_delta_step << std::endl;

            current_delta_step *= NEXT_DELTA_STEP_FACTOR;
        }
        else
        {
            //std::cerr << "next succeeded for " << current_delta_step << std::endl;

            trim_s12 = next_s12;
            trim_s23 = next_s23;
            trim_s31 = next_s31;

            return true;
        }
    }

    return false;
}

void Spin_qsic_trim_3::normalize(double &p0, double &p1, double &p2)
{
    double len = std::sqrt(p0 * p0 + p1 * p1 + p2 * p2);
    p0 /= len;
    p1 /= len;
    p2 /= len;
}

bool Spin_qsic_trim_3::try_refine_begin_step(double s12, double s23, double s31, double radius,
                                             double &refined_epsilon,
                                             double &refined_s12, double &refined_s23, double &refined_s31) const
{
    double squared_radius = radius * radius;

    double h00, h01, h02;
    double h10, h11, h12;
    double h20, h21, h22;
    double g0, g1, g2;

    m_spin_reduced_quadric_a.gradient(s12, s23, s31, h00, h01, h02);
    m_spin_reduced_quadric_b.gradient(s12, s23, s31, h10, h11, h12);
    h20 = 2 * s12;
    h21 = 2 * s23;
    h22 = 2 * s31;

    g0 = m_spin_reduced_quadric_a.evaluate(s12, s23, s31);
    g1 = m_spin_reduced_quadric_b.evaluate(s12, s23, s31);
    g2 = s12 * s12 + s23 * s23 + s31 * s31 - squared_radius;

    refined_epsilon = std::sqrt(g0 * g0 + g1 * g1 + g2 * g2);

    double x0, x1, x2;

    if (!linear_solve(h00, h01, h02,
                      h10, h11, h12,
                      h20, h21, h22,
                      -g0, -g1, -g2,
                      x0, x1, x2))
    {
        return false;
    }

    refined_s12 = s12 + x0;
    refined_s23 = s23 + x1;
    refined_s31 = s31 + x2;
    return true;
}

bool Spin_qsic_trim_3::try_refine_next_step(double s12, double s23, double s31,
                                            double base_s12, double base_s23, double base_s31,
                                            double dir_s12, double dir_s23, double dir_s31, double delta_step,
                                            double &refined_epsilon,
                                            double &refined_s12, double &refined_s23, double &refined_s31) const
{
    double h00, h01, h02;
    double h10, h11, h12;
    double h20, h21, h22;
    double g0, g1, g2;

    m_spin_reduced_quadric_a.gradient(s12, s23, s31, h00, h01, h02);
    m_spin_reduced_quadric_b.gradient(s12, s23, s31, h10, h11, h12);
    h20 = -dir_s12;
    h21 = -dir_s23;
    h22 = -dir_s31;

    g0 = m_spin_reduced_quadric_a.evaluate(s12, s23, s31);
    g1 = m_spin_reduced_quadric_b.evaluate(s12, s23, s31);
    g2 = (base_s12 - s12) * dir_s12 + (base_s23 - s23) * dir_s23 + (base_s31 - s31) * dir_s31 - delta_step;

    refined_epsilon = std::sqrt(g0 * g0 + g1 * g1 + g2 * g2);

    double x0, x1, x2;

    if (!linear_solve(h00, h01, h02,
                      h10, h11, h12,
                      h20, h21, h22,
                      -g0, -g1, -g2,
                      x0, x1, x2))
    {
        return false;
    }

    refined_s12 = s12 + x0;
    refined_s23 = s23 + x1;
    refined_s31 = s31 + x2;
    return true;
}

bool Spin_qsic_trim_3::try_refine_begin(double s12, double s23, double s31, double radius,
                                        double &refined_s12, double &refined_s23, double &refined_s31) const
{
    double current_s12 = s12, current_s23 = s23, current_s31 = s31;
    double next_s12, next_s23, next_s31;
    double next_epsilon = 10e10;
    const double TARGET_REFINE_EPSILON = 10e-6;
    const int MAXIMUM_NUMBER_OF_REFINE_STEPS = 20;
    int current_number_of_refine_step = 0;

    while (++current_number_of_refine_step < MAXIMUM_NUMBER_OF_REFINE_STEPS
           && try_refine_begin_step(current_s12, current_s23, current_s31, radius,
                                    next_epsilon,
                                    next_s12, next_s23, next_s31))
    {
        current_s12 = next_s12;
        current_s23 = next_s23;
        current_s31 = next_s31;

        // std::cout << "eps: " << neps << std::endl; // FIXME: remove
        // std::cout << "cs12: " << cs12 << ", cs23: " << cs23 << ", cs31: " << cs31 << std::endl; // FIXME: remove

        if (next_epsilon < TARGET_REFINE_EPSILON)
        {
            refined_s12 = current_s12;
            refined_s23 = current_s23;
            refined_s31 = current_s31;

            return true;
        }
    }

    return false;
}

bool Spin_qsic_trim_3::try_refine_next(double s12, double s23, double s31,
                                       double dir_s12, double dir_s23, double dir_s31, double delta_step,
                                       double &refined_s12, double &refined_s23, double &refined_s31) const
{
    double current_s12 = s12, current_s23 = s23, current_s31 = s31;
    double next_s12, next_s23, next_s31;
    double next_epsilon = 10e10;
    const double TARGET_REFINE_EPSILON = 10e-6;
    const int MAXIMUM_NUMBER_OF_REFINE_STEPS = 20;
    int current_number_of_refine_step = 0;

    while (++current_number_of_refine_step < MAXIMUM_NUMBER_OF_REFINE_STEPS
           && try_refine_next_step(current_s12, current_s23, current_s31,
                                   s12, s23, s31,
                                   dir_s12, dir_s23, dir_s31, delta_step,
                                   next_epsilon,
                                   next_s12, next_s23, next_s31))
    {
        current_s12 = next_s12;
        current_s23 = next_s23;
        current_s31 = next_s31;

        // std::cout << "eps: " << neps << std::endl; // FIXME: remove
        // std::cout << "cs12: " << cs12 << ", cs23: " << cs23 << ", cs31: " << cs31 << std::endl; // FIXME: remove

        if (next_epsilon < TARGET_REFINE_EPSILON)
        {
            refined_s12 = current_s12;
            refined_s23 = current_s23;
            refined_s31 = current_s31;

            return true;
        }
    }

    return false;
}
} // namespace CS
