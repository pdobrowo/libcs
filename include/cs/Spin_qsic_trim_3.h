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
#ifndef LIBCS_SPIN_QSIC_TRIM_3_H
#define LIBCS_SPIN_QSIC_TRIM_3_H

#include "Spin_inexact_kernel_3.h"
#include "Spin_reduced_quadric_3.h"

namespace CS
{
// Netwon iterations
class Spin_qsic_trim_3
{
public:
    typedef CS::Spin_reduced_quadric_3<Default_inexact_kernel> Spin_reduced_quadric_3;

    Spin_qsic_trim_3(const Spin_reduced_quadric_3 &spin_reduced_quadric_a, const Spin_reduced_quadric_3 &spin_reduced_quadric_b);

    bool begin(double &trim_s12, double &trim_s23, double &trim_s31) const;
    bool next(double &trim_s12, double &trim_s23, double &trim_s31) const;

private:
    const Spin_reduced_quadric_3 &m_spin_reduced_quadric_a;
    const Spin_reduced_quadric_3 &m_spin_reduced_quadric_b;

    static void normalize(double &p0, double &p1, double &p2);

    bool try_refine_begin_step(double s12, double s23, double s31, double radius,
                               double &refined_epsilon,
                               double &refined_s12, double &refined_s23, double &refined_s31) const;

    bool try_refine_next_step(double s12, double s23, double s31,
                              double base_s12, double base_s23, double base_s31,
                              double dir_s12, double dir_s23, double dir_s31, double delta_step,
                              double &refined_epsilon,
                              double &refined_s12, double &refined_s23, double &refined_s31) const;

    bool try_refine_begin(double s12, double s23, double s31, double radius,
                          double &refined_s12, double &refined_s23, double &refined_s31) const;

    bool try_refine_next(double s12, double s23, double s31,
                         double dir_s12, double dir_s23, double dir_s31, double delta_step,
                         double &refined_s12, double &refined_s23, double &refined_s31) const;

};
} // namespace CS

#endif // LIBCS_SPIN_QSIC_TRIM_3_H
