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
#ifndef LIBCS_SPIN_QSIC_3_H
#define LIBCS_SPIN_QSIC_3_H

#include "Benchmark.h"
#include <libqi.h>

namespace CS
{
// Spin_qsic_3:
//    F(u, v) = G dot s / || G ||
//
//    (u, v) in P^1
//    F(u, v) in Spin(3)
//
//    G(u, v) = qsic(P, Q)
//    G(u, v) in R^4
//
//    P, Q - Spin_quadric_3 in P^3
//    s0 - homogenous component
//
//    s = [s12; s23, s31; s0]
//
//    s12^2 + s23^2 + s31^2 + s0^2 = 1
template<class Kernel_>
class Spin_qsic_3
{
    typedef typename Kernel_::RT                 RT;
    typedef typename Kernel_::Vector_3           Vector_3;

    typedef typename Kernel_::Spin_quadric_3     Spin_quadric_3;

    typedef typename Kernel_::Qsic               Qsic;
    typedef typename Kernel_::Qsic_component     Qsic_component;

    typedef typename Kernel_::Spin_qsip_point    Spin_qsip_point;

public:
    typedef Kernel_                              Kernel;

    Spin_qsic_3(const Spin_quadric_3 &q1, const Spin_quadric_3 &q2);

    const Qsic &                qsic() const;
    const Qsic_component &      component(size_t index) const;

    const Spin_quadric_3 &      q1() const;
    const Spin_quadric_3 &      q2() const;

    size_t                      num_components() const;
    bool                        is_component_visible(size_t index) const;
    std::string                 component_to_string(size_t index) const;

    bool                        is_rational() const; // in Q(sqrt(delta))
    bool                        is_smooth() const;

    int                         component_dimension(size_t index) const;

    template<typename OutputIterator>
    void                        extract_self_intersections(OutputIterator output) const
    {
        // this extracts all self-intersections of 1-dimensional components
    }

private:
    Spin_quadric_3      m_q1;
    Spin_quadric_3      m_q2;
    Qsic                m_qsic;
};
} // namespace CS

#include "Spin_qsic_3.ipp"

#endif // LIBCS_SPIN_QSIC_3_H
