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
#include "Spin_3.h"

namespace CS
{
template<class FT>
void check_spin_3_norm(const FT &s12, const FT &s23, const FT &s31, const FT &s0)
{
    // exact checking
    assert(s12 * s12 + s23 * s23 + s31 * s31 + s0 * s0 == FT(1));
}

template<class FT_>
Spin_3<FT_>::Spin_3()
    : m_s12(FT(0)),
      m_s23(FT(0)),
      m_s31(FT(0)),
      m_s0(FT(1))
{
    // check
    check_spin_3_norm<FT>(m_s12, m_s23, m_s31, m_s0);
    premultiply();
}

template<class FT_>
Spin_3<FT_>::Spin_3(
        const FT &s12,
        const FT &s23,
        const FT &s31,
        const FT &s0)
    : m_s12(s12),
      m_s23(s23),
      m_s31(s31),
      m_s0(s0)
{
    // check
    check_spin_3_norm(m_s12, m_s23, m_s31, m_s0);
    premultiply();
}

template<class FT_>
bool operator == (const Spin_3<FT_> &lhs, const Spin_3<FT_> &rhs)
{
    return lhs.m_s12 == rhs.m_s12 &&
           lhs.m_s23 == rhs.m_s23 &&
           lhs.m_s31 == rhs.m_s31 &&
           lhs.m_s0 == rhs.m_s0;
}

template<class FT_>
Spin_3<FT_> Spin_3<FT_>::operator-() const
{
    return Spin_3(-m_s12, -m_s23, -m_s31, -m_s0);
}

template<class FT_>
const typename Spin_3<FT_>::FT &Spin_3<FT_>::s12() const
{
    return m_s12;
}

template<class FT_>
const typename Spin_3<FT_>::FT &Spin_3<FT_>::s23() const
{
    return m_s23;
}

template<class FT_>
const typename Spin_3<FT_>::FT &Spin_3<FT_>::s31() const
{
    return m_s31;
}

template<class FT_>
const typename Spin_3<FT_>::FT &Spin_3<FT_>::s0() const
{
    return m_s0;
}

template<class FT_>
const typename Spin_3<FT_>::FT &Spin_3<FT_>::s12s12() const
{
    return m_s12s12;
}

template<class FT_>
const typename Spin_3<FT_>::FT &Spin_3<FT_>::s23s23() const
{
    return m_s23s23;
}

template<class FT_>
const typename Spin_3<FT_>::FT &Spin_3<FT_>::s31s31() const
{
    return m_s31s31;
}

template<class FT_>
const typename Spin_3<FT_>::FT &Spin_3<FT_>::s0s0() const
{
    return m_s0s0;
}

template<class FT_>
const typename Spin_3<FT_>::FT &Spin_3<FT_>::s12s23() const
{
    return m_s12s23;
}

template<class FT_>
const typename Spin_3<FT_>::FT &Spin_3<FT_>::s12s31() const
{
    return m_s12s31;
}

template<class FT_>
const typename Spin_3<FT_>::FT &Spin_3<FT_>::s12s0() const
{
    return m_s12s0;
}

template<class FT_>
const typename Spin_3<FT_>::FT &Spin_3<FT_>::s23s31() const
{
    return m_s23s31;
}

template<class FT_>
const typename Spin_3<FT_>::FT &Spin_3<FT_>::s23s0() const
{
    return m_s23s0;
}

template<class FT_>
const typename Spin_3<FT_>::FT &Spin_3<FT_>::s31s0() const
{
    return m_s31s0;
}

template<class FT>
void Spin_3<FT>::premultiply()
{
    m_s12s12 = m_s12 * m_s12;
    m_s23s23 = m_s23 * m_s23;
    m_s31s31 = m_s31 * m_s31;
    m_s0s0   = m_s0  * m_s0;
    m_s12s23 = m_s12 * m_s23;
    m_s12s31 = m_s12 * m_s31;
    m_s12s0  = m_s12 * m_s0;
    m_s23s31 = m_s23 * m_s31;
    m_s23s0  = m_s23 * m_s0;
    m_s31s0  = m_s31 * m_s0;
}

template<class FT>
FT imaginary_squared_distance(const Spin_3<FT> &lhs, const Spin_3<FT> &rhs)
{
    return (lhs.s12() - rhs.s12()) * (lhs.s12() - rhs.s12()) +
           (lhs.s23() - rhs.s23()) * (lhs.s23() - rhs.s23()) +
           (lhs.s31() - rhs.s31()) * (lhs.s31() - rhs.s31());
}

template<class FT>
Spin_3<FT> slerp(const Spin_3<FT> &spin_0, const Spin_3<FT> &spin_1, double time)
{
    // taken from Shoemake's paper
    FT dot = spin_0.s12() * spin_1.s12() + spin_0.s23() * spin_1.s23() + spin_0.s31() * spin_1.s31() + spin_0.s0()  * spin_1.s0();
    FT angle, sin_angle, sin_time_angle, sin_inv_time_angle, coefficient_0, coefficient_1;

    angle = std::acos(dot);

    if (std::fabs(angle) > 0.0)
    {
        sin_angle           = std::sin(angle);
        sin_time_angle      = std::sin(time * angle);
        sin_inv_time_angle  = std::sin((1.0 - time) * angle);
        coefficient_0       = sin_inv_time_angle / sin_angle;
        coefficient_1       = sin_time_angle / sin_angle;

        double result_s12 = coefficient_0 * spin_0.s12() + coefficient_1 * spin_1.s12();
        double result_s23 = coefficient_0 * spin_0.s23() + coefficient_1 * spin_1.s23();
        double result_s31 = coefficient_0 * spin_0.s31() + coefficient_1 * spin_1.s31();
        double result_s0  = coefficient_0 * spin_0.s0()  + coefficient_1 * spin_1.s0();

        double norm = std::sqrt(result_s12 * result_s12 + result_s23 * result_s23 + result_s31 * result_s31 + result_s0 * result_s0);

        return Spin_3<FT>(result_s12 / norm, result_s23 / norm, result_s31 / norm, result_s0 / norm);
    }
    else
    {
        return spin_0;
    }
}

template<class FT>
std::ostream &operator <<(std::ostream &stream, const Spin_3<FT> &spin)
{
    stream << spin.m_s12 << " * e_12 + "
           << spin.m_s23 << " * e_23 + "
           << spin.m_s31 << " * e_31 + "
           << spin.m_s0 << " * e_0";

    return stream;
}

template<class FT_>
Diff_spin_3<FT_>::Diff_spin_3()
    : m_ds12(FT(0)),
      m_ds23(FT(0)),
      m_ds31(FT(0)),
      m_ds0(FT(0))
{
    // no check
}

template<class FT_>
Diff_spin_3<FT_>::Diff_spin_3(
        const FT &ds12,
        const FT &ds23,
        const FT &ds31,
        const FT &ds0)
    : m_ds12(ds12),
      m_ds23(ds23),
      m_ds31(ds31),
      m_ds0(ds0)
{
    // no check
}

template<class FT_>
bool operator == (const Diff_spin_3<FT_> &lhs, const Diff_spin_3<FT_> &rhs)
{
    return lhs.m_ds12 == rhs.m_ds12 &&
           lhs.m_ds23 == rhs.m_ds23 &&
           lhs.m_ds31 == rhs.m_ds31 &&
           lhs.m_ds0 == rhs.m_ds0;
}

template<class FT_>
Diff_spin_3<FT_> Diff_spin_3<FT_>::operator-() const
{
    return Diff_spin_3(-m_ds12, -m_ds23, -m_ds31, -m_ds0);
}

template<class FT_>
const typename Diff_spin_3<FT_>::FT &Diff_spin_3<FT_>::ds12() const
{
    return m_ds12;
}

template<class FT_>
const typename Diff_spin_3<FT_>::FT &Diff_spin_3<FT_>::ds23() const
{
    return m_ds23;
}

template<class FT_>
const typename Diff_spin_3<FT_>::FT &Diff_spin_3<FT_>::ds31() const
{
    return m_ds31;
}

template<class FT_>
const typename Diff_spin_3<FT_>::FT &Diff_spin_3<FT_>::ds0() const
{
    return m_ds0;
}

template<class FT>
std::ostream &operator <<(std::ostream &stream, const Diff_spin_3<FT> &spin)
{
    stream << spin.m_ds12 << " * d e_12 + "
           << spin.m_ds23 << " * d e_23 + "
           << spin.m_ds31 << " * d e_31 + "
           << spin.m_ds0 << " * d e_0";

    return stream;
}
} // namespace CS
