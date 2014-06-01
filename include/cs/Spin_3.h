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
#ifndef LIBCS_SPIN_3_H
#define LIBCS_SPIN_3_H

#include <cassert>
#include <ostream>
#include <cmath>

namespace CS
{
template<class FT_>
void check_spin_3_norm(const FT_ &s12, const FT_ &s23, const FT_ &s31, const FT_ &s0);

template<>
void check_spin_3_norm<float>(const float &s12, const float &s23, const float &s31, const float &s0);

template<>
void check_spin_3_norm<double>(const double &s12, const double &s23, const double &s31, const double &s0);

#ifdef Gmpfr // TODO: Fixme
template<>
void check_spin_3_norm<CGAL::Gmpfr>(const CGAL::Gmpfr &s12, const CGAL::Gmpfr &s23, const CGAL::Gmpfr &s31, const CGAL::Gmpfr &s0);
#endif // Gmpfr

template<class FT_>
class Spin_3;

template<class FT_>
bool operator ==(const Spin_3<FT_> &lhs, const Spin_3<FT_> &rhs);

template<class FT_>
std::ostream &operator <<(std::ostream &stream, const Spin_3<FT_> &spin);

template<class FT_>
class Spin_3
{
    typedef FT_                 FT;

public:
    Spin_3();
    Spin_3(const FT &s12,
           const FT &s23,
           const FT &s31,
           const FT &s0);

    friend bool operator == <>(const Spin_3 &lhs, const Spin_3 &rhs);

    Spin_3 operator-() const;

    const FT &s12() const;
    const FT &s23() const;
    const FT &s31() const;
    const FT &s0() const;

    // premultiplied
    const FT &s12s12() const;
    const FT &s23s23() const;
    const FT &s31s31() const;
    const FT &s0s0() const;
    const FT &s12s23() const;
    const FT &s12s31() const;
    const FT &s12s0() const;
    const FT &s23s31() const;
    const FT &s23s0() const;
    const FT &s31s0() const;

    friend std::ostream &operator << <>(std::ostream &stream, const Spin_3<FT> &spin);

private:
    FT m_s12, m_s23, m_s31, m_s0;

    // premultiplied base: for fast quadric sign evaluation
    void premultiply();

    FT m_s12s12;
    FT m_s23s23;
    FT m_s31s31;
    FT m_s0s0;
    FT m_s12s23;
    FT m_s12s31;
    FT m_s12s0;
    FT m_s23s31;
    FT m_s23s0;
    FT m_s31s0;
};

template<class FT>
FT imaginary_squared_distance(const Spin_3<FT> &lhs, const Spin_3<FT> &rhs);

template<class FT>
Spin_3<FT> slerp(const Spin_3<FT> &s1, const Spin_3<FT> &s2, double t);

// spin differential
template<class FT>
class Diff_spin_3;

template<class FT>
bool operator ==(const Diff_spin_3<FT> &lhs, const Diff_spin_3<FT> &rhs);

template<class FT>
std::ostream &operator <<(std::ostream &stream, const Diff_spin_3<FT> &spin);

template<class FT_>
class Diff_spin_3
{
    typedef FT_                 FT;

public:
    Diff_spin_3();
    Diff_spin_3(const FT &ds12,
                const FT &ds23,
                const FT &ds31,
                const FT &ds0);

    friend bool operator == <>(const Diff_spin_3 &lhs, const Diff_spin_3 &rhs);

    Diff_spin_3 operator-() const;

    const FT &ds12() const;
    const FT &ds23() const;
    const FT &ds31() const;
    const FT &ds0() const;

    friend std::ostream &operator << <>(std::ostream &stream, const Diff_spin_3<FT> &spin);

private:
    FT m_ds12, m_ds23, m_ds31, m_ds0;
};
} // namespace CS

#include "Spin_3.ipp"

#endif // LIBCS_SPIN_3_H
