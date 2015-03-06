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
#ifndef LIBCS_PREDICATE_BB_3_H
#define LIBCS_PREDICATE_BB_3_H

#include <ostream>

namespace CS
{
// BB3: ball-ball collision predicate
//
//      A, r - rotating ball with center in point A and radius r
//      B, l - stable ball with center in point B and radius l
//
// predicate evaluates to true if and only if
// there is a collision of two balls
//
template<class Kernel_>
class Predicate_bb_3
{
    typedef typename Kernel_::RT            RT;
    typedef typename Kernel_::FT            FT;
    typedef typename Kernel_::Vector_3      Vector_3;
    typedef typename Kernel_::Plane_3       Plane_3;
    typedef typename Kernel_::Ball_3        Ball_3;
    typedef typename Kernel_::Predicate_h_3 Predicate_h_3;

public:
    typedef Kernel_                         Kernel;

    Predicate_bb_3();

    Predicate_bb_3(
            const Vector_3 &B,
            const RT &l,
            const Vector_3 &A,
            const RT &r);

    Predicate_bb_3(
            const Ball_3 &b,
            const Ball_3 &a);

    const Vector_3 &    b() const;
    const RT &          l() const;
    const Vector_3 &    a() const;
    const RT &          r() const;

    // access to sub-predicates
    //
    // return value is a pointer to a table with one element
    static const size_t SUB_PREDICATE_COUNT = 1;
    typedef Predicate_h_3 Sub_predicate;

    const Predicate_h_3 *sub_predicates() const;
    bool evaluate(const bool *signs) const;

private:
    Vector_3    m_b;
    RT          m_l;
    Vector_3    m_a;
    RT          m_r;

    Predicate_h_3 m_predicate;

    void construct(
            const Vector_3 &B,
            const RT &l,
            const Vector_3 &A,
            const RT &r);

};
} // namespace CS

#include "Predicate_bb_3.ipp"

#endif // LIBCS_PREDICATE_BB_3_H
