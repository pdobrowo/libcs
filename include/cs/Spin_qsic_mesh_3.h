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
#ifndef LIBCS_SPIN_QSIC_MESH_3
#define LIBCS_SPIN_QSIC_MESH_3

#include "Spin_3.h"
#include <CGAL/Kernel/global_functions.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <iterator>
#include <deque>
#include <cmath>
#include <lidia/bigfloat.h>
#include <lidia/math_vector.h>
#include <lidia/bigint.h>

namespace CS
{
// Spin_qsic_mesh_3:
//     mesher
//
// This is a decorator. It does not take ownership of the object.
template<class Kernel_>
class Spin_qsic_mesh_3
{
    typedef typename Kernel_::RT                 RT;
    typedef typename Kernel_::Qsic_component     Qsic_component;
    typedef typename Kernel_::Qsic_curve         Qsic_curve;

    typedef typename Kernel_::Hom_polynomial     Hom_polynomial;

public:
    typedef CS::Spin_3<double> Spin_3;

    typedef std::deque<Spin_3> Spin_list_3;

private:
    typedef LiDIA::bigfloat                     bigfloat;
    typedef bigfloat                            bigfloat_vector[4];
    typedef LiDIA::math_vector<LiDIA::bigint>   bigint_vector;

    class scoped_bigfloat_precision
    {
    public:
        scoped_bigfloat_precision(long precision)
        {
            m_old_precision = bigfloat::get_precision();
            bigfloat::set_precision(precision);
        }

        ~scoped_bigfloat_precision()
        {
            bigfloat::set_precision(m_old_precision);
        }

    private:
        long m_old_precision;
    };

    static void fromVb(bigfloat_vector &out, const bigint_vector &v)
    {
        for (size_t i = 0; i < 4; ++i)
            out[i] = v[i];
    }

    static void add(bigfloat_vector &out, const bigfloat_vector &a, const bigfloat_vector &b)
    {
        for (size_t i = 0; i < 4; ++i)
            out[i] = a[i] + b[i];
    }

    static void mul(bigfloat_vector &out, const bigfloat_vector &a, const bigfloat &m)
    {
        for (size_t i = 0; i < 4; ++i)
            out[i] = a[i] * m;
    }

    static Spin_3 projectQsic(const bigfloat_vector &gamma)
    {
        bigfloat norm = LiDIA::sqrt(gamma[0] * gamma[0] + gamma[1] * gamma[1] + gamma[2] * gamma[2] + gamma[3] * gamma[3]);
        bigfloat_vector normedGamma;

        if (norm == 0)
            return Spin_3();

        for (size_t i = 0; i < 4; ++i)
            normedGamma[i] = gamma[i] / norm;

        double x, y, z, w;
        normedGamma[0].doublify(x); // e12
        normedGamma[1].doublify(y); // e23
        normedGamma[2].doublify(z); // e31
        normedGamma[3].doublify(w); // 1

        return Spin_3(x, y, z, w);
    }

    static bool evaluate_spin_qsic_component(bigfloat_vector &result, const Qsic_component *component, const bigint &a, const bigint &b)
    {
        // The structure used to store each component of the intersection curve
        // In the worst case, the format of the component is:
        //
        //   c[0] + sqrt(d[0])*c[1] + (c[2] + sqrt(d[0])*c[3])*sqrt(d[1] + sqrt(d[0])*d[2]),
        //
        // where the cp's are curve_params, d[0] is a bigint, and d[1], d[2] are hom_polys
        //
        // - type: type of real component (1: smooth quartic, 2: nodal quartic,
        //     3: cuspidal quartic, 4: cubic, 5: conic, 6: line, 7: lines
        //     with constraint (four concurrent lines), 8: point, 9: smooth
        //     quadric, 10: cone, 11: pair of planes, 12: plane, 13: universe
        // - nb_cp: number of cp's in the representation of the parameterization
        // - cp's: the cp[] in the representation above
        // - d's: the d[] in the representation above
        // - m's: matrices when the output is a surface
        // - h: used for bihomogeneous equation when smooth quartic
        Qsic_curve cpC1(4), cpC2(4), cpC3(4), cpC4(4);
        Hom_polynomial cpD, cpD1, cpD2;

        if (component->nb_cp >= 1)
        {
            cpC1 = component->c[0];
        }

        if (component->nb_cp >= 2)
        {
            cpC2 = component->c[1];
            cpD = component->d[0];
        }

        if (component->nb_cp >= 3)
        {
            cpC3 = component->c[2];
            cpD1 = component->d[1];
        }

        if (component->nb_cp >= 4)
        {
            cpC4 = component->c[3];
            cpD2 = component->d[2];
        }

        // parametrization: c1 + sqrt(D) * c2 + (c3 + sqrt(D) * c4) * sqrt(D1 + D2 * sqrt(D))
        bigint_vector c1    = cpC1.eval(a, b);
        bigint_vector c2    = cpC2.eval(a, b);
        bigint_vector c3    = cpC3.eval(a, b);
        bigint_vector c4    = cpC4.eval(a, b);
        bigint D            = cpD.eval(a, b);
        bigint D1           = cpD1.eval(a, b);
        bigint D2           = cpD2.eval(a, b);

        // increase precision (the default is 5 digits only)
        scoped_bigfloat_precision prec(128);

        bigfloat_vector fc1, fc2, fc3, fc4;

        fromVb(fc1, c1);
        fromVb(fc2, c2);
        fromVb(fc3, c3);
        fromVb(fc4, c4);

        bigfloat fD(D);
        bigfloat fD1(D1);
        bigfloat fD2(D2);

        // out of bounds
        if (fD.is_lt_zero())
            return false;

        bigfloat fDsqrt = LiDIA::sqrt(fD);

        bigfloat inner = fD1 + fD2 * fDsqrt;

        // out of bounds
        if (inner.is_lt_zero())
            return false;

        inner = LiDIA::sqrt(inner);

        bigfloat_vector right;

        mul(right, fc4, fDsqrt);
        add(right, right, fc3);

        bigfloat_vector left;

        mul(left, fc2, fDsqrt);
        add(left, left, fc1);

        // final
        mul(result, right, inner);
        add(result, result, left);

        return true;
    }

public:
    typedef Kernel_                      Kernel;

    Spin_qsic_mesh_3(const Spin_qsic_3<Kernel> &spinQsic);

    size_t          size_of_components() const;
    bool            evaluate_component(Spin_3 &outSpin, size_t component, const RT &a, const RT &b) const;

    static bool     evaluate_component(Spin_3 &outSpin, const Qsic_component &component, const RT &a, const RT &b);

    void            mesh_component(Spin_list_3 &outCurve, size_t component, double radiusBound);

    bool            is_component_visible(size_t component) const;

private:
    const Spin_qsic_3<Kernel>        &m_spinQsic;

    template<typename OutputIterator>
    bool            mesh_component_walk(OutputIterator outputIterator,  // output curve points
                                        size_t component,               // which component is being traced
                                        double radiusBound,             // meshing target maximum radius bound
                                        const bigint &minimumStep,      // 1 / minimumStep allowed to execute, otherwise give up
                                        const bigint &start_a,          // parameter to begin: numerator
                                        const bigint &start_b,          // parameter to begin: denominator
                                        bool positive_direction)        // which parametrization direction take
    {
//#define SHOW_WALKING
        Spin_qsic_mesh_3::Spin_3 currentSpin, lastSpin;

#ifdef SHOW_WALKING
        //std::cerr << "qsic spin: begin" << std::endl;
#endif

        // accepted parameters
        bigint a = start_a;
        bigint b = start_b;

        // next jump: start with 1/2
        bigint da = 0; // 0 indicates a value below 1
        bigint db = 2;

        // next parameter try
        bigint na;
        bigint nb;

        // initial position
        evaluate_component(lastSpin, component, a, b);

        Spin_3 veryFirstSpin = lastSpin;
        bool hasVeryFirstSpin = true;

        //
        // note: the parametrization has the nice property that we will never go round and round, because
        //       the parametrization is on P^1 which is not homeomorphic to S^1
        //
        for (;;)
        {
            // new jump
            if (da == 0)
            {
                // na / nb = a / b + 1 / db = (a db + b) / (b db)
                if (positive_direction)
                    na = a * db + b;
                else
                    na = a * db - b;

                nb = b * db;
            }
            else
            {
                // na / nb = a / b + da / db = (a db + b da) / (b db)
                if (positive_direction)
                    na = a * db + b * da;
                else
                    na = a * db - b * da;

                nb = b * db;
            }

#ifdef SHOW_WALKING
            //std::cerr << "_check: na / nb: " << na << " / " << nb << std::endl;
#endif

            // evaluate the curve at the selected parameter
            bool black_hole = false;

            if (!evaluate_component(currentSpin, component, na, nb))
            {
                // jumped into a hole, treat it as a failure or abandon
                if (hasVeryFirstSpin)
                {
                    // we did not start in a good place - a black hole is a bad place to start
                    return false;
                }
                else
                {
                    // we just ended up in a black hole - try to get closer to it
                    black_hole = true;
                }
            }

            // check bound
            double bound;

            if (!black_hole)
            {
                // calculate bound if we did not get into a blackhole
                bound = std::sqrt(imaginary_squared_distance(lastSpin, currentSpin));
            }

            // did we succedd ?
            if (!black_hole &&
                bound <= radiusBound)
            {
                // add next spin point
                if (hasVeryFirstSpin)
                {
                    // this is added with a delay because, if we totally fail, we do not want to mess output
                    // iterator with a single spin
                    *(outputIterator++) = veryFirstSpin;
                    hasVeryFirstSpin = false;
                }

                *(outputIterator++) = currentSpin;

                // we might have gone very good
                // ...too good
                if (lastSpin == currentSpin)
                    return true;

                // save last
                lastSpin = currentSpin;
                a = na;
                b = nb;

                // correct jump
                if (da == 0)
                {
                    db /= 2;

                    if (db == 1)
                        da = 1;
                }
                else
                {
                    da *= 2;
                }

#ifdef SHOW_WALKING
                //std::cerr << "qsic spin (" << a << " / " << b << " ) : " << currentSpin.s12 << " e12 + " << currentSpin.s23 << " e23 + " << currentSpin.s31 << " e31 + " << currentSpin.s0 << std::endl;
#endif
            }
            else
            {
                // correct jump
                if (da == 0)
                {
                    db *= 2;

                    if (db >= minimumStep)
                        return true; // BREAK
                }
                else
                {
                    da /= 2;

                    if (da == 0)
                    {
                        db = 2;
                    }
                }
            }

#ifdef SHOW_WALKING
            //std::cerr << "_new_jump: da / db: " << da << " / " << db << std::endl;
#endif
        }
    }

    void mesh_component_optimize(Spin_list_3 &outCurve, const Spin_list_3 &dirtyCurve);
};

template<class Kernel_>
Spin_qsic_mesh_3<Kernel_>::Spin_qsic_mesh_3(const Spin_qsic_3<Kernel> &spinQsic)
    : m_spinQsic(spinQsic)
{
}

template<class Kernel_>
size_t Spin_qsic_mesh_3<Kernel_>::size_of_components() const
{
    // forward
    return m_spinQsic.num_components();
}

template<class Kernel_>
bool Spin_qsic_mesh_3<Kernel_>::is_component_visible(size_t component) const
{
    // forward
    return m_spinQsic.is_component_visible(component);
}

template<class Kernel_>
bool Spin_qsic_mesh_3<Kernel_>::evaluate_component(Spin_3 &outSpin, size_t component, const RT &a, const RT &b) const
{
    bigfloat_vector qsicValue;

    // evaluate QSIC and project onto spin space
    if (evaluate_spin_qsic_component(qsicValue, m_spinQsic.qsic().cc + component, a, b))
    {
        // argument in component domain
        outSpin = projectQsic(qsicValue);
        return true;
    }

    // out of domain
    return false;
}

template<class Kernel_>
bool Spin_qsic_mesh_3<Kernel_>::evaluate_component(Spin_3 &outSpin, const Qsic_component &component, const RT &a, const RT &b)
{
    bigfloat_vector qsicValue;

    // evaluate QSIC and project onto spin space
    if (evaluate_spin_qsic_component(qsicValue, &component, a, b))
    {
        // argument in component domain
        outSpin = projectQsic(qsicValue);
        return true;
    }

    // out of domain
    return false;
}

template<class Kernel_>
void Spin_qsic_mesh_3<Kernel_>::mesh_component(Spin_list_3 &outCurve, size_t component, double radiusBound)
{
    // bounds
    const bigint minimumStep = 1 << 24; // maximum number of subdivisions

    // initial parameter: try not to hit any essential point
    bigint initial_parameter_a, initial_parameter_b;
    Spin_list_3 dirtyCurve;
    bool direction;

    // try to init parametrization near zero
    initial_parameter_a = 0;
    initial_parameter_b = 1;
    direction = true;

    if (!mesh_component_walk(std::back_inserter(dirtyCurve), component, radiusBound, minimumStep, initial_parameter_a, initial_parameter_b, direction))
    {
        // now, try to init parametrization near minus infinity
        initial_parameter_a = -1000000;
        initial_parameter_b = 1;
        direction = true;

        if (!mesh_component_walk(std::back_inserter(dirtyCurve), component, radiusBound, minimumStep, initial_parameter_a, initial_parameter_b, direction))
        {
            // now, try to init parametrization near plus infinity
            initial_parameter_a = 1000000;
            initial_parameter_b = 1;
            direction = false;

            if (!mesh_component_walk(std::back_inserter(dirtyCurve), component, radiusBound, minimumStep, initial_parameter_a, initial_parameter_b, direction))
            {
                // failed everywhere - resign
                return;
            }
        }
    }

    // go negative from the starting point
    mesh_component_walk(std::front_inserter(dirtyCurve), component, radiusBound, minimumStep, initial_parameter_a, initial_parameter_b, !direction);

    // cleanup dity curve - remove duplicated nodes
    mesh_component_optimize(outCurve, dirtyCurve);
}

template<class Kernel_>
void Spin_qsic_mesh_3<Kernel_>::mesh_component_optimize(Spin_list_3 &outCurve, const Spin_list_3 &dirtyCurve)
{
    // cleanup dity curve - remove duplicated nodes
    outCurve.clear();

    if (dirtyCurve.empty())
        return;

    outCurve.push_back(dirtyCurve.front());

    typename Spin_list_3::const_iterator curr, prev, last;
    last = dirtyCurve.end();
    prev = dirtyCurve.begin();
    curr = prev;
    ++curr;

    const double MIN_DIST = 0.000001;

    while (curr != last)
    {
        if (std::sqrt(imaginary_squared_distance(*curr, *prev)) <= MIN_DIST)
        {
            // it is too near, ignore this node
            ++curr;
            continue;
        }

        outCurve.push_back(*curr);
        prev = curr;
        ++curr;
    }
}
} // namespace CS

#endif // LIBCS_SPIN_QSIC_MESH_3
