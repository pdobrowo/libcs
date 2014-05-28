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
#include "Spin_qsip_3.h"

namespace CS
{
template<class Kernel_>
Spin_qsip_3<Kernel_>::Spin_qsip_3(const Spin_quadric_3 &q, const Spin_qsic_3 &c)
{
    CS_BENCHMARK_POINT();

    // intersect qsic depending on whether the parametrization is rational or not
    if (c.is_smooth())
    {
        // a random intersection is smooth with probability 1
        // this is the general intersetion case
        intersect_smooth_qsic(q, c);
    }
    else
    {
        // rational intersections are all special cases
        // intersect each rational component
        for (size_t index = 0; index < c.num_components(); ++index)
        {
            // discard non-curve components (only 1-dimensional components please)
            if (c.component_dimension(index) != 1)
            {
                // FIXME: Handle 0-dimensional qsics as well
                continue;
            }

            // do intersect
            intersect_rational_component(q, c.component(index));
        }
    }
}

template<class Kernel_>
size_t Spin_qsip_3<Kernel_>::size_of_points() const
{
    return m_points.size();
}

template<class Kernel_>
typename Spin_qsip_3<Kernel_>::Point_const_iterator Spin_qsip_3<Kernel_>::points_begin() const
{
    return m_points.begin();
}

template<class Kernel_>
typename Spin_qsip_3<Kernel_>::Point_const_iterator Spin_qsip_3<Kernel_>::points_end() const
{
    return m_points.end();
}

template<class Kernel_>
const typename Spin_qsip_3<Kernel_>::Spin_qsip_point & Spin_qsip_3<Kernel_>::point_at(size_t index) const
{
    return m_points[index];
}

template<class Kernel_>
void Spin_qsip_3<Kernel_>::mul3(
        Hom_polynomial_with_sqrt &out,
        const Hom_polynomial_with_sqrt &a,
        const Hom_polynomial_with_sqrt &b,
        const Hom_polynomial_with_sqrt &c) const
{
    CS_BENCHMARK_POINT();

    Hom_polynomial_with_sqrt tmp;
    multiply(tmp, a, b);
    multiply(out, tmp, c);
}

template<class Kernel_>
void Spin_qsip_3<Kernel_>::add10(
        Hom_polynomial_with_sqrt &out,
        const Hom_polynomial_with_sqrt a[10]) const
{
    CS_BENCHMARK_POINT();

    Hom_polynomial_with_sqrt b[5];
    add(b[0], a[0], a[1]);
    add(b[1], a[2], a[3]);
    add(b[2], a[4], a[5]);
    add(b[3], a[6], a[7]);
    add(b[4], a[8], a[9]);

    Hom_polynomial_with_sqrt c[2];
    add(c[0], b[0], b[1]);
    add(c[1], b[2], b[3]);

    Hom_polynomial_with_sqrt d[1];
    add(d[0], c[0], c[1]);

    add(out, b[4], d[0]);
}

template<class Kernel_>
typename Spin_qsip_3<Kernel_>::Hom_polynomial_with_sqrt
Spin_qsip_3<Kernel_>::extend_hom_polynomial(const Hom_polynomial &a0)
{
    CS_BENCHMARK_POINT();

    Hom_polynomial_with_sqrt result;
    result.set_degree(a0.degree());

    for (int i = 0; i <= result.degree(); ++i)
        result[i] = Sqrt_extension(a0[i]);

    return result;
}

template<class Kernel_>
typename Spin_qsip_3<Kernel_>::Hom_polynomial_with_sqrt
Spin_qsip_3<Kernel_>::extend_hom_polynomial(
        const Hom_polynomial &a0,
        const Hom_polynomial &a1,
        const RT &root)
{
    CS_BENCHMARK_POINT();

    Hom_polynomial_with_sqrt result;
    result.set_degree(std::max(a0.degree(), a1.degree()));

    for (int i = 0; i <= result.degree(); ++i)
    {
        result[i] = Sqrt_extension(i <= a0.degree() ? a0[i] : RT(0),
                                   i <= a1.degree() ? a1[i] : RT(0),
                                   root);
    }

    return result;
}

template<class Kernel_>
void Spin_qsip_3<Kernel_>::intersect_rational_component(
        const Spin_quadric_3 &quadric,
        const Qsic_component &component)
{
    CS_BENCHMARK_POINT();

    // method:
    //
    // assumes that component is rational
    //
    // spin-qsic = N[gamma(u, v)]
    // quadric = x^T Q x >= 0
    //
    // sols = N[gamma]^T Q N[gamma] >= 0 <=> gamma^T Q gamma >= 0
    // for nodal quartic there may be two parameter solutions for the same point
    Hom_polynomial_with_sqrt c[4];
    Sqrt_extension a11, a22, a33, a44, a12, a13, a14, a23, a24, a34;

    switch (component.nb_cp)
    {
        case 1:
        {
            // parametrization: c1
            Qsic_curve c1 = component.c[0];

            for (size_t i = 0; i<4; ++i)
                c[i] = extend_hom_polynomial(c1[i]);

            a11 = Sqrt_extension(quadric.a11());
            a22 = Sqrt_extension(quadric.a22());
            a33 = Sqrt_extension(quadric.a33());
            a44 = Sqrt_extension(quadric.a44());
            a12 = Sqrt_extension(quadric.a12());
            a13 = Sqrt_extension(quadric.a13());
            a14 = Sqrt_extension(quadric.a14());
            a23 = Sqrt_extension(quadric.a23());
            a24 = Sqrt_extension(quadric.a24());
            a34 = Sqrt_extension(quadric.a34());
        }
        break;

        case 2:
        {        
            assert(component.d[0].degree() == 0);

            // parametrization: c1 + sqrt(D)*c2
            Qsic_curve c1 = component.c[0];
            Qsic_curve c2 = component.c[1];
            RT root = component.d[0][0];

            for (size_t i = 0; i<4; ++i)
                c[i] = extend_hom_polynomial(c1[i], c2[i], root);

            a11 = Sqrt_extension(quadric.a11(), RT(0), root);
            a22 = Sqrt_extension(quadric.a22(), RT(0), root);
            a33 = Sqrt_extension(quadric.a33(), RT(0), root);
            a44 = Sqrt_extension(quadric.a44(), RT(0), root);
            a12 = Sqrt_extension(quadric.a12(), RT(0), root);
            a13 = Sqrt_extension(quadric.a13(), RT(0), root);
            a14 = Sqrt_extension(quadric.a14(), RT(0), root);
            a23 = Sqrt_extension(quadric.a23(), RT(0), root);
            a24 = Sqrt_extension(quadric.a24(), RT(0), root);
            a34 = Sqrt_extension(quadric.a34(), RT(0), root);
        }
        break;

        default:
        {
            // only these may be non-zero: c[0], d[0], c[1]
            assert(component.nb_cp <= 2);
        }
        break;
    }

    // calculate sign
    Hom_polynomial_with_sqrt characteristic;
    Hom_polynomial_with_sqrt coeffs[10];

    mul3(coeffs[0], a11, c[0], c[0]);
    mul3(coeffs[1], a22, c[1], c[1]);
    mul3(coeffs[2], a33, c[2], c[2]);
    mul3(coeffs[3], a44, c[3], c[3]);
    mul3(coeffs[4], 2 * a12, c[0], c[1]);
    mul3(coeffs[5], 2 * a13, c[0], c[2]);
    mul3(coeffs[6], 2 * a14, c[0], c[3]);
    mul3(coeffs[7], 2 * a23, c[1], c[2]);
    mul3(coeffs[8], 2 * a24, c[1], c[3]);
    mul3(coeffs[9], 2 * a34, c[2], c[3]);

    add10(characteristic, coeffs);

    // extract intersections as a zeroes of the characteristic polynomial
    //std::cerr << "rational component characteristic: " << characteristic << std::endl;

    populate_component_intersections(characteristic, component);
}

template<class Kernel_>
void Spin_qsip_3<Kernel_>::populate_component_intersections(
        const Hom_polynomial_with_sqrt &characteristic,
        const Qsic_component &component)
{
    CS_BENCHMARK_POINT();

    std::vector<Hom_root> roots;

    // extract isolated roots
    extract_isolated_hom_roots<Kernel_>(characteristic, std::back_inserter(roots));

    // create qsip points
    for (size_t i = 0; i < roots.size(); ++i)
    {
        // add new points: a conjugate pair
        m_points.push_back(Spin_qsip_point(roots[i], component));
    }
}

template<class Kernel_>
void Spin_qsip_3<Kernel_>::calc_implicit_qsic_equation(
        const bigint_matrix &quadric,           // first quadric
        const Qsic_surface &s1,                 // bilinear parametrization of second quadric: first component of Q[sqrt(xi)]
        const Qsic_surface &s2,                 // bilinear parametrization of second quadric: second component of Q[sqrt(xi)]
        const bigint &delta,                    // bilinear parametrization of second quadric: xi
        Hom_hom_polynomial &out_implicit_s1,    // out hom_hom polynomial: first component of Q[sqrt(xi)]
        Hom_hom_polynomial &out_implicit_s2)    // out hom_hom polynomial: second component of Q[sqrt(xi)]
{
    CS_BENCHMARK_POINT();

    Hom_hom_polynomial tmp;

    // the implicit equation is:
    // s1 q s1 + d s2 q s2 + 2 sqrt(d) s1 q s2
    out_implicit_s1= plug_param_in_quadric(s1, quadric, s1);
    tmp = plug_param_in_quadric(s2, quadric, s2);
    multiply(tmp, tmp, delta);
    add(out_implicit_s1, out_implicit_s1, tmp);

    out_implicit_s2 = plug_param_in_quadric(s1, quadric, s2);
    add(out_implicit_s2, out_implicit_s2, out_implicit_s2); // * 2
}

template<class Kernel_>
void Spin_qsip_3<Kernel_>::intersect_smooth_qsic(const Spin_quadric_3 &quadric, const Spin_qsic_3 &qsic)
{
    CS_BENCHMARK_POINT();

    // note:
    //     parameter D is a constant polynomial (which is xi) or
    //     a polynomial of degree 4. In the second case, the xi is equal to zero.

    // Dupont quadric is: Qr = s1 + sqrt(xi) s2
    Qsic_surface    s1 = qsic.qsic().s1;
    Qsic_surface    s2 = qsic.qsic().s2;
    bigint          d;

    if (qsic.qsic().cc[0].d[0].degree() == 0)
        d = qsic.qsic().cc[0].d[0][0]; // xi is a constant
    else
        d = 0;

    //std::cerr << "Dupont quadric parametrization is:" << std::endl;
    //std::cerr << "  s1 = " << s1 << std::endl;
    //std::cerr << "  s2 = " << s2 << std::endl;
    //std::cerr << "  d = " << d << std::endl;

    //
    // the spin-qsic curve is: C = Qs ^ Qt
    // Dupont quadric is: Qr
    //
    // the spin-qsip point is calculated with the following formula:
    //
    // C ^ Qu = Qs ^ Qt ^ Qu = Qs ^ Qr ^ Qu = (Qr ^ Qs ) ^ (Qr ^ Qu)
    // taken from Hemmer et. al.
    //
    Hom_hom_polynomial fs1, fs2; // implicit f = Qr ^ Qs
    Hom_hom_polynomial gs1, gs2; // implicit g = Qr ^ Qu

    calc_implicit_qsic_equation(Kernel::to_lidia_matrix(qsic.q1().matrix()), s1, s2, d, fs1, fs2);

    if (fs1.is_zero() && fs2.is_zero())
    {
        // ops, the Dupont quadric is the first quadric, try with the other one
        calc_implicit_qsic_equation(Kernel::to_lidia_matrix(qsic.q2().matrix()), s1, s2, d, fs1, fs2);

        if (fs1.is_zero() && fs2.is_zero())
        {
            // something gone wild
            assert(0);
        }
    }

    calc_implicit_qsic_equation(Kernel::to_lidia_matrix(quadric.matrix()), s1, s2, d, gs1, gs2);

    //std::cerr << "implicit curves are:" << std::endl;
    //std::cerr << "  fs1 = " << fs1 << std::endl;
    //std::cerr << "  fs2 = " << fs2 << std::endl;
    //std::cerr << "  and" << std::endl;
    //std::cerr << "  gs1 = " << gs1 << std::endl;
    //std::cerr << "  gs2 = " << gs2 << std::endl;

    // read parameters
    Hom_polynomial_with_sqrt a0, a1, a2, b0, b1, b2;
    Hom_polynomial_with_sqrt a0b1, a1b0, a0b2, a2b0, a1b2, a2b1;
    Hom_polynomial_with_sqrt s01, s02, s12;
    Hom_polynomial_with_sqrt s02s02, s01s12, resultant;
    Hom_polynomial_with_sqrt discriminant;
    Hom_polynomial_with_sqrt a1a1, a0a2;

    a0 = get_coefficient(0, fs1, fs2, d);
    a1 = get_coefficient(1, fs1, fs2, d);
    a2 = get_coefficient(2, fs1, fs2, d);
    b0 = get_coefficient(0, gs1, gs2, d);
    b1 = get_coefficient(1, gs1, gs2, d);
    b2 = get_coefficient(2, gs1, gs2, d);

    multiply(a0b1, a0, b1);
    multiply(a0b2, a0, b2);
    multiply(a1b2, a1, b2);
    multiply(a1b0, a1, b0);
    multiply(a2b0, a2, b0);
    multiply(a2b1, a2, b1);
    subtract(s01, a0b1, a1b0);
    subtract(s02, a0b2, a2b0);
    subtract(s12, a1b2, a2b1);
    multiply(s02s02, s02, s02);
    multiply(s01s12, s01, s12);
    subtract(resultant, s02s02, s01s12);
    multiply(a1a1, a1, a1);
    multiply(a0a2, a0, a2);
    multiply(a0a2, a0a2, d.is_zero() ? Sqrt_extension(4) : Sqrt_extension(4, 0, d));
    subtract(discriminant, a1a1, a0a2);

    //std::cerr << "characteristic (rational resultant) = " << resultant << std::endl;
    //std::cerr << "discriminant (rational resultant) = " << discriminant << std::endl;

    // extract roots
    std::vector<Hom_root> roots;
    extract_isolated_hom_roots<Kernel_>(resultant, std::back_inserter(roots));

    // find out which components are these roots from
    for (size_t i = 0; i < roots.size(); ++i)
    {
        CGAL::Sign sign_discriminant = hom_polynomial_with_sqrt_sign_at<Kernel_>(discriminant, roots[i]);

        // Hemmer: theorem 2
        if (sign_discriminant == CGAL::NEGATIVE)
        {
            // two complex intersection points
            continue;
        }
        else if (sign_discriminant == CGAL::ZERO)
        {
            // a real endpoint of both arcs, take first one
            m_points.push_back(Spin_qsip_point(roots[i], qsic.component(0)));

            //std::cerr << "  extracted point on component 0" << std::endl;
        }
        else
        {
            // otherwise
            CGAL::Sign sign_s01 = hom_polynomial_with_sqrt_sign_at<Kernel_>(s01, roots[i]);

            if (sign_s01 == CGAL::ZERO)
            {
                // two real points, one on each arc
                m_points.push_back(Spin_qsip_point(roots[i], qsic.component(0)));
                m_points.push_back(Spin_qsip_point(roots[i], qsic.component(1)));

                //std::cerr << "  extracted point on component 0 and 1" << std::endl;
            }
            else
            {
                Hom_polynomial_with_sqrt arc_poly;

                Hom_polynomial_with_sqrt tmp1;
                multiply(tmp1, a0, s02);
                multiply(tmp1, tmp1, d.is_zero() ? Sqrt_extension(2) : Sqrt_extension(2, 0, d));

                Hom_polynomial_with_sqrt tmp2;
                multiply(tmp2, a1, s01);

                subtract(arc_poly, tmp1, tmp2);

                // one real intersection point on arc E
                CGAL::Sign sign = -sign_s01 * hom_polynomial_with_sqrt_sign_at<Kernel_>(arc_poly, roots[i]);

                if (sign == CGAL::NEGATIVE)
                {
                    m_points.push_back(Spin_qsip_point(roots[i], qsic.component(0)));

                    //std::cerr << "  extracted point on component 0" << std::endl;
                }
                else
                {
                    m_points.push_back(Spin_qsip_point(roots[i], qsic.component(1)));

                    //std::cerr << "  extracted point on component 1" << std::endl;
                }
            }
        }
    }
}

template<class Kernel_>
typename Spin_qsip_3<Kernel_>::Hom_polynomial_with_sqrt Spin_qsip_3<Kernel_>::get_coefficient(
        int index,
        const Hom_hom_polynomial &a,
        const Hom_hom_polynomial &b,
        const RT &root)
{
    CS_BENCHMARK_POINT();

    Hom_polynomial coeff_a;
    if (a.degree() >= index) coeff_a = a[index];

    if (root.is_zero())
    {
        // extended polynomilal with zero root
        return extend_hom_polynomial(coeff_a);
    }
    else
    {
        // extended polynomial with proper root
        Hom_polynomial coeff_b;
        if (b.degree() >= index) coeff_b = b[index];

        return extend_hom_polynomial(coeff_a, coeff_b, root);
    }
}
} // namespace CS
