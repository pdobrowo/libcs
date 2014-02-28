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
#include "Spin_quadric_3.h"

namespace CS
{
template<class Kernel_>
bool operator <(const Spin_quadric_3<Kernel_> &lhs, const Spin_quadric_3<Kernel_> &rhs)
{
    #define CS_CMP(Coeff) if (lhs.m_##Coeff != rhs.m_##Coeff) return lhs.m_##Coeff < rhs.m_##Coeff;

    CS_CMP(a11);
    CS_CMP(a22);
    CS_CMP(a33);
    CS_CMP(a44);
    CS_CMP(a12);
    CS_CMP(a13);
    CS_CMP(a14);
    CS_CMP(a23);
    CS_CMP(a24);
    CS_CMP(a34);

    #undef CS_CMP

    return false;
}

template<class Kernel_>
Spin_quadric_3<Kernel_>::Spin_quadric_3()
    : m_a11(RT(0)),
      m_a22(RT(0)),
      m_a33(RT(0)),
      m_a44(RT(0)),
      m_a12(RT(0)),
      m_a13(RT(0)),
      m_a14(RT(0)),
      m_a23(RT(0)),
      m_a24(RT(0)),
      m_a34(RT(0))
{
}

template<class Kernel_>
Spin_quadric_3<Kernel_>::Spin_quadric_3(const Predicate_h_3<Kernel_> &h3)
{
    construct(Predicate_g_3<Kernel_>(h3));
}

template<class Kernel_>
Spin_quadric_3<Kernel_>::Spin_quadric_3(const Predicate_s_3<Kernel_> &s3)
{
    construct(Predicate_g_3<Kernel_>(s3));
}

template<class Kernel_>
Spin_quadric_3<Kernel_>::Spin_quadric_3(const Predicate_g_3<Kernel_> &g3)
{
    construct(g3);
}

template<class Kernel_>
const typename Spin_quadric_3<Kernel_>::RT &Spin_quadric_3<Kernel_>::a11() const
{
    return m_a11;
}

template<class Kernel_>
const typename Spin_quadric_3<Kernel_>::RT &Spin_quadric_3<Kernel_>::a22() const
{
    return m_a22;
}

template<class Kernel_>
const typename Spin_quadric_3<Kernel_>::RT &Spin_quadric_3<Kernel_>::a33() const
{
    return m_a33;
}

template<class Kernel_>
const typename Spin_quadric_3<Kernel_>::RT &Spin_quadric_3<Kernel_>::a44() const
{
    return m_a44;
}

template<class Kernel_>
const typename Spin_quadric_3<Kernel_>::RT &Spin_quadric_3<Kernel_>::a12() const
{
    return m_a12;
}

template<class Kernel_>
const typename Spin_quadric_3<Kernel_>::RT &Spin_quadric_3<Kernel_>::a13() const
{
    return m_a13;
}

template<class Kernel_>
const typename Spin_quadric_3<Kernel_>::RT &Spin_quadric_3<Kernel_>::a14() const
{
    return m_a14;
}

template<class Kernel_>
const typename Spin_quadric_3<Kernel_>::RT &Spin_quadric_3<Kernel_>::a23() const
{
    return m_a23;
}

template<class Kernel_>
const typename Spin_quadric_3<Kernel_>::RT &Spin_quadric_3<Kernel_>::a24() const
{
    return m_a24;
}

template<class Kernel_>
const typename Spin_quadric_3<Kernel_>::RT &Spin_quadric_3<Kernel_>::a34() const
{
    return m_a34;
}

template<class Kernel_>
void Spin_quadric_3<Kernel_>::construct(const Predicate_g_3<Kernel_> &g3)
{
#if 0
    // PROBLEM: This implementation may contain a bug (this is not derived from the Mathematica package)

    // parameters
    Vector_3 ab = g3.a() - g3.b();
    Vector_3 kl = g3.k() - g3.l();

    RT abxy = g3.a().x() * g3.b().y() - g3.a().y() * g3.b().x();
    RT abyz = g3.a().y() * g3.b().z() - g3.a().z() * g3.b().y();
    RT abzx = g3.a().z() * g3.b().x() - g3.a().x() * g3.b().z();

    RT klxy = g3.k().x() * g3.l().y() - g3.k().y() * g3.l().x();
    RT klyz = g3.k().y() * g3.l().z() - g3.k().z() * g3.l().y();
    RT klzx = g3.k().z() * g3.l().x() - g3.k().x() * g3.l().z();

    RT kxxy = kl.x() * abxy; RT kxyz = kl.x() * abyz; RT kxzx = kl.x() * abzx;
    RT kyxy = kl.y() * abxy; RT kyyz = kl.y() * abyz; RT kyzx = kl.y() * abzx;
    RT kzxy = kl.z() * abxy; RT kzyz = kl.z() * abyz; RT kzzx = kl.z() * abzx;

    RT axxy = ab.x() * klxy; RT axyz = ab.x() * klyz; RT axzx = ab.x() * klzx;
    RT ayxy = ab.y() * klxy; RT ayyz = ab.y() * klyz; RT ayzx = ab.y() * klzx;
    RT azxy = ab.z() * klxy; RT azyz = ab.z() * klyz; RT azzx = ab.z() * klzx;

    // unreduced quadric
    m_a11 = -kxyz - kyzx + kzxy - axyz - ayzx + azxy; // s12 s12
    m_a22 =  kxyz - kyzx - kzxy + axyz - ayzx - azxy; // s23 s23
    m_a33 = -kxyz + kyzx - kzxy - axyz + ayzx - azxy; // s31 s31
    m_a44 =  kxyz + kyzx + kzxy + azxy + ayzx + axyz; // s0 s0
    m_a12 =  kxxy + kzyz + axxy + azyz; // s12 s23
    m_a13 =  kyxy + kzzx + ayxy + azzx; // s12 s31
    m_a14 =  kxzx - kyyz - axzx + ayyz; // s12 s0
    m_a23 =  kxzx + kyyz + axzx + ayyz; // s23 s31
    m_a24 =  kyxy - kzzx - ayxy + azzx; // s23 s0
    m_a34 = -kxxy + kzyz + axxy - azyz; // s31 s0

#else

    RT kx = g3.k().x();
    RT ky = g3.k().y();
    RT kz = g3.k().z();

    RT lx = g3.l().x();
    RT ly = g3.l().y();
    RT lz = g3.l().z();

    RT ax = g3.a().x();
    RT ay = g3.a().y();
    RT az = g3.a().z();

    RT bx = g3.b().x();
    RT by = g3.b().y();
    RT bz = g3.b().z();

    // TODO: Optimize more
    RT kxmlx = kx - lx;
    RT kymly = ky - ly;
    RT kzmlz = kz - lz;

    RT axmbx = ax - bx;
    RT aymby = ay - by;
    RT azmbz = az - bz;

    m_a11 = (az * by - ay * bz) * kxmlx + (-az * bx + ax * bz) * kymly + azmbz * (-ky * lx + kx * ly) + (-ay * bx + ax * by) * kzmlz - aymby * (kz * lx - kx * lz) - axmbx * (-kz * ly + ky *lz);  // s12^2
    m_a22 = (-az * by + ay * bz) * kxmlx + (-az * bx + ax * bz) * kymly - azmbz * (-ky * lx + kx * ly) + (ay * bx - ax * by) * kzmlz  - aymby * (kz * lx - kx * lz) + axmbx * (-kz * ly + ky * lz); // s23^2
    m_a33 = (az * by - ay * bz) * kxmlx + (az * bx - ax * bz) * kymly  - azmbz * (-ky * lx + kx * ly) + (ay * bx - ax * by) * kzmlz  + aymby * (kz * lx - kx * lz) - axmbx * (-kz * ly + ky * lz); // s31^2
    m_a44 = (-az * by + ay * bz) * kxmlx + (az * bx - ax * bz) * kymly  + azmbz * (-ky * lx + kx * ly) + (-ay * bx + ax * by) * kzmlz + aymby * (kz * lx - kx * lz) + axmbx * (-kz * ly + ky * lz); // s0^2
    m_a12 = (-ay * bx + ax * by) * kxmlx + axmbx * (-ky * lx + kx * ly) + (-az * by + ay * bz) * kzmlz + azmbz * (-kz * ly + ky * lz); // s23 s12
    m_a13 = (-ay * bx + ax * by) * kymly + aymby * (-ky * lx + kx * ly) + (az * bx - ax * bz) * kzmlz  + azmbz * (kz * lx - kx * lz);  // s31 s12
    m_a14 = (az * bx - ax * bz) * kxmlx - (-az * by + ay * bz) * kymly - axmbx * (kz * lx - kx * lz)  + aymby * (-kz * ly + ky * lz); // s0 s12
    m_a23 = (az * bx - ax * bz) * kxmlx + (-az * by + ay * bz) * kymly + axmbx * (kz * lx - kx * lz)  + aymby * (-kz * ly + ky * lz); // s23 s31
    m_a24 = (-ay * bx + ax * by) * kymly - aymby * (-ky * lx + kx * ly) - (az * bx - ax * bz) * kzmlz  + azmbz * (kz * lx - kx * lz);  // s23 s0
    m_a34 = (ay * bx - ax * by) * kxmlx + axmbx * (-ky * lx + kx * ly) + (-az * by + ay * bz) * kzmlz - azmbz * (-kz * ly + ky * lz); // s0 s31

#endif

    // reduced quadric
    RT a55 = g3.c();

    m_a11 += a55;
    m_a22 += a55;
    m_a33 += a55;
    m_a44 += a55;
}

template<class Kernel_>
typename Spin_quadric_3<Kernel_>::Matrix Spin_quadric_3<Kernel_>::matrix() const
{
    // TODO: This may be cached instead of raw values
    Matrix m(4, 4);

    m.sto(0, 0, m_a11);
    m.sto(0, 1, m_a12);
    m.sto(0, 2, m_a13);
    m.sto(0, 3, m_a14);
    m.sto(1, 0, m_a12);
    m.sto(1, 1, m_a22);
    m.sto(1, 2, m_a23);
    m.sto(1, 3, m_a24);
    m.sto(2, 0, m_a13);
    m.sto(2, 1, m_a23);
    m.sto(2, 2, m_a33);
    m.sto(2, 3, m_a34);
    m.sto(3, 0, m_a14);
    m.sto(3, 1, m_a24);
    m.sto(3, 2, m_a34);
    m.sto(3, 3, m_a44);

    return m;
}

template<class Kernel_>
typename Spin_quadric_3<Kernel_>::Matrix Spin_quadric_3<Kernel_>::ellipsoid_matrix() const
{
    // TODO: This may be cached instead of raw values
    Matrix m(4, 4);

    m.sto(0, 0, m_a11 + 1);
    m.sto(0, 1, m_a12);
    m.sto(0, 2, m_a13);
    m.sto(0, 3, m_a14);
    m.sto(1, 0, m_a12);
    m.sto(1, 1, m_a22 + 1);
    m.sto(1, 2, m_a23);
    m.sto(1, 3, m_a24);
    m.sto(2, 0, m_a13);
    m.sto(2, 1, m_a23);
    m.sto(2, 2, m_a33 + 1);
    m.sto(2, 3, m_a34);
    m.sto(3, 0, m_a14);
    m.sto(3, 1, m_a24);
    m.sto(3, 2, m_a34);
    m.sto(3, 3, m_a44 + 1);

    return m;
}

template<class Kernel_>
std::ostream &operator <<(std::ostream &os, const Spin_quadric_3<Kernel_> &q)
{
    return (os << "[" << q.a11() << ";" << q.a22() << ";" << q.a33() << ";" << q.a44() << ";"
                      << q.a12() << ";" << q.a13() << ";" << q.a14() << ";" << q.a23() << ";"
                      << q.a24() << ";" << q.a34() << "]");
}

template<class Kernel_>
std::string Spin_quadric_3<Kernel_>::to_string() const
{
    //return QI::quad2string(matrix(), true);

    // copied and adopted from libqi: kernel/QIElem.cc
    bool homogeneous = true;
    RT v[10] = { m_a11, 2 * m_a12, 2 * m_a13, 2 * m_a14, m_a22, 2 * m_a23, 2 * m_a24, m_a33, 2 * m_a34, m_a44 };

    const char *hom_labels[] = {"x^2","x*y","x*z","x*w","y^2","y*z","y*w","z^2","z*w","w^2"};
    const char *aff_labels[] = {"x^2","x*y","x*z","x","y^2","y*z","y","z^2","z"};
    const char **labels = homogeneous ? hom_labels : aff_labels;

    bool first_coefficient = true;
    std::stringstream result(std::stringstream::in | std::stringstream::out);
    RT coefficient;

    for (int k = 0; k < 10; k ++)
    {
        if (v[k] == 0)
            continue;

        if (!first_coefficient && v[k] > 0)
            result << "+";
        else
            if (v[k] < 0)
                result << "-";

        coefficient = abs(v[k]);

        if (coefficient != 1 || (k > 8 && !homogeneous))
            result << coefficient;

        if (homogeneous || (k <= 8))
        {
            if (coefficient != 1)
                result << "*";

            result << labels[k];
        }

        first_coefficient = false;
    }

    std::string result_string = result.str();
    return result_string.length() > 0 ? result_string : "0";

}

template<class Kernel_>
template<class NT>
NT Spin_quadric_3<Kernel_>::evaluate(const CS::Spin_3<NT> &spin) const
{
#if 0
    // old code
    return m_a11 * spin.s12() * spin.s12() +
           m_a22 * spin.s23() * spin.s23() +
           m_a33 * spin.s31() * spin.s31() +
           m_a44 * spin.s0() * spin.s0() +
           NT(2) * (m_a12 * spin.s12() * spin.s23() +
                    m_a13 * spin.s12() * spin.s31() +
                    m_a14 * spin.s12() * spin.s0() +
                    m_a23 * spin.s23() * spin.s31() +
                    m_a24 * spin.s23() * spin.s0() +
                    m_a34 * spin.s31() * spin.s0());
#else
    // use premultiplied spin base
    return m_a11 * spin.s12s12() +
           m_a22 * spin.s23s23() +
           m_a33 * spin.s31s31() +
           m_a44 * spin.s0s0() +
           NT(2) * (m_a12 * spin.s12s23() +
                    m_a13 * spin.s12s31() +
                    m_a14 * spin.s12s0() +
                    m_a23 * spin.s23s31() +
                    m_a24 * spin.s23s0() +
                    m_a34 * spin.s31s0());
#endif
}

template<class Kernel_>
template<class NT>
CGAL::Sign Spin_quadric_3<Kernel_>::evaluate_sign(const Spin_3<NT> &spin) const
{
#if 1

    return CGAL::sign(evaluate(spin));

#else

    RT value = evaluate(spin);
    std::cout << "evaluate() = " << value << std::endl;
    return CGAL::sign(value);

#endif
}

template<class Kernel_>
typename Spin_quadric_3<Kernel_>::RT Spin_quadric_3<Kernel_>::max_abs_coefficient() const
{
    RT result = RT(0);

    result = std::max(result, std::max(m_a11, -m_a11));
    result = std::max(result, std::max(m_a22, -m_a22));
    result = std::max(result, std::max(m_a33, -m_a33));
    result = std::max(result, std::max(m_a44, -m_a44));
    result = std::max(result, std::max(m_a12, -m_a12));
    result = std::max(result, std::max(m_a13, -m_a13));
    result = std::max(result, std::max(m_a14, -m_a14));
    result = std::max(result, std::max(m_a23, -m_a23));
    result = std::max(result, std::max(m_a24, -m_a24));
    result = std::max(result, std::max(m_a34, -m_a34));

    return result;
}

template<class Kernel_>
void Spin_quadric_3<Kernel_>::scale(const RT &s)
{
    assert(s > RT(0));

    m_a11 *= s;
    m_a22 *= s;
    m_a33 *= s;
    m_a44 *= s;
    m_a12 *= s;
    m_a13 *= s;
    m_a14 *= s;
    m_a23 *= s;
    m_a24 *= s;
    m_a34 *= s;
}

template<class Kernel_>
typename Spin_quadric_3<Kernel_>::RT max_abs_difference(const Spin_quadric_3<Kernel_> &lhs, const Spin_quadric_3<Kernel_> &rhs)
{
    typename Spin_quadric_3<Kernel_>::RT result =
        typename Spin_quadric_3<Kernel_>::RT(0);

    result = std::max(result, std::max(lhs.a11() - rhs.a11(), rhs.a11() - lhs.a11()));
    result = std::max(result, std::max(lhs.a22() - rhs.a22(), rhs.a22() - lhs.a22()));
    result = std::max(result, std::max(lhs.a33() - rhs.a33(), rhs.a33() - lhs.a33()));
    result = std::max(result, std::max(lhs.a44() - rhs.a44(), rhs.a44() - lhs.a44()));
    result = std::max(result, std::max(lhs.a12() - rhs.a12(), rhs.a12() - lhs.a12()));
    result = std::max(result, std::max(lhs.a13() - rhs.a13(), rhs.a13() - lhs.a13()));
    result = std::max(result, std::max(lhs.a14() - rhs.a14(), rhs.a14() - lhs.a14()));
    result = std::max(result, std::max(lhs.a23() - rhs.a23(), rhs.a23() - lhs.a23()));
    result = std::max(result, std::max(lhs.a24() - rhs.a24(), rhs.a24() - lhs.a24()));
    result = std::max(result, std::max(lhs.a34() - rhs.a34(), rhs.a34() - lhs.a34()));

    return result;
}

template<class Kernel_>
void Spin_quadric_3<Kernel_>::inverse()
{
    m_a11 = -m_a11;
    m_a22 = -m_a22;
    m_a33 = -m_a33;
    m_a44 = -m_a44;
    m_a12 = -m_a12;
    m_a13 = -m_a13;
    m_a14 = -m_a14;
    m_a23 = -m_a23;
    m_a24 = -m_a24;
    m_a34 = -m_a34;
}

template<class Kernel_>
Spin_quadric_3<Kernel_> Spin_quadric_3<Kernel_>::inversed() const
{
    Spin_quadric_3<Kernel_> result = *this;
    result.inverse();
    return result;
}
} // namespace CS
