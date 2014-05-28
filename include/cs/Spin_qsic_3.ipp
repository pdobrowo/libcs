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
#include "Spin_qsic_3.h"

namespace CS
{
template<class Kernel_>
Spin_qsic_3<Kernel_>::Spin_qsic_3(const Spin_quadric_3 &q1, const Spin_quadric_3 &q2)
    : m_q1(q1),
      m_q2(q2)
{
    CS_BENCHMARK_POINT();

    // use QI to intersect
    std::ostringstream logger;
    int opt_level = 0;

    m_qsic = QI::intersection(Kernel::to_lidia_matrix(m_q1.matrix()), Kernel::to_lidia_matrix(m_q2.matrix()), opt_level, logger);

    // spin qsic quadric is: N[m_qsic]
}

template<class Kernel_>
const typename Spin_qsic_3<Kernel_>::Qsic &Spin_qsic_3<Kernel_>::qsic() const
{
    return m_qsic;
}

template<class Kernel_>
const typename Spin_qsic_3<Kernel_>::Spin_quadric_3 &Spin_qsic_3<Kernel_>::q1() const
{
    return m_q1;
}

template<class Kernel_>
const typename Spin_qsic_3<Kernel_>::Spin_quadric_3 &Spin_qsic_3<Kernel_>::q2() const
{
    return m_q2;
}

template<class Kernel_>
const typename Spin_qsic_3<Kernel_>::Qsic_component &Spin_qsic_3<Kernel_>::component(size_t index) const
{
    return m_qsic.cc[index];
}

template<class Kernel_>
size_t Spin_qsic_3<Kernel_>::num_components() const
{
    return static_cast<size_t>(m_qsic.nb_cc);
}

template<class Kernel_>
bool Spin_qsic_3<Kernel_>::is_component_visible(size_t index) const
{
    return m_qsic.cc[index].isInRealAffineSpace();
}

template<class Kernel_>
std::string Spin_qsic_3<Kernel_>::component_to_string(size_t index) const
{
    // output components
    // FIXME: Should we cache it ?
    QIOutputter outputter;

    //outputter.setAffineInputQuadrics();
    //outputter.setAffineParametrizations();
    //outputter.setOmitImaginaryParametrizations();
    //outputter.showInputEuclideanType();
    //outputter.showCutParams();
    outputter.useLaTeX(false);
    outputter.setVerbosityLevel(VERBOSITY_EXHAUSTIVE);
    outputter.disableMultiline();

    outputter.output(m_qsic, m_q1.matrix(), m_q2.matrix());

    // return parametrization
    std::string msg =
        std::string("mode: ") + outputter.getOutput()->parametrizations[index].optimality + "\n" +
        std::string("type: ") + outputter.getOutput()->parametrizations[index].label + "\n" +
        std::string("parametrization: gamma(u, v) = N") + outputter.getOutput()->parametrizations[index].param;

    if (!outputter.getOutput()->parametrizations[index].delta_param.empty())
        msg += "\n" + std::string("delta: ") + outputter.getOutput()->parametrizations[index].delta_param;

    msg += "\n" + std::string("is rational: ") + (is_rational() ? std::string("yes") : std::string("no"));

    return msg;
}

template<class Kernel_>
bool Spin_qsic_3<Kernel_>::is_smooth() const
{
    return m_qsic.ctype == 1 && (m_qsic.rtype == 2 || m_qsic.rtype == 3 || m_qsic.rtype == 4);
}

template<class Kernel_>
bool Spin_qsic_3<Kernel_>::is_rational() const
{
    return !is_smooth();
}

template<class Kernel_>
int Spin_qsic_3<Kernel_>::component_dimension(size_t index) const
{
    switch (m_qsic.cc[index].type)
    {
    case INTER_TYPE_UNDEFINED : return -1;
    case INTER_TYPE_SMOOTH_QUARTIC_BRANCH_1 : return 1;
    case INTER_TYPE_SMOOTH_QUARTIC_BRANCH_2 : return 1;
    case INTER_TYPE_NODAL_QUARTIC : return 1;
    case INTER_TYPE_CUSPIDAL_QUARTIC : return 1;
    case INTER_TYPE_CUBIC : return 1;
    case INTER_TYPE_CONIC : return 1;
    case INTER_TYPE_LINE : return 1;
    case INTER_TYPE_LINES_WITH_CONSTRAINT : return 1;
    case INTER_TYPE_POINT : return 0;
    case INTER_TYPE_SMOOTH_QUADRIC : return 2;
    case INTER_TYPE_CONE : return 2;
    case INTER_TYPE_PAIR_OF_PLANES : return 2;
    case INTER_TYPE_PLANE : return 2;
    case INTER_TYPE_UNIVERSE : return 3;
    default: return -1;
    }
}
} // namespace CS
