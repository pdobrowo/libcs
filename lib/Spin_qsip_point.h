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
#ifndef LIBCS_SPIN_QSIP_POINT_H
#define LIBCS_SPIN_QSIP_POINT_H

namespace CS
{
template<class Kernel_>
class Spin_qsip_point
{
public:
    typedef Kernel_                              R;
    typedef typename Kernel_::Spin_qsic_3        Spin_qsic_3;
    typedef typename Kernel_::Qsic_component     Qsic_component;

    typedef typename Kernel_::Hom_root           Hom_root;

    Spin_qsip_point(const Hom_root &root, const Qsic_component &component);

    const Hom_root &        root() const;
    const Qsic_component &  component() const;

private:
    Hom_root            m_root;
    Qsic_component      m_component;
};
} // namespace CS

#include "Spin_qsip_point.ipp"

#endif // LIBCS_SPIN_QSIP_POINT_H
