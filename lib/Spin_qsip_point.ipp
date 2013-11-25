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
#include "Spin_qsip_point.h"

namespace CS
{
template<class Kernel_>
Spin_qsip_point<Kernel_>::Spin_qsip_point(const Hom_root &root, const Qsic_component &component)
    : m_root(root),
      m_component(component)
{
}

template<class Kernel_>
const typename Spin_qsip_point<Kernel_>::Hom_root &Spin_qsip_point<Kernel_>::root() const
{
    return m_root;
}

template<class Kernel_>
const typename Spin_qsip_point<Kernel_>::Qsic_component &Spin_qsip_point<Kernel_>::component() const
{
    return m_component;
}
} // namespace CS
