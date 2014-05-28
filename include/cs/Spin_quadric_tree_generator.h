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
#ifndef LIBCS_SPIN_QUADRIC_TREE_GENERATOR_H
#define LIBCS_SPIN_QUADRIC_TREE_GENERATOR_H

#include <cs/Spin_quadric_tree_3.h>

namespace CS
{
template<class Kernel_>
class Triangle_3
{
    typedef typename Kernel_::Vector_3 Vector_3;

public:
    Triangle_3(const Vector_3 &a,
               const Vector_3 &b,
               const Vector_3 &c)
        : m_a(a),
          m_b(b),
          m_c(c)
    {
    }

    const Vector_3 &a() const
    {
        return m_a;
    }

    const Vector_3 &b() const
    {
        return m_b;
    }

    const Vector_3 &c() const
    {
        return m_c;
    }

private:
    Vector_3 m_a, m_b, m_c;
};

template<class Kernel_>
void spin_quadric_tree_from_tritri(
        typename Kernel_::Spin_quadric_tree_3 &tree, // out
        const Triangle_3<Kernel_> &movable, // movable
        const Triangle_3<Kernel_> &obstacle) // obstacle
{
    typedef typename Kernel_::Spin_quadric_tree_3 Spin_quadric_tree_3;
    Spin_quadric_tree_3 cObstacle;

    (void)movable;
    (void)obstacle;

    tree = tree | cObstacle;
}

// generate quadric tree from scene
template<class Kernel_>
void spin_quadric_tree_from_scene(
        typename Kernel_::Spin_quadric_tree_3 &tree, // out
        const std::vector<Triangle_3<Kernel_> > &movable,
        const std::vector<Triangle_3<Kernel_> > &obstracles)
{
    // add each triangle collision to sentence
    for (typename std::vector<Triangle_3<Kernel_> >::const_iterator i = movable.begin(); i != movable.end(); ++i)
        for (typename std::vector<Triangle_3<Kernel_> >::const_iterator j = obstracles.begin(); j != obstracles.end(); ++j)
            spin_quadric_tree_from_tritri(tree, *i, *j);
}
} // namespace CS

#include "Spin_quadric_tree_generator.ipp"

#endif // LIBCS_SPIN_QUADRIC_TREE_GENERATOR_H
