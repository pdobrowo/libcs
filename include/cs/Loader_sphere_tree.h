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
#ifndef LIBCS_LOADER_SPHERE_TREE_H
#define LIBCS_LOADER_SPHERE_TREE_H

#include "Ball_3.h"
#include <cs/Logger.h>
#include <vector>
#include <cstddef>
#include <fstream>

namespace CS
{
template<class Kernel_>
class Loader_sphere_tree
{
    typedef typename Kernel_::Ball_3        Ball_3;

    typedef std::vector<Ball_3>             Level;
    typedef std::vector<Level>              Levels;

    typedef typename Kernel_::RT            RT;

    // module
    static constexpr char const * const MODULE = "CS.Loader_sphere_tree";

public:
    typedef typename Level::const_iterator  const_iterator;

    Loader_sphere_tree();

    bool load_from_file(const char *file_name, bool normalize = true);

    size_t number_of_levels() const;
    size_t level_degree() const;

    const_iterator level_begin(size_t level) const;
    const_iterator level_end(size_t level) const;

private:
    size_t      m_number_of_levels;
    size_t      m_level_degree;
    Levels      m_levels;
};
} // namespace CS

#include "Loader_sphere_tree.ipp"

#endif // LIBCS_LOADER_SPHERE_TREE_H
