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
#include "Loader_sphere_tree.h"

namespace CS
{
template<class Kernel_>
Loader_sphere_tree<Kernel_>::Loader_sphere_tree()
    : m_number_of_levels(0),
      m_level_degree(0),
      m_logger(log4cxx::Logger::getLogger("CS.Loader_sphere_tree"))
{
}

template<class Kernel_>
bool Loader_sphere_tree<Kernel_>::load_from_file(const char *file_name, bool normalize)
{
    std::ifstream file(file_name);

    if (!file.is_open())
        return false;

    std::string line;

    if (!std::getline(file, line))
        return false;

    int number_of_levels, level_degree;

    std::istringstream header(line);

    header >> number_of_levels >> level_degree;

    if (!header)
        return false;

    if (number_of_levels <= 0 || level_degree <= 0)
        return false;

    m_number_of_levels = static_cast<size_t>(number_of_levels);
    m_level_degree = static_cast<size_t>(level_degree);

    int number_of_nodes = 1;
    double max_absolute_coordinate = 0;

    for (int k = 0; k < number_of_levels; k++)
    {
        Level level;

        for (int i = 0; i < number_of_nodes; i++)
        {
            if (!std::getline(file, line))
                return false;

            //  read the sphere
            double x, y, z, r, dummy;
            std::istringstream parser(line);

            parser >> x >> y >> z >> r >> dummy;

            if (parser)
            {
                level.push_back(Ball_3(RT(x), RT(y), RT(z), RT(r)));

                if (fabs(x) > max_absolute_coordinate) max_absolute_coordinate = fabs(x);
                if (fabs(y) > max_absolute_coordinate) max_absolute_coordinate = fabs(y);
                if (fabs(z) > max_absolute_coordinate) max_absolute_coordinate = fabs(z);
            }
        }

        // accumulate level
        m_levels.push_back(level);

        // next degree
        number_of_nodes *= level_degree;
    }

    // normalize if needed
    if (normalize)
    {
        double scale = 1.0 / max_absolute_coordinate;

        for (typename Levels::iterator level_iterator = m_levels.begin(); level_iterator != m_levels.end(); ++level_iterator)
            for (typename Level::iterator ball_iterator = level_iterator->begin(); ball_iterator != level_iterator->end(); ++ball_iterator)
                ball_iterator->scale(scale);
    }

    LOG4CXX_INFO(m_logger, "Loaded sphere tree: " << file_name);
    LOG4CXX_INFO(m_logger, "  Sphere tree: " << m_number_of_levels << " levels");
    LOG4CXX_INFO(m_logger, "  Sphere tree: " << m_level_degree << " degree");

    return true;
}

template<class Kernel_>
size_t Loader_sphere_tree<Kernel_>::number_of_levels() const
{
    return m_number_of_levels;
}

template<class Kernel_>
size_t Loader_sphere_tree<Kernel_>::level_degree() const
{
    return m_level_degree;
}

template<class Kernel_>
typename Loader_sphere_tree<Kernel_>::const_iterator Loader_sphere_tree<Kernel_>::level_begin(size_t level) const
{
    return m_levels[level].begin();
}

template<class Kernel_>
typename Loader_sphere_tree<Kernel_>::const_iterator Loader_sphere_tree<Kernel_>::level_end(size_t level) const
{
    return m_levels[level].end();
}
} // namespace CS
