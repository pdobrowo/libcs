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
#ifndef LIBCS_LOADER_SCENE_H
#define LIBCS_LOADER_SCENE_H

#include <CGAL/Vector_3.h>
#include <vector>
#include <string>
#include <fstream>

namespace CS
{
// Spin_configuration_space_3:
//
class Face
{
public:
    Face(size_t a, size_t b, size_t c);

    size_t vertex(size_t i) const;

private:
    size_t m_vertex[3];
};

template<class Kernel_, typename VertexOutputIterator_, typename FaceOutputIterator_>
bool load_scene_object(const char *fileName, VertexOutputIterator_ vertexIterator, FaceOutputIterator_ faceIterator,
                       double scale = 1.0, const CGAL::Vector_3<Kernel_> &translation = CGAL::Vector_3<Kernel_>(0, 0, 0));

template<class Kernel_, typename TriangleOutputIterator_>
bool load_scene(const char *path, TriangleOutputIterator_ outRobotFaces, TriangleOutputIterator_ outObstacleFaces,
                double scale = 1.0, const CGAL::Vector_3<Kernel_> &scaledRobotTranslation = CGAL::Vector_3<Kernel_>(0, 0, 0));
} // namespace CS

#include "Loader_scene.ipp"

#endif // LIBCS_LOADER_SCENE_H
