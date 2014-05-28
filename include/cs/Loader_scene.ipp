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
#include "Loader_scene.h"

namespace CS
{
template<class Kernel_, typename VertexOutputIterator_, typename FaceOutputIterator_>
bool load_scene_object(const char *fileName, VertexOutputIterator_ vertexIterator, FaceOutputIterator_ faceIterator,
                       double scale, const CGAL::Vector_3<Kernel_> &translation)
{
    typedef typename Kernel_::RT      RT;
    typedef CGAL::Vector_3<Kernel_>   Vector_3;

    std::ifstream file(fileName);

    if (!file.is_open())
        return false;

    // parse file
    size_t num_vertices;

    if (!(file >> num_vertices))
        return false;

    for (size_t i = 0; i < num_vertices; ++i)
    {
        // data is stored in double format
        double x, y, z;

        if (!(file >> x >> y >> z))
            return false;

        // scale, and convert to target type
        RT sx(scale * x);
        RT sy(scale * y);
        RT sz(scale * z);

        // push back vertex
        (*vertexIterator)++ = Vector_3(sx, sy, sz) + translation;
    }

    size_t num_faces;

    if (!(file >> num_faces))
        return false;

    for (size_t i = 0; i < num_faces; ++i)
    {
        // indexes are stored as unsigned integers
        size_t a, b, c;

        if (!(file >> a >> b >> c))
            return false;

        (*faceIterator)++ = Face(a, b, c);
    }

    return true;
}

template<class Kernel_, typename TriangleOutputIterator_>
bool load_scene(const char *path, TriangleOutputIterator_ outRobotFaces, TriangleOutputIterator_ outObstacleFaces,
                       double scale, const CGAL::Vector_3<Kernel_> &scaledRobotTranslation)
{
    typedef typename Kernel_::Vector_3 Vector_3;
    typedef typename Kernel_::Point_3 Point_3;
    typedef typename Kernel_::Triangle_3 Triangle_3;

    std::string obstacleFileName = std::string(path) + "/" + "obstacle.txt";
    std::string robotFileName = std::string(path) + "/" + "robot.txt";

    std::vector<Vector_3> obstacleVertices;
    std::vector<Face> obstacleFaces;

    std::vector<Vector_3> robotVertices;
    std::vector<Face> robotFaces;

    if (!load_scene_object<Kernel_>(obstacleFileName.c_str(), std::back_inserter(obstacleVertices), std::back_inserter(obstacleFaces), scale) ||
        !load_scene_object<Kernel_>(robotFileName.c_str(), std::back_inserter(robotVertices), std::back_inserter(robotFaces), scale))
    {
        // failed to read one of definition files
        return false;
    }

    // create triangle lists
    for (size_t i = 0; i < robotFaces.size(); ++i)
    {
        (*outRobotFaces)++ = Triangle_3(
            Point_3(robotVertices[robotFaces[i].vertex(0)].x(),
                    robotVertices[robotFaces[i].vertex(0)].y(),
                    robotVertices[robotFaces[i].vertex(0)].z()),
            Point_3(robotVertices[robotFaces[i].vertex(1)].x(),
                    robotVertices[robotFaces[i].vertex(1)].y(),
                    robotVertices[robotFaces[i].vertex(1)].z()),
            Point_3(robotVertices[robotFaces[i].vertex(2)].x(),
                    robotVertices[robotFaces[i].vertex(2)].y(),
                    robotVertices[robotFaces[i].vertex(2)].z()));
    }

    for (size_t i = 0; i < obstacleFaces.size(); ++i)
    {
        // robot object is not translated in order to preserve rotation center
        // a translation is realised by making an inverse translation of scene

        (*outObstacleFaces)++ = Triangle_3(
            Point_3(obstacleVertices[obstacleFaces[i].vertex(0)].x() - scaledRobotTranslation.x(),
                    obstacleVertices[obstacleFaces[i].vertex(0)].y() - scaledRobotTranslation.y(),
                    obstacleVertices[obstacleFaces[i].vertex(0)].z() - scaledRobotTranslation.z()),
            Point_3(obstacleVertices[obstacleFaces[i].vertex(1)].x() - scaledRobotTranslation.x(),
                    obstacleVertices[obstacleFaces[i].vertex(1)].y() - scaledRobotTranslation.y(),
                    obstacleVertices[obstacleFaces[i].vertex(1)].z() - scaledRobotTranslation.z()),
            Point_3(obstacleVertices[obstacleFaces[i].vertex(2)].x() - scaledRobotTranslation.x(),
                    obstacleVertices[obstacleFaces[i].vertex(2)].y() - scaledRobotTranslation.y(),
                    obstacleVertices[obstacleFaces[i].vertex(2)].z() - scaledRobotTranslation.z()));
    }

    return true;
}
} // namespace CS
