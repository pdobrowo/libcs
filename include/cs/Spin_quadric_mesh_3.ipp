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
#include "Spin_quadric_mesh_3.h"

namespace CS
{
/*
 * VTK stub code
 *
ImplicitSphere.cxx

#include <vtkSmartPointer.h>

#include <vtkSampleFunction.h>
#include <vtkContourFilter.h>
#include <vtkOutlineFilter.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkImageData.h>

#include <vtkSphere.h>

int main (int argc, char *argv[])
{
  vtkSmartPointer<vtkSphere> sphere = vtkSmartPointer<vtkSphere>::New();

  // Sample the function
  vtkSmartPointer<vtkSampleFunction> sample = vtkSmartPointer<vtkSampleFunction>::New();
  sample->SetSampleDimensions(50,50,50);
  sample->SetImplicitFunction(sphere);
  double value = 2.0;
  double xmin = -value, xmax = value,
         ymin = -value, ymax = value,
         zmin = -value, zmax = value;
  sample->SetModelBounds(xmin, xmax, ymin, ymax, zmin, zmax);

  // Create the 0 isosurface
  vtkSmartPointer<vtkContourFilter> contours = vtkSmartPointer<vtkContourFilter>::New();
  contours->SetInputConnection(sample->GetOutputPort());
  contours->GenerateValues(1, 1, 1);

  // Map the contours to graphical primitives
  vtkSmartPointer<vtkPolyDataMapper> contourMapper =  vtkSmartPointer<vtkPolyDataMapper>::New();
  contourMapper->SetInputConnection(contours->GetOutputPort());
  contourMapper->ScalarVisibilityOff();

  // Create an actor for the contours
  vtkSmartPointer<vtkActor> contourActor = vtkSmartPointer<vtkActor>::New();
  contourActor->SetMapper(contourMapper);

  // Visualize
  vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
  vtkSmartPointer<vtkRenderWindow> renderWindow = vtkSmartPointer<vtkRenderWindow>::New();
  renderWindow->AddRenderer(renderer);
  vtkSmartPointer<vtkRenderWindowInteractor> interactor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
  interactor->SetRenderWindow(renderWindow);

  renderer->AddActor(contourActor);
  renderer->SetBackground(.2, .3, .4);

  renderWindow->Render();
  interactor->Start();

  return EXIT_SUCCESS;
}
*/
template<class Kernel_>
typename Spin_quadric_mesh_3<Kernel_>::FT Spin_quadric_mesh_3<Kernel_>::evaluate_with_sign(Point_3 p, FT sign)
{
    FT squared_length = p.x() * p.x() + p.y() * p.y() + p.z() * p.z();

    if (squared_length > FT(1))
        return FT(1); // out of bounds (any positive value)

    FT pw = sign * static_cast<FT>(std::sqrt(1.0 - static_cast<double>(squared_length)));

    // minus changes the orientation of the spin-surface so that enclosed is the forbidden part
    // CGAL encloses negative part, not positive as defined by a predicate
    return -(m_a11 * p.x() * p.x() +
             m_a22 * p.y() * p.y() +
             m_a33 * p.z() * p.z() +
             m_a44 * pw * pw +
             FT(2) * (m_a12 * p.x() * p.y() +
                      m_a13 * p.x() * p.z() +
                      m_a14 * p.x() * pw +
                      m_a23 * p.y() * p.z() +
                      m_a24 * p.y() * pw +
                      m_a34 * p.z() * pw));
}

template<class Kernel_>
typename Spin_quadric_mesh_3<Kernel_>::FT Spin_quadric_mesh_3<Kernel_>::evaluate_left(Point_3 p)
{
    return evaluate_with_sign(p, FT(-1));
}

template<class Kernel_>
typename Spin_quadric_mesh_3<Kernel_>::FT Spin_quadric_mesh_3<Kernel_>::evaluate_right(Point_3 p)
{
    return evaluate_with_sign(p, FT(1));
}

template<class Kernel_>
typename Spin_quadric_mesh_3<Kernel_>::Vector_3 Spin_quadric_mesh_3<Kernel_>::evaluate_gradient(Point_3 p, FT sign)
{
    FT x = p.x();
    FT y = p.y();
    FT z = p.z();
    FT squared_length = x * x + y * y + z * z;

    // fix floating arithmetic near space border at radius 1
    if (squared_length >= FT(0.8))
    {
        // essential point, use theoretical limit instead

        // return FT(2) * Vector_3(m_a11 * x + m_a12 * y + m_a13 * z,
        //                         m_a22 * y + m_a12 * x + m_a23 * z,
        //                         m_a33 * z + m_a13 * x + m_a23 * y);
        return FT(2) * Vector_3(x, y, z);
    }
    else
    {
        // standard computations: w component is non-zero
        FT w = static_cast<FT>(std::sqrt(1.0 - static_cast<double>(squared_length)));

        return FT(2) * Vector_3(x * (m_a11 - m_a44) + m_a12 * y + m_a13 * z + m_a14 * sign * (w - x * x / w),
                                y * (m_a22 - m_a44) + m_a12 * x + m_a23 * z + m_a24 * sign * (w - y * y / w),
                                z * (m_a33 - m_a44) + m_a13 * x + m_a23 * y + m_a34 * sign * (w - z * z / w));
    }
}

template<class Kernel_>
template<typename OutputIterator>
void Spin_quadric_mesh_3<Kernel_>::mesh_triangle_soup(
        OutputIterator output_iterator_left,
        OutputIterator output_iterator_right,
        double angular_bound,
        double radius_bound,
        double distance_bound)
{
    Surface_mesh_triangulation_3 triangulation_left, triangulation_right;
    Complex_2_in_triangulation_3 c2t3_left(triangulation_left);
    Complex_2_in_triangulation_3 c2t3_right(triangulation_right);

    // mesh to a triangulation
    mesh_internal_complex_in_triangulation(c2t3_left, c2t3_right, angular_bound, radius_bound, distance_bound);

    // convert to a well oriented triangle soup
    mesh_internal_triangle_soup(output_iterator_left, triangulation_left, c2t3_left, FT(-1));
    mesh_internal_triangle_soup(output_iterator_right, triangulation_right, c2t3_right, FT(1));
}

template<class Kernel_>
template<typename OutputIterator>
void Spin_quadric_mesh_3<Kernel_>::mesh_internal_triangle_soup(
        OutputIterator output_iterator,
        Surface_mesh_triangulation_3 &triangulation,
        Complex_2_in_triangulation_3 &c2t3,
        const FT &sign)
{
    // convert to a well oriented triangle soup
    typedef typename Surface_mesh_triangulation_3::Vertex_handle Vertex_handle;
    typedef typename Surface_mesh_triangulation_3::Point Point;
    typedef typename Surface_mesh_triangulation_3::Geom_traits::Vector_3 Vector;
    typedef typename Surface_mesh_triangulation_3::Edge Edge;
    typedef typename Surface_mesh_triangulation_3::Facet Facet;
    typedef typename Surface_mesh_triangulation_3::Finite_vertices_iterator Finite_vertices_iterator;
    typedef typename Surface_mesh_triangulation_3::Finite_facets_iterator Finite_facets_iterator;

    const typename Surface_mesh_triangulation_3::size_type number_of_facets = c2t3.number_of_facets();

    // Finite vertices coordinates.
    std::map<Vertex_handle, int> vertex_to_index;
    std::vector<Point_3> points;
    int index = 0;

    for (Finite_vertices_iterator vertex_iterator = triangulation.finite_vertices_begin();
         vertex_iterator != triangulation.finite_vertices_end(); ++vertex_iterator)
    {
        vertex_to_index[vertex_iterator] = index++;
        points.push_back(static_cast<Point>(vertex_iterator->point()));
    }

    Finite_facets_iterator facet_iterator = triangulation.finite_facets_begin();
    std::set<Facet> oriented_set;
    std::stack<Facet> stack;

    CGAL_assertion_code(typename Surface_mesh_triangulation_3::size_type nb_facets = 0;)

    while (oriented_set.size() != number_of_facets)
    {
        while (facet_iterator->first->is_facet_on_surface(facet_iterator->second) == false ||
               oriented_set.find(*facet_iterator) != oriented_set.end() ||
               oriented_set.find(c2t3.opposite_facet(*facet_iterator)) != oriented_set.end() )
        {
            ++facet_iterator;
        }

        oriented_set.insert(*facet_iterator);

        stack.push(*facet_iterator);

        while(! stack.empty() )
        {
            Facet f = stack.top();
            stack.pop();

            for (int ih = 0 ; ih < 3 ; ++ih)
            {
                int i1  = triangulation.vertex_triple_index(f.second, triangulation.cw(ih));
                int i2  = triangulation.vertex_triple_index(f.second, triangulation.ccw(ih));

                if (c2t3.face_status(Edge(f.first, i1, i2)) == Complex_2_in_triangulation_3::REGULAR)
                {
                    Facet fn = c2t3.neighbor(f, ih);

                    if (oriented_set.find(fn) == oriented_set.end() &&
                        oriented_set.find(c2t3.opposite_facet(fn)) == oriented_set.end())
                    {
                        oriented_set.insert(fn);
                        stack.push(fn);
                    }
                }
            }
        }
    }

    // Orients the whole mesh towards outside
    // - find the facet with max z
    typename std::set<Facet>::const_iterator top_facet = oriented_set.begin();

    for (typename std::set<Facet>::const_iterator fit = oriented_set.begin(); fit != oriented_set.end(); ++fit)
    {
        double top_z = (top_facet->first->vertex(triangulation.vertex_triple_index(top_facet->second, 0))->point().z()
                      + top_facet->first->vertex(triangulation.vertex_triple_index(top_facet->second, 1))->point().z()
                      + top_facet->first->vertex(triangulation.vertex_triple_index(top_facet->second, 2))->point().z())/3.;
        double z = (fit->first->vertex(triangulation.vertex_triple_index(fit->second, 0))->point().z()
                  + fit->first->vertex(triangulation.vertex_triple_index(fit->second, 1))->point().z()
                  + fit->first->vertex(triangulation.vertex_triple_index(fit->second, 2))->point().z())/3.;

        if (top_z < z)
            top_facet = fit;
    }

    // - orient the facet with max z towards +Z axis
    Vertex_handle v0 = top_facet->first->vertex(triangulation.vertex_triple_index(top_facet->second, 0));
    Vertex_handle v1 = top_facet->first->vertex(triangulation.vertex_triple_index(top_facet->second, 1));
    Vertex_handle v2 = top_facet->first->vertex(triangulation.vertex_triple_index(top_facet->second, 2));
    Vector normal = cross_product(v1->point()-v0->point(), v2->point()-v1->point());

    const Vector z_vector(0, 0, 1);
    bool regular_orientation = (z_vector * normal >= 0);

    for (typename std::set<Facet>::const_iterator fit = oriented_set.begin(); fit != oriented_set.end(); ++fit)
    {
        int indices[3];
        int index = 0;

        for (int i = 0; i < 3; i++)
            indices[index++] = vertex_to_index[fit->first->vertex(triangulation.vertex_triple_index(fit->second, i))];

        Triangle_3 triangle(points[indices[0]],
                            regular_orientation ? points[indices[1]] : points[indices[2]],
                            regular_orientation ? points[indices[2]] : points[indices[1]]);

        Vector_3 gradients[3] =
        {
            evaluate_gradient(triangle[0], sign),
            evaluate_gradient(triangle[1], sign),
            evaluate_gradient(triangle[2], sign)
        };

        (*output_iterator)++ = Smooth_triangle_3(triangle,
                                                 gradients[0],
                                                 gradients[1],
                                                 gradients[2]);

        CGAL_assertion_code(++nb_facets);
    }

    CGAL_assertion(nb_facets == number_of_facets);
}

template<class Kernel_>
Spin_quadric_mesh_3<Kernel_>::Spin_quadric_mesh_3(const Spin_quadric_3<Kernel_> &spin_quadric)
    : m_spin_quadric(spin_quadric)
{
    // doublify quadric coefficients
    m_a11 = to_double(m_spin_quadric.a11());
    m_a12 = to_double(m_spin_quadric.a12());
    m_a13 = to_double(m_spin_quadric.a13());
    m_a14 = to_double(m_spin_quadric.a14());
    m_a22 = to_double(m_spin_quadric.a22());
    m_a23 = to_double(m_spin_quadric.a23());
    m_a24 = to_double(m_spin_quadric.a24());
    m_a33 = to_double(m_spin_quadric.a33());
    m_a34 = to_double(m_spin_quadric.a34());
    m_a44 = to_double(m_spin_quadric.a44());
}

template<class Kernel_>
void Spin_quadric_mesh_3<Kernel_>::mesh_polyhedron(
        Polyhedron_3 &polyhedron_left,
        Polyhedron_3 &polyhedron_right,
        double angular_bound,
        double radius_bound,
        double distance_bound)
{
    Surface_mesh_triangulation_3 triangulation_left, triangulation_right;
    Complex_2_in_triangulation_3 c2t3_left(triangulation_left);
    Complex_2_in_triangulation_3 c2t3_right(triangulation_right);

    // mesh to a triangulation
    mesh_surface_triangulation(c2t3_left, c2t3_right, angular_bound, radius_bound, distance_bound);

    // convert to poyhedron
    CGAL::facets_in_complex_2_to_triangle_mesh(c2t3_left, polyhedron_left);
    CGAL::facets_in_complex_2_to_triangle_mesh(c2t3_right, polyhedron_right);
}

template<class Kernel_>
void Spin_quadric_mesh_3<Kernel_>::mesh_surface_triangulation(
        Surface_mesh_triangulation_3 &triangulation_left,
        Surface_mesh_triangulation_3 &triangulation_right,
        double angular_bound,
        double radius_bound,
        double distance_bound)
{
    Complex_2_in_triangulation_3 c2t3_left(triangulation_left);
    Complex_2_in_triangulation_3 c2t3_right(triangulation_right);

    mesh_internal_complex_in_triangulation(c2t3_left, c2t3_right, angular_bound, radius_bound, distance_bound);
}

template<class Kernel_>
void Spin_quadric_mesh_3<Kernel_>::mesh_internal_complex_in_triangulation(
        Complex_2_in_triangulation_3 &c2t3_left,
        Complex_2_in_triangulation_3 &c2t3_right,
        double angular_bound,
        double radius_bound,
        double distance_bound)
{
    Point_3 leftOrigin = CGAL::ORIGIN;
    Point_3 rightOrigin = CGAL::ORIGIN;

    // Implicit surface
    using std::placeholders::_1;

    Implicit_surface_3 surface_left(
        std::bind(&Spin_quadric_mesh_3<Kernel_>::evaluate_left, this, _1),
        Sphere_3(leftOrigin, 4.)); // 4 - square of maximum sphere radius for correct evaluation

    Implicit_surface_3 surface_right(
        std::bind(&Spin_quadric_mesh_3<Kernel_>::evaluate_right, this, _1),
        Sphere_3(rightOrigin, 4.));

    // meshing criteria
    CGAL::Surface_mesh_default_criteria_3<Surface_mesh_triangulation_3>
        criteria_left(angular_bound, radius_bound, distance_bound);

    CGAL::Surface_mesh_default_criteria_3<Surface_mesh_triangulation_3>
        criteria_right(angular_bound, radius_bound, distance_bound);

    // mesh surfaces
    CGAL::make_surface_mesh(c2t3_left, surface_left, criteria_left, CGAL::Non_manifold_tag());
    CGAL::make_surface_mesh(c2t3_right, surface_right, criteria_right, CGAL::Non_manifold_tag());
}
} // namespace CS
