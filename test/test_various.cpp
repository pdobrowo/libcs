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
#include <cs/Predicate_g_3.h>
#include <cs/Spin_inexact_kernel_3.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Gmpq.h>
#include <log4cxx/logger.h>

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("Test.test_various"));

void test_various()
{
    typedef CGAL::Gmpq RT;
    typedef CGAL::Simple_cartesian<RT> Base_kernel;
    typedef CS::Spin_inexact_kernel_3<Base_kernel> Kernel;
    typedef Kernel::Predicate_g_3 Predicate_g_3;
    typedef Kernel::Predicate_h_3 Predicate_h_3;
    typedef Kernel::Vector_3 Vector_3;
    typedef Kernel::Plane_3 Plane_3;
    typedef Kernel::Point_3 Point_3;
    typedef Predicate_g_3::Parametrization Parametrization;

    Predicate_h_3 h(Vector_3(1, 2, 3),
                    Plane_3(Point_3(0, 0, 0), Vector_3(1, 0, 0)));

    Predicate_g_3 g(h);

    Parametrization p = g.parametrization();

    CS::Spin_3<double> v = p.approx_evaluate(0, 0);

}
