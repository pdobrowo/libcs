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
#include "Benchmark.h"

namespace CS
{
template<class Kernel_>
typename Benchmark<Kernel_>::RT Benchmark<Kernel_>::random_point(int cube_size)
{
    return random_int(-cube_size, cube_size);
}

template<class Kernel_>
typename Benchmark<Kernel_>::Vector_3 Benchmark<Kernel_>::random_vector(int cube_size)
{
    // take vector from cube: [0; cube_size]^3
    int x = random_int(-cube_size, cube_size);
    int y = random_int(-cube_size, cube_size);
    int z = random_int(-cube_size, cube_size);

    return Vector_3(x, y, z);
}

template<class Kernel_>
typename Benchmark<Kernel_>::Plane_3 Benchmark<Kernel_>::random_plane(int cube_size)
{
    // take plane coefficients from cube: [0; cube_size]^3
    int a = random_int(-cube_size, cube_size);
    int b = random_int(-cube_size, cube_size);
    int c = random_int(-cube_size, cube_size);
    int d = random_int(-cube_size, cube_size);

    return Plane_3(a, b, c, d);
}

template<class Kernel_>
void Benchmark<Kernel_>::random_H3_intersection_test(size_t count, Report_proc proc) const
{
    // get predicates
    unsigned long long checkPredicatesBegin = get_tick_count();
    unsigned long long now;

    std::vector<Predicate_g_3> predicates;

    for (size_t i = 0; i < count; ++i)
        predicates.push_back(Predicate_h_3(random_vector(), random_plane()));

    now = get_tick_count();
    {
        std::ostringstream out;
        out << predicates.size() << " predicates created in " << now - checkPredicatesBegin << " ms; "
            << "average " << (now - checkPredicatesBegin) / predicates.size() << " ms per predicate";
        proc(out.str());
    }

    // default intersection test
    intersection_test(predicates, proc);

    // cleanup
    predicates.clear();
}

template<class Kernel_>
void Benchmark<Kernel_>::intersection_test(const std::vector<Predicate_g_3> &predicates, Report_proc proc) const
{
    unsigned long long now;

    // create spin quadrics
    unsigned long long checkQuadricsBegin = get_tick_count();

    std::vector<Spin_quadric_3 *> quadrics;

    for (size_t i = 0; i < predicates.size(); ++i)
        quadrics.push_back(new Spin_quadric_3(predicates[i]));

    now = get_tick_count();
    {
        std::ostringstream out;
        out << quadrics.size() << " quadrics created in " << now - checkQuadricsBegin << " ms; "
            << "average " << (now - checkQuadricsBegin) / quadrics.size() << " ms per quadric";
        proc(out.str());
    }

    // create qsics
    unsigned long long checkQsicsBegin = get_tick_count();

    std::vector<Spin_qsic_3 *> qsics;

    for (size_t i = 0; i < quadrics.size(); ++i)
        for (size_t j = 0; j < quadrics.size(); ++j)
            qsics.push_back(new Spin_qsic_3(*quadrics[i], *quadrics[j]));

    now = get_tick_count();
    {
        std::ostringstream out;
        out << qsics.size() << " qsics created in " << now - checkQsicsBegin << " ms; "
            << "average " << (now - checkQsicsBegin) / qsics.size() << " ms per qsic";
        proc(out.str());
    }

    // create qsips
    unsigned long long checkQsipsBegin = get_tick_count();

    std::vector<Spin_qsip_3 *> qsips;

    for (size_t i = 0; i < quadrics.size(); ++i) // Quadric i
        for (size_t j = 0; j < qsics.size(); ++j) // Qsic j
            qsips.push_back(new Spin_qsip_3(*quadrics[i], *qsics[j]));

    now = get_tick_count();
    {
        std::ostringstream out;
        out << qsips.size() << " qsips created in " << now - checkQsipsBegin << " ms; "
            << "average " << (now - checkQsipsBegin) / qsips.size() << " ms per qsip";
        proc(out.str());
    }

    // cleanup
    for (typename std::vector<Spin_quadric_3 *>::iterator iterator = quadrics.begin(); iterator != quadrics.end(); ++iterator)
        delete *iterator;
    quadrics.clear();

    for (typename std::vector<Spin_qsic_3 *>::iterator iterator = qsics.begin(); iterator != qsics.end(); ++iterator)
        delete *iterator;
    qsics.clear();

    for (typename std::vector<Spin_qsip_3 *>::iterator iterator = qsips.begin(); iterator != qsips.end(); ++iterator)
        delete *iterator;
    qsips.clear();
}
} // namespace CS
