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
#ifndef LIBCS_BENCHMARK_H
#define LIBCS_BENCHMARK_H

#include "Random.h"
#include <functional>
#include <sstream>
#include <vector>
#include <cstdlib>

#ifdef __linux__
#include <sys/time.h>
#endif // __linux__

#if defined(__linux__) || defined(__FreeBSD__) || defined(__OpenBSD__)
#define PRINTF_FORMAT_LLU "%llu"
#endif // defined(__linux__) || defined(__FreeBSD__) || defined(__OpenBSD__)

#ifdef _WIN32
#define PRINTF_FORMAT_LLU "%I64u"
#endif // _WIN32

#ifdef _WIN32
// import by hand, because windows.h messes up everything
#ifndef GetTickCount
extern "C" unsigned long __stdcall GetTickCount();
#endif // GetTickCount
#endif // _WIN32

namespace CS
{
// tick counter
unsigned long long get_tick_count();

// benchmarking
template<class Kernel_>
class Benchmark
{
public:
    typedef typename Kernel_::RT                     RT;
    typedef typename Kernel_::Vector_3               Vector_3;
    typedef typename Kernel_::Plane_3                Plane_3;
    typedef typename Kernel_::Spin_quadric_3         Spin_quadric_3;
    typedef typename Kernel_::Spin_qsic_3            Spin_qsic_3;
    typedef typename Kernel_::Spin_qsip_3            Spin_qsip_3;
    typedef typename Kernel_::Predicate_h_3          Predicate_h_3;
    typedef typename Kernel_::Predicate_s_3          Predicate_s_3;
    typedef typename Kernel_::Predicate_g_3          Predicate_g_3;

    typedef std::function<void (const std::string &)> Report_proc;

public:
    static RT       random_point(int cube_size = 10);
    static Vector_3 random_vector(int cube_size = 10);
    static Plane_3  random_plane(int cube_size = 10);

    void random_H3_intersection_test(size_t count, Report_proc proc) const;
    void intersection_test(const std::vector<Predicate_g_3> &predicates, Report_proc proc) const;
};

// soft/hard benchmark points
class HardPoint
{
private:
    const char *        m_file;
    const char *        m_function;
    int                 m_line;

    // stats
    unsigned long       m_calls;

    unsigned long long  m_min_time;
    unsigned long long  m_max_time;

    unsigned long long  m_total_time;

    friend class SoftPoint;

public:
    HardPoint(const char *file, const char *function, int line);
    ~HardPoint();
};

class SoftPoint
{
private:
    unsigned long long m_start;
    //unsigned long long m_end;
    HardPoint *m_hard_point;

public:
    SoftPoint(HardPoint *hard_point);
    ~SoftPoint();
};
} // namespace CS

//#define FORCE_BENCHMARK

//#if !defined(QT_NO_DEBUG) || defined(FORCE_BENCHMARK)
#if 0
#define CS_BENCHMARK_POINT()                                                \
    static CS::HardPoint _hard_point(__FILE__, __FUNCTION__, __LINE__);     \
    CS::SoftPoint _soft_point(&_hard_point);
#else
#   define CS_BENCHMARK_POINT()
#endif

#include "Benchmark.ipp"

#endif // LIBCS_BENCHMARK_H
