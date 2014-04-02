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
#include <cs/Benchmark.h>
#include <cstdio>
#if defined(__linux__) || defined(__FreeBSD__) || defined(__OpenBSD__)
#include <sys/time.h>
#endif // defined(__linux__) || defined(__FreeBSD__) || defined(__OpenBSD__)

namespace CS
{
int random_int(int min, int max)
{
    return static_cast<int>(
                static_cast<long long>(min) +
                    static_cast<long long>(rand()) * static_cast<long long>(max - min) / static_cast<long long>(RAND_MAX));
}

double random_double(double min, double max)
{
    return min + (max - min) * (static_cast<double>(rand()) / static_cast<double>(RAND_MAX));
}

unsigned long long get_tick_count()
{    
#if defined(__linux__) || defined(__FreeBSD__) || defined(__OpenBSD__)
    struct timeval tv;
    gettimeofday(&tv, 0);
    return static_cast<unsigned long long>(tv.tv_sec) * 1000 + static_cast<unsigned long long>(tv.tv_usec) / 1000;
#endif // defined(__linux__) || defined(__FreeBSD__) || defined(__OpenBSD__)

#ifdef _WIN32
    // GetTickCount declared as extern in header
    return static_cast<unsigned long long>(GetTickCount());
#endif // _WIN32
}

HardPoint::HardPoint(const char *file, const char *function, int line)
    : m_file(file),
      m_function(function),
      m_line(line),
      m_calls(0),
      m_min_time(static_cast<unsigned long long>(-1)),
      m_max_time(0),
      m_total_time(0)

{
}

HardPoint::~HardPoint()
{
    const char *file = strrchr(m_file, '/');

    if (file)
        ++file;
    else
        file = m_file;

    // print stats
    fprintf(stderr, "@ HardPoint: %s:%i [%s]:\n", file, m_line, m_function);
    fprintf(stderr, "    Total calls: %lu\n", m_calls);
    fprintf(stderr, "    Min call time: " PRINTF_FORMAT_LLU "ms\n", m_min_time);
    fprintf(stderr, "    Max call time: " PRINTF_FORMAT_LLU "ms\n", m_max_time);
    fprintf(stderr, "    Average call time: " PRINTF_FORMAT_LLU "ms\n", m_total_time / m_calls);

    if (m_total_time > 3000)
        fprintf(stderr, "    Total time: \033[1;31m" PRINTF_FORMAT_LLU "ms\033[0m\n", m_total_time);
    else if (m_total_time > 2000)
        fprintf(stderr, "    Total time: \033[0;31m" PRINTF_FORMAT_LLU "ms\033[0m\n", m_total_time);
    else if (m_total_time > 1000)
        fprintf(stderr, "    Total time: \033[0;33m" PRINTF_FORMAT_LLU "ms\033[0m\n", m_total_time);
    else
        fprintf(stderr, "    Total time: " PRINTF_FORMAT_LLU "ms\n", m_total_time);
}

SoftPoint::SoftPoint(HardPoint *hard_point)
    : m_hard_point(hard_point)
{
    ++m_hard_point->m_calls;
    m_start = get_tick_count();
}

SoftPoint::~SoftPoint()
{
    unsigned long long t = get_tick_count() - m_start;

    if (t < m_hard_point->m_min_time)
        m_hard_point->m_min_time = t;

    if (t > m_hard_point->m_max_time)
        m_hard_point->m_max_time = t;

    m_hard_point->m_total_time += t;
}
} // namespace CS
