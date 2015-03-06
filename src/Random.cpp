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
#include <cs/Random.h>
#include <cstdlib>
#ifdef __linux__
#include <unistd.h>
#endif // __linux__
#if defined(__FreeBSD__) || defined(__OpenBSD__)
#include <sys/types.h>
#include <unistd.h>
#endif // defined(__FreeBSD__) || defined(__OpenBSD__)
#ifdef _WIN32
#include <windows.h>
int getpid() { return static_cast<int>(GetCurrentProcessId()); }
#endif // _WIN32

int init_srand()
{
    srand(static_cast<unsigned int>(getpid()));
    return 0;
}

static int init_srand_var = init_srand();

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
} // namespace CS
