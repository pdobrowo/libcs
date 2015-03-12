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
#include <cs/Logger.h>
#include <cstdio>

#if defined(__linux__) || defined(__FreeBSD__) || defined(__OpenBSD__)

#include <sys/time.h>

static unsigned long long sys_time_usec()
{
    struct timeval tv;
    (void)gettimeofday(&tv, 0);
    return static_cast<unsigned long long>(tv.tv_sec) * 1000000 + static_cast<unsigned long long>(tv.tv_usec);
}

#endif // defined(__linux__) || defined(__FreeBSD__) || defined(__OpenBSD__)

#ifdef _WIN32

#ifndef GetTickCount
extern "C" unsigned long __stdcall GetTickCount();
#endif // GetTickCount

static unsigned long long sys_time_usec()
{
    return static_cast<unsigned long long>(GetTickCount()) * 1000;
}

#endif // _WIN32

static unsigned long long g_logger_base_timestamp = sys_time_usec();

static void logger_print(const std::string &type, const std::string &module, const std::string &msg)
{
    std::printf("[%8llu us] [%7s] %s: %s\n", sys_time_usec() - g_logger_base_timestamp, type.c_str(), module.c_str(), msg.c_str());
}

namespace CS
{
namespace Logger
{
void debug(const std::string &module, const std::string &msg)
{
    logger_print("DEBUG", module, msg);
}

void info(const std::string &module, const std::string &msg)
{
    logger_print("INFO", module, msg);
}

void warning(const std::string &module, const std::string &msg)
{
    logger_print("WARNING", module, msg);
}
} // namespace Logger
} // namespace CS
