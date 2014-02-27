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
#include <cstdlib>
#include <iostream>
#include <log4cxx/logger.h>
#include <log4cxx/basicconfigurator.h>
#include <log4cxx/propertyconfigurator.h>
#include <log4cxx/helpers/exception.h>
#include <cs/Build.h>
#include "test_various.h"
#include "test_exact_perf.h"
#include "test_cell_perf.h"

static log4cxx::LoggerPtr logger(log4cxx::Logger::getLogger("Test.main"));

const char *build_string()
{
#if defined(_DEBUG) || !defined(NDEBUG)
    return "test " __TIMESTAMP__ " Debug";
#else
    return "test " __TIMESTAMP__ " Release";
#endif
}

int main(int argc, char *argv[])
{
    (void)argc;
    (void)argv;

    // configure logger
    log4cxx::PropertyConfigurator::configure("logger.conf");

    LOG4CXX_INFO(logger, "$lib cs  : " << CS::build_string());
    LOG4CXX_INFO(logger, "$bin test: " << build_string());

/*
    int param;

    if (argc != 2)
    {
        std::cout << "test [param]: invalid parameters" << std::endl;
        return 0;
    }

    param = atoi(argv[1]);

    test_various(param);
*/

    //test_exact_perf();
    test_cell_perf();

    return 0;
}
