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
#include <log4cxx/logger.h>
#include <log4cxx/propertyconfigurator.h>

//#include "test_various.h"
//#include "test_exact_perf.h"
//#include "test_cell_perf.h"
#include "test_inexact_perf.h"

int main()
{
    // setup global logger
    log4cxx::LoggerPtr g_logger(log4cxx::Logger::getLogger("Test.main"));

    // configure logger
    log4cxx::PropertyConfigurator::configure("logger.conf");

    //test_exact_perf();
    //test_cell_perf();
    test_inexact_perf();

    return 0;
}
