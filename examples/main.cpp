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
#include "example_simple_raster.h"
#include "example_gearbox.h"
#include "example_ball_knob.h"
#include "example_find_path.h"
#include <log4cxx/propertyconfigurator.h>

int main()
{
    log4cxx::PropertyConfigurator::configure("logger.conf");

    //example_simple_raster();
    //example_gearbox();
    //example_ball_knob();
    example_find_path();

    return 0;
}
