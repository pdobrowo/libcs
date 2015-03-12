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
#ifndef LIBCS_LOGGER_H
#define LIBCS_LOGGER_H

#include <sstream>
#include <string>

namespace CS
{
namespace Logger
{
void debug(const std::string &module, const std::string &msg);
void info(const std::string &module, const std::string &msg);
void warning(const std::string &module, const std::string &msg);
} // namespace Logger
} // namespace CS

#define CS_logger_debug(Module, Msg) do { std::ostringstream logger_stream; logger_stream << Msg; CS::Logger::debug(Module, logger_stream.str()); } while (false)
#define CS_logger_info(Module, Msg) do { std::ostringstream logger_stream; logger_stream << Msg; CS::Logger::info(Module, logger_stream.str()); } while (false)
#define CS_logger_warning(Module, Msg) do { std::ostringstream logger_stream; logger_stream << Msg; CS::Logger::warning(Module, logger_stream.str()); } while (false)

#endif // LIBCS_LOGGER_H
