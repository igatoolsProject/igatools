#-+--------------------------------------------------------------------
# Igatools a general purpose Isogeometric analysis library.
# Copyright (C) 2012-2015  by the igatools authors (see authors.txt).
#
# This file is part of the igatools library.
#
# The igatools library is free software: you can use it, redistribute
# it and/or modify it under the terms of the GNU General Public
# License as published by the Free Software Foundation, either
# version 3 of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#-+--------------------------------------------------------------------

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Find Boost library (Required)
#-------------------------------------------------------------------------------
macro(find_boost)

    message("Serialization ${SERIALIZATION}")

  if (SERIALIZATION)
    find_package(Boost 1.54.0 REQUIRED COMPONENTS serialization)
    include_directories(${Boost_INCLUDE_DIRS})
    find_library (Boost_serialization boost_serialization ${Boost_LIBRARY_DIRS})
    message("${Boost_serialization}")
    set(Boost_LIBRARIES "${Boost_serialization}")
  else ()
    message("${Boost_LIBRARY_DIRS}")
    find_package(Boost 1.54.0 REQUIRED)
    include_directories(${Boost_INCLUDE_DIRS})
  endif ()
  
endmacro(find_boost)
