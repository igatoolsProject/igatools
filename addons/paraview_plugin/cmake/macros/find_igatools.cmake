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
# Find Igatools library (Required)
#-------------------------------------------------------------------------------
macro(find_igatools)

  if (IGATOOLS_PREFIX)
    set(igatools_PREFIX ${IGATOOLS_PREFIX} CACHE LOCATION 
      "Location where igatools library is installed")
  else ()
    set(igatools_PREFIX $ENV{IGATOOLS_PREFIX} CACHE LOCATION 
      "Location where igatools library is installed")
  endif ()

  set(CMAKE_PREFIX_PATH ${igatools_PREFIX}/lib ${CMAKE_PREFIX_PATH} )

  find_package(igatools 1.0.0 REQUIRED PATHS ${igatools_PREFIX})
  link_directories(${IGATOOLS_LIBRARY_DIR})
  include_directories(${IGATOOLS_INCLUDE_DIRS})

endmacro(find_igatools)

