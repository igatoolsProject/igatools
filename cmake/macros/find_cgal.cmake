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
# Find CGAL library (Optional)
#-------------------------------------------------------------------------------
macro(find_cgal)
  set(CGAL_PREFIX $ENV{CGAL_DIR} CACHE LOCATION 
    "Location where CGAL library is installed")
  find_package(CGAL REQUIRED PATHS ${CGAL_PREFIX})
  message(STATUS "Found CGAL:  version ${CGAL_VERSION}.")

  # add the location of CGAL headers to the include directories
  include_directories(${CGAL_INCLUDE_DIRS})

endmacro(find_cgal)
