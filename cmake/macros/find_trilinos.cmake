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
# Find Trilinos library (Optional)
#-------------------------------------------------------------------------------
macro(find_trilinos)
  set(Trilinos_PREFIX $ENV{TRILINOS_PREFIX} CACHE LOCATION 
    "Location where Trilinos library is installed")
  find_package(Trilinos 11 REQUIRED PATHS ${Trilinos_PREFIX})
  message(STATUS "Found Trilinos:  version ${Trilinos_VERSION}.")
  if (NOT (Trilinos_VERSION VERSION_GREATER 11.6))
    message(FATAL_ERROR "Trilinos 11.6.1 or greater is required.")
  endif()

  # add the location of Trilinos headers to the include directories
  include_directories( ${Trilinos_INCLUDE_DIRS} ${Trilinos_TPL_INCLUDE_DIRS})

  # Check that individual required Trilinos packages are available
  set(tri_required_packages Tpetra Epetra Belos Amesos Amesos2)
  foreach(package ${tri_required_packages})
    list(FIND Trilinos_PACKAGE_LIST  ${package} package_index)
    if (package_index EQUAL -1)
      message(FATAL_ERROR "Trilinos ${package} package not found.")
    endif()
  endforeach(package)
endmacro(find_trilinos)
