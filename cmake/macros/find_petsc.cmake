#-+--------------------------------------------------------------------
# Igatools a general purpose Isogeometric analysis library.
# Copyright (C) 2012-2014  by the igatools authors (see authors.txt).
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
# Find PETSc library (Optional)
#-------------------------------------------------------------------------------
macro(find_petsc)
  set(Petsc_PREFIX $ENV{PETSC_PREFIX} CACHE LOCATION 
    "Location where PETSc library is installed")
  find_package(PETSc 11 REQUIRED PATHS ${PETSC_PREFIX})
#  message(STATUS "Found PETSc:  version ${Petsc_VERSION}.")
#  if (NOT (Petsc_VERSION VERSION_GREATER 11.6))
#    message(FATAL_ERROR "Petsc 11.6.1 or greater is required.")
#  endif()

  # add the location of PETSc headers to the include directories
  include_directories( ${Petsc_INCLUDE_DIRS} )

  # Check that individual required PETSc packages are available
#  set(tri_required_packages Tpetra Belos)
#  foreach(package ${tri_required_packages})
#    list(FIND Petsc_PACKAGE_LIST  ${package} package_index)
#    if (package_index EQUAL -1)
#      message(FATAL_ERROR "Petsc ${package} package not found.")
#    endif()
#  endforeach(package)
endmacro(find_petsc)
