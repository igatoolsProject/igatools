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
# Find PETSc library (Optional)
#-------------------------------------------------------------------------------
macro(find_petsc)
  set(Petsc_INCLUDE_DIRS $ENV{PETSC_INCLUDE_DIRS} CACHE LOCATION 
    "Location where PETSc headers are installed")

  set(Petsc_LIBRARY_DIRS $ENV{PETSC_LIBRARY_DIRS} CACHE LOCATION 
    "Location where PETSc libraries are installed")
   
  include_directories(${Petsc_INCLUDE_DIRS} ${Petsc_INCLUDE_DIRS}/$ENV{PETSC_ARCH}/include $ENV{PETSC_DIR}/$ENV{PETSC_ARCH}/include)
  link_directories(${Petsc_LIBRARY_DIRS})

  set(Petsc_LIBRARIES libpetsc.so)
    
  message("-- Found PETSc library in ${Petsc_LIBRARY_DIRS}")
    
endmacro(find_petsc)
