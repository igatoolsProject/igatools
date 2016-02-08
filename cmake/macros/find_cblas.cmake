#-+--------------------------------------------------------------------
# Igatools a general purpose Isogeometric analysis library.
# Copyright (C) 2012-2016  by the igatools authors (see authors.txt).
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
# Find CBLAS library (Required)
#-------------------------------------------------------------------------------
macro(find_cblas)
  set(CBLAS_SEARCH_PATHS /usr/lib)
  find_library(CBLAS_LIBRARY NAMES cblas PATHS ${CBLAS_SEARCH_PATHS})
  
  if (CBLAS_LIBRARY)
    set(CBLAS_FOUND ON)
    set(CBLAS_LIBRARY_DIR ${CBLAS_SEARCH_PATHS})
  endif (CBLAS_LIBRARY)

  if (CBLAS_FOUND)
    message(STATUS "Found CBLAS in ${CBLAS_LIBRARY}")
  else (CBLAS_FOUND)
    message(FATAL_ERROR "Could not find CBLAS libraries.")
  endif (CBLAS_FOUND)

#  mark_as_advanced(CBLAS_LIBRARY)
endmacro(find_cblas)
