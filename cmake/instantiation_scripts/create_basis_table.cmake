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
# Creates a file with a table of the spaces that the library
# is suppossed to be compiled for.
#
macro(create_basis_table)
  message(STATUS "Generating physical basis table.")
  execute_process(COMMAND 
    ${PYTHON_EXECUTABLE} 
    -B ${PROJECT_SOURCE_DIR}/cmake/instantiation_scripts/generate_inst_table.py 
    dim_ref=${dim}
    range_ref=${range}
    rank_ref=${rank}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    RESULT_VARIABLE res) 
endmacro(create_basis_table)
