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

# Ensures that we do an out of build install

macro( MACRO_ENSURE_OUT_OF_BUILD_INSTALL )
    string( COMPARE EQUAL "${CMAKE_BINARY_DIR}" "${CMAKE_INSTALL_PREFIX}" inbuild )
    if( inbuild )
        MESSAGE(FATAL_ERROR "The install has to be different from the build directory" )
    endif( inbuild )
endmacro(MACRO_ENSURE_OUT_OF_BUILD_INSTALL)

