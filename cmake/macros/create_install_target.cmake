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
# Create the install target
macro(create_install_target)
  install(TARGETS   ${lib_name}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
    LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
  install(DIRECTORY   ${PROJECT_SOURCE_DIR}/include/
    DESTINATION ${CMAKE_INSTALL_PREFIX}/include/
    PATTERN ".*" EXCLUDE
    PATTERN "*.in" EXCLUDE)
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/include/${lib_name}/base/config.h 
    DESTINATION ${CMAKE_INSTALL_PREFIX}/include/${lib_name}/base)
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${lib_name}Config.cmake 
    ${CMAKE_CURRENT_BINARY_DIR}/${lib_name}ConfigVersion.cmake
    DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
endmacro(create_install_target)