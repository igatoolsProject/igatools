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
# Create a target to generate the online documentation 
# using Doxygen
#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
macro(create_doc_target)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/doc/Doxyfile.in 
    ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/doc)
  add_custom_command(OUTPUT ${CMAKE_BINARY_DIR}/doc/html/index.html
    COMMAND ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Generating API documentation with Doxygen" VERBATIM)
  add_custom_target(doc COMMAND ${CMAKE_COMMAND} -P install_doc.cmake
    DEPENDS ${CMAKE_BINARY_DIR}/doc/html/index.html 
    ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
    COMMENT "The documentation is found at: ${CMAKE_INSTALL_PREFIX}/doc" 
    VERBATIM)
  file(WRITE ${CMAKE_BINARY_DIR}/install_doc.cmake
    "file( INSTALL ${CMAKE_BINARY_DIR}/doc/html
    DESTINATION ${CMAKE_INSTALL_PREFIX}/doc)\n"
    "file( REMOVE_RECURSE ${CMAKE_BINARY_DIR}/doc/html)")
endmacro(create_doc_target)
