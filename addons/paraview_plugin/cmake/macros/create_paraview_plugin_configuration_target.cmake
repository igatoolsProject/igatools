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
# Create the configuration target
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

macro(create_paraview_plugin_configuration_target)

  # igatools Paraview plugin manager source files
  set (manager_source_files_xml "${CMAKE_SOURCE_DIR}/IgatoolsParaViewReader.xml")
  set (manager_source_files_cxx "${CMAKE_SOURCE_DIR}/IgatoolsParaViewReader.cpp")

  # igatools Paraview plugin source files
  set (source_dirs "source")
  foreach(dir ${source_dirs})
    file(GLOB source ${CMAKE_SOURCE_DIR}/${dir}/*.cpp)
    list(APPEND source_files_cxx ${source})
  endforeach()

  include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include
                      ${CMAKE_CURRENT_BINARY_DIR}/include)

  add_paraview_plugin (${igatools_paraview_lib_name} "${IGATOOLS_PARAVIEW_VERSION}"
    SERVER_MANAGER_XML       ${manager_source_files_xml}
    SERVER_MANAGER_SOURCES   ${manager_source_files_cxx}
    SERVER_SOURCES           ${source_files_cxx}
    )

  set_property(TARGET ${igatools_paraview_lib_name} PROPERTY VERSION ${IGATOOLS_PARAVIEW_VERSION})

  target_link_libraries(${igatools_paraview_lib_name}
                        LINK_PRIVATE
                        ${IGATOOLS_LIBRARIES})

endmacro(create_paraview_plugin_configuration_target)
