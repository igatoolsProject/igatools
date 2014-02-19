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

project(tutorial)

file(GLOB dirs "${CMAKE_CURRENT_SOURCE_DIR}/example*")

foreach(dir ${dirs})
  set(IGA_TUTORIAL_DIRS "${IGA_TUTORIAL_DIRS} ${dir}")
endforeach()
set(IGA_TUTORIAL_DIRS "${IGA_TUTORIAL_DIRS}" PARENT_SCOPE)

message("Configuring Tutorial")

foreach(dir ${dirs})
  get_filename_component(name ${dir} NAME_WE)
  set(example_name ${name})
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cmake_template.txt
    ${CMAKE_CURRENT_BINARY_DIR}/${name}/CMakeLists.txt @ONLY)
  
  execute_process(COMMAND ${SED_EXECUTABLE} -e  \\,//,d
    INPUT_FILE ${dir}/${name}.cpp
    OUTPUT_FILE ${dir}/${name}_no_tag.cpp)
    
  configure_file(${dir}/${name}_no_tag.cpp 
    ${CMAKE_CURRENT_BINARY_DIR}/${name}/${name}.cpp)

  # Add data files if any
  file(GLOB files "${dir}/*.xml")
  file(COPY ${files} DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/${name})
  add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/${name}
              ${CMAKE_CURRENT_BINARY_DIR}/${name})
endforeach()
