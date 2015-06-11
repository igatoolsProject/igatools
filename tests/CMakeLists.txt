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
# cmakefile for the igatools library unit testing configuration
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

project(tests)

message("Configuring unit tests")

#--------------------------------------------------------------
# Find the igatools library
set(CMAKE_PREFIX_PATH ${iga_BINARY_DIR} ${CMAKE_PREFIX_PATH})

find_package(igatools REQUIRED)

include_directories(${IGATOOLS_INCLUDE_DIRS})
link_directories(${IGATOOLS_LIBRARY_DIR})


#--------------------------------------------------------------
# Add a target for each .cpp file in test_dirs/*.cpp
#
set (test_dirs base basis_functions utils
  equations geometry linear_algebra io 
  tutorial templates assemble functions refinement)

foreach(dir ${test_dirs})
  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${dir})
  file(GLOB files "${CMAKE_CURRENT_SOURCE_DIR}/${dir}/*.cpp")
  foreach(filename ${files})
    get_filename_component(name ${filename} NAME_WE)
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${dir}/${name})
    add_executable(${dir}/${name}.exe EXCLUDE_FROM_ALL ${dir}/${name}.cpp)
    target_link_libraries (${dir}/${name}.exe  ${IGATOOLS_LIBRARIES} )
  endforeach()
endforeach()
#--------------------------------------------------------------



# Create the tests 
if (submit_test)
  include(CTest)
endif()
enable_testing()
foreach(dir ${test_dirs})
  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${dir})
  file(GLOB files "${CMAKE_CURRENT_SOURCE_DIR}/${dir}/*.cpp")
  foreach(filename ${files})
    get_filename_component(name ${filename} NAME_WE)
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/${dir}/${name}/)
    file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/${dir}/${name}
      DESTINATION ${CMAKE_CURRENT_BINARY_DIR}/${dir}/
      FILES_MATCHING PATTERN "*.xml")
    add_test(
      NAME ${dir}-${name}
      COMMAND ${CMAKE_COMMAND}
      -DEXE_TARGET=${dir}/${name}.exe
      -DTEST_PROG=${CMAKE_CURRENT_BINARY_DIR}/${dir}/${name}.exe
      -DSOURCE_DIR=${CMAKE_CURRENT_SOURCE_DIR}/${dir}/${name}
      -DWORK_DIR=${CMAKE_CURRENT_BINARY_DIR}/${dir}/${name}
      -P ${CMAKE_CURRENT_SOURCE_DIR}/runtest.cmake)
  endforeach()
endforeach()


 
 
