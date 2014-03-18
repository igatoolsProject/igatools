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
# Main cmake file for the igatools library
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#TODO: Define the configure time variable for dimension compilation
#TODO: split this CMakeList into modules to make it more understandable

cmake_minimum_required(VERSION 2.8.10)

project(iga CXX)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Include all local macros
file(GLOB_RECURSE macro_list ${PROJECT_SOURCE_DIR}/cmake/*.cmake)
foreach(item ${macro_list})
  include (${item})
endforeach(item)
#
#-------------------------------------------------------------------------------

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Make sure the build directory is not the source or
# the install directory
#
macro_ensure_out_of_source_build()
macro_ensure_out_of_build_install()
#
#-------------------------------------------------------------------------------


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# If the prefix is not provided by the user but the
# enviroment variable IGA_PREFIX exists the we use it
# instead of the default prefix
if((DEFINED ENV{IGA_PREFIX}) AND (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT))
set(CMAKE_INSTALL_PREFIX 
   $ENV{IGA_PREFIX} CACHE PATH "Installation prefix" FORCE)
endif()
#
#-------------------------------------------------------------------------------


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Set the igatools version 
#
set(IGA_MAJOR_VERSION 0)
set(IGA_MINOR_VERSION 3)
set(IGA_PATCH_VERSION 0)
set(IGA_VERSION
  ${IGA_MAJOR_VERSION}.${IGA_MINOR_VERSION}.${IGA_PATCH_VERSION})
#
#-------------------------------------------------------------------------------


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Compiler and compile flags
init_cxx_flags()
#
#-------------------------------------------------------------------------------


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Find Required programs
#
find_sed()
find_python()
find_doxygen()
#
#-------------------------------------------------------------------------------


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Find Required libraries
#
find_trilinos()
find_boost()
#
#-------------------------------------------------------------------------------


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Optional features for building  igatools
#
# Choose the precision for the real values
option(REAL_IS_LONG_DOUBLE "Enable quadruple precision" OFF)

option(submit_test "Allow test submition" OFF)

#
#-------------------------------------------------------------------------------


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#  spaces dimensions and ranks to compile the library for
#
set(dim   CACHE STRINGS "list of domain dimensions")
set(codim CACHE STRINGS "list of domain codimension")
set(range CACHE STRINGS "list of reference space range dim")
set(rank  CACHE STRINGS "list of reference space rank")
set(trans_type  CACHE STRINGS "list of transformation type")
set(max_der_order 2 CACHE STRINGS "maximum order of derivatives")
#
#-------------------------------------------------------------------------------


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Configuration of the library 
#
set(lib_name igatools)

# igatools source files
set (source_dirs base utils geometry basis_functions linear_algebra io contrib)
foreach(dir ${source_dirs})
  file(GLOB source ${PROJECT_SOURCE_DIR}/source/${dir}/*.cpp)
  list(APPEND source_files ${source})
endforeach()

# igatools include directories
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include 
                    ${CMAKE_CURRENT_BINARY_DIR}/include)

# the igatools library targent
add_library(${lib_name} SHARED  ${source_files})
target_link_libraries(${lib_name} 
  ${Trilinos_LIBRARIES} 
  ${Trilinos_TPL_LIBRARIES} 
  ${Boost_LIBRARIES})
set_property(TARGET ${lib_name} PROPERTY VERSION ${IGA_VERSION})


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Generate configure dependent source files (in the build directory)
configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/include/${lib_name}/base/config.h.in"
  "${CMAKE_CURRENT_BINARY_DIR}/include/${lib_name}/base/config.h")
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/config/${lib_name}Config.cmake.in 
  ${CMAKE_CURRENT_BINARY_DIR}/${lib_name}Config.cmake @ONLY)
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/config/${lib_name}ConfigVersion.cmake.in 
  ${CMAKE_CURRENT_BINARY_DIR}/${lib_name}ConfigVersion.cmake @ONLY)





# Create the instantiation files
message("Configuring instantiations")
create_space_table()
generate_instantiations()
#
#-------------------------------------------------------------------------------


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Create install target
create_install_target()
#
#-------------------------------------------------------------------------------


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Create target to generate documentation
create_doc_target()
#
#-------------------------------------------------------------------------------


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Tutorial example section
add_subdirectory(${PROJECT_SOURCE_DIR}/tutorial)
#
#-------------------------------------------------------------------------------


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Unit Test configuration
add_subdirectory(${PROJECT_SOURCE_DIR}/tests)
#
#-------------------------------------------------------------------------------


print_final_message()