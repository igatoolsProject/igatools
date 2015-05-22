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
# Main cmake file for the igatools library
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#TODO: Define the configure time variable for dimension compilation

cmake_minimum_required(VERSION 2.8.10)

project(iga CXX C)

list (APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/Modules")
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
# enviroment variable IGATOOLS_PREFIX exists then we use it
# instead of the default prefix
if((DEFINED ENV{IGATOOLS_PREFIX}) AND (CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT))
set(CMAKE_INSTALL_PREFIX 
   $ENV{IGATOOLS_PREFIX} CACHE PATH "Installation prefix" FORCE)
endif()
#
#-------------------------------------------------------------------------------


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Set the igatools version 
#
set(IGATOOLS_MAJOR_VERSION 1)
set(IGATOOLS_MINOR_VERSION 0)
set(IGATOOLS_PATCH_VERSION 0)
set(IGATOOLS_VERSION
  ${IGATOOLS_MAJOR_VERSION}.${IGATOOLS_MINOR_VERSION}.${IGATOOLS_PATCH_VERSION})
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
# Optional features for building  igatools
#

option(USE_TRILINOS "Use the linear algebra packages from Trilinos" ON)

option(NURBS "Enable NURBS support" ON)
if (NURBS)
    message("-- NURBS support is enabled.")
else(NURBS)
    message("-- NURBS support is not enabled.")    
endif(NURBS)

option(REAL_IS_LONG_DOUBLE "Enable quadruple precision" OFF)
if (REAL_IS_LONG_DOUBLE)
    message("-- Quadruple precision support (EXPERIMENTAL) is enabled.")
else(REAL_IS_LONG_DOUBLE)
    message("-- Quadruple precision support (EXPERIMENTAL) is not enabled.")    
endif(REAL_IS_LONG_DOUBLE)

option(MPATCH "Enable multipatch features of the library (experimental)" OFF)
if (MPATCH)
    message("-- Multipatch features (EXPERIMENTAL) are enabled.")
else(MPATCH)
    message("-- Multipatch features (EXPERIMENTAL) are not enabled.")    
endif(MPATCH)

option(HIERARCHICAL "Enable Hierarchical splines support (experimental)" OFF)
if (HIERARCHICAL)
    message("-- Hierarchical splines support (EXPERIMENTAL) is enabled.")
else(HIERARCHICAL)
    message("-- Hierarchical splines support (EXPERIMENTAL) is not enabled.")    
endif(HIERARCHICAL)

option(SERIALIZATION "Enable serialization support" ON)
if (SERIALIZATION)
    message("-- Serialization support (EXPERIMENTAL) is enabled.")
else(SERIALIZATION)
    message("-- Serialization support (EXPERIMENTAL) is not enabled.")    
endif(SERIALIZATION)

option(USE_PARAVIEW_PLUGIN "Enable the building of the ParaView plugin" OFF)
if (USE_PARAVIEW_PLUGIN)
    if (NOT SERIALIZATION)
        message (FATAL_ERROR "USE_PARAVIEW_PLUGIN variable set to ON is not "
        "compatible with SERIALIZATION set to OFF. Serialization must be "
        "activated.")
    endif ()
    message("-- Paraview plugin building (EXPERIMENTAL) is enabled.")
else(USE_PARAVIEW_PLUGIN)
    message("-- Paraview plugin building (EXPERIMENTAL) is not enabled.")
endif(USE_PARAVIEW_PLUGIN)


option(use_other_templates "Enable a mechanism to use uninstantiated templates" ON)
option(submit_test "Allow test submition" OFF)

#
#-------------------------------------------------------------------------------


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Find Required libraries
#
find_boost(${SERIALIZATION})
#
#-------------------------------------------------------------------------------


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Find optional libraries

if (USE_TRILINOS)
  find_trilinos()
endif()

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
set(num_sub_elem 1 CACHE STRINGS "number of sub element to consider 1 for faces")
#
#-------------------------------------------------------------------------------


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Configuration of the library 
#
set(lib_name igatools)

# igatools source files
set (source_dirs base utils geometry basis_functions linear_algebra io contrib functions)
foreach(dir ${source_dirs})
  file(GLOB source ${PROJECT_SOURCE_DIR}/source/${dir}/*.cpp)
  list(APPEND source_files ${source})
endforeach()

# igatools include directories
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include 
                    ${CMAKE_CURRENT_BINARY_DIR}/include)

# the igatools library target
add_library(${lib_name} SHARED  ${source_files})
target_link_libraries(${lib_name} 
  ${Trilinos_TPL_LIBRARIES} 
  ${Trilinos_LIBRARIES} 
  ${Petsc_LIBRARIES} 
  ${Boost_LIBRARIES})
set_property(TARGET ${lib_name} PROPERTY VERSION ${IGATOOLS_VERSION})


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
create_install_target(${use_other_templates})
#
#-------------------------------------------------------------------------------


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Create target to generate documentation
create_doc_target()
#
#-------------------------------------------------------------------------------

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# ParaView plugin section
if (USE_PARAVIEW_PLUGIN)
  add_subdirectory(${PROJECT_SOURCE_DIR}/addons/paraview_plugin)
endif()
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
