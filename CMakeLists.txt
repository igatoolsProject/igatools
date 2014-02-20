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
# Main cmake file for the igatool library
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#TODO: Define the configure time variable for dimension compilation
#TODO: split this CMakeList into modules to make it more understandable

cmake_minimum_required(VERSION 2.8)

project(iga CXX)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Include all local macros
file(GLOB macro_list ${PROJECT_SOURCE_DIR}/cmake/macros/*.cmake)
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
set(IGA_MINOR_VERSION 2)
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
option(build_tutorial "Build tutorial" ON)
option(build_test "Build test" ON)
option(submit_test "Allow test submition" OFF)

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Igatools specific configuration 
#



# Choose the precision for the real values
option(REAL_IS_LONG_DOUBLE "Enable quadruple precision" OFF)

set(lib_name igatools)

configure_file (
  "${CMAKE_CURRENT_SOURCE_DIR}/include/${lib_name}/base/config.h.in"
  "${CMAKE_CURRENT_BINARY_DIR}/include/${lib_name}/base/config.h")

# Collect all igatools source files, contained in the directory list below

set (source_dirs base)
foreach(dir ${source_dirs})
  file(GLOB source ${PROJECT_SOURCE_DIR}/source/${dir}/*.cpp)
  list(APPEND source_files ${source})
endforeach()


# Pass the igatools include directories
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include 
                    ${CMAKE_CURRENT_BINARY_DIR}/include
                    ${addons_include_dir})

# Add the igatools library
add_library(${lib_name} SHARED  ${source_files})
target_link_libraries( ${lib_name} 
  ${Trilinos_LIBRARIES} ${Trilinos_TPL_LIBRARIES} ${Boost_LIBRARIES}
  ${TBB_LIBRARIES})
set_property( TARGET ${libname} PROPERTY VERSION ${IGA_VERSION} )

# todo: put in a CMakelist.txt in the config directory
# generate the cmake scripts for external applications to find igatools
# install will put them in the install/lib directory
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/config/${lib_name}Config.cmake.in 
  ${CMAKE_CURRENT_BINARY_DIR}/${lib_name}Config.cmake @ONLY)
configure_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/config/${lib_name}ConfigVersion.cmake.in 
  ${CMAKE_CURRENT_BINARY_DIR}/${lib_name}ConfigVersion.cmake @ONLY)
#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Dealing with template instatiations

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Specifying the spaces for which we wish to compile the library
#
set(dim_ref   CACHE STRINGS "list of reference space dom dim")
set(range_ref CACHE STRINGS "list of reference space range dim")
set(rank_ref  CACHE STRINGS "list of reference space rank")
set(dim_phys  CACHE STRINGS "list of physical embedding space dim")
set(trans_type  CACHE STRINGS "list of transformation type")
set(max_der_order 2 CACHE STRINGS "maximum order of derivatives")

#TODO: these properties should be marked as advenced
#set_property(CACHE dims1 PROPERTY ADVANCED)


message("Configuring instantiations")
message(STATUS "Generating physical space table to create instantiations.")
execute_process(COMMAND 
  ${PYTHON_EXECUTABLE} 
  -B ${PROJECT_SOURCE_DIR}/cmake/instantiation_scripts/generate_inst_table.py 
  dim_ref=${dim_ref}
  range_ref=${range_ref}
  rank_ref=${rank_ref}
  WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
  RESULT_VARIABLE res) 

include (${PROJECT_SOURCE_DIR}/cmake/instantiation_scripts/generate_instantiations.cmake)
generate_instantiations()



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# set the directory where the lib and include directories will be installed 
install( TARGETS   ${lib_name}
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
#
#-------------------------------------------------------------------------------

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Tutorial example section

if(build_tutorial)
  add_subdirectory(${PROJECT_SOURCE_DIR}/tutorial)
endif(build_tutorial)
#
#-------------------------------------------------------------------------------


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Unit Test configuration

add_subdirectory(${PROJECT_SOURCE_DIR}/tests)
#
#-------------------------------------------------------------------------------



#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Create a target to generate the online documentation 
# using Doxygen
find_package(Doxygen 1.8.1 REQUIRED)
if(DOXYGEN_FOUND)
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
    COMMENT "The documentation is found at: ${CMAKE_INSTALL_PREFIX}/doc" VERBATIM)
  file(WRITE ${CMAKE_BINARY_DIR}/install_doc.cmake
    "file( INSTALL ${CMAKE_BINARY_DIR}/doc/html
    DESTINATION ${CMAKE_INSTALL_PREFIX}/doc)\n"
    "file( REMOVE_RECURSE ${CMAKE_BINARY_DIR}/doc/html)")
endif(DOXYGEN_FOUND)
#
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


print_final_message()