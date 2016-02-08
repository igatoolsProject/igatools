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

#------------------------------------------------------------------------------
# Find and Use ParaView
#------------------------------------------------------------------------------

cmake_minimum_required(VERSION 2.8.10)

set(PLUGIN_SOURCE_DIR ${PROJECT_SOURCE_DIR})

# Redefining paths: pointing to the igatools source directory.
set(PROJECT_SOURCE_DIR ${PLUGIN_SOURCE_DIR}/../..)

#-------------------------------------------------------------------------------
# Set the igatools-ParaView version 
#
set(IGATOOLS_PARAVIEW_MAJOR_VERSION 1)
set(IGATOOLS_PARAVIEW_MINOR_VERSION 0)
set(IGATOOLS_PARAVIEW_PATCH_VERSION 0)
set(IGATOOLS_PARAVIEW_VERSION
  ${IGATOOLS_PARAVIEW_MAJOR_VERSION}.${IGATOOLS_PARAVIEW_MINOR_VERSION}.${IGATOOLS_PARAVIEW_PATCH_VERSION})
#
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Set the igatools required version 
#
set(IGATOOLS_REQUIRED_MAJOR_VERSION 1)
set(IGATOOLS_REQUIRED_MINOR_VERSION 0)
set(IGATOOLS_REQUIRED_PATCH_VERSION 0)
set(IGATOOLS_REQUIRED_VERSION
  ${IGATOOLS_REQUIRED_MAJOR_VERSION}.${IGATOOLS_REQUIRED_MINOR_VERSION}.${IGATOOLS_REQUIRED_PATCH_VERSION})
#
#-------------------------------------------------------------------------------

set(igatools_paraview_name IgatoolsParaViewReader)
set(igatools_paraview_lib_name igatools_paraview)
set(igatools_paraview_include_name paraview_plugin)
# project (${igatools_paraview_name} CXX)

message("Configuring igatools-ParaView plugin")

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Include all local macros from the paraview plugin (and also from igatools)

file(GLOB_RECURSE macro_list ${PLUGIN_SOURCE_DIR}/cmake/*.cmake
                             ${PROJECT_SOURCE_DIR}/cmake/*.cmake)
foreach(item ${macro_list})
  include (${item})
endforeach(item)

#-------------------------------------------------------------------------------


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Make sure the build directory is not the source or
# the install directory

macro_ensure_out_of_source_build()
macro_ensure_out_of_build_install()

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
# Required libraries and programs

find_igatools (${IGATOOLS_REQUIRED_VERSION})

find_paraview ()

find_sed()

find_python()
#-------------------------------------------------------------------------------


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Compiler and compile flags

if(CMAKE_BUILD_TYPE MATCHES Release)
  message("-- Igatools ParaView plugin build type: RELEASE")
else ()
  message("-- Igatools ParaView plugin type: DEBUG")
endif()

init_cxx_flags()

#-------------------------------------------------------------------------------

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Create the configuration, install and documentation targets and finishing

create_paraview_plugin_configuration_target ()

set(max_der_order 2 CACHE STRINGS "maximum order of derivatives")
message("Configuring instantiations")
create_basis_table()
generate_plugin_instantiations()

create_paraview_plugin_install_target()

print_paraview_plugin_final_message()
#-------------------------------------------------------------------------------
