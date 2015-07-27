
#------------------------------------------------------------------------------
# Find and Use ParaView
#------------------------------------------------------------------------------

cmake_minimum_required(VERSION 2.8.10)

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
# project (${igatools_paraview_name} CXX)

message("Configuring igatools-ParaView plugin")

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Include all local macros from the paraview plugin (and also from igatools)

file(GLOB_RECURSE macro_list ${PROJECT_SOURCE_DIR}/cmake/*.cmake
                             ${PROJECT_SOURCE_DIR}/../../cmake/*.cmake)
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

create_paraview_plugin_install_target()

print_paraview_plugin_final_message()
#-------------------------------------------------------------------------------
