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
# cmakefile for the igatools-ParaView plugin
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


#-------------------------------------------------------------------------------
# Set the igatools-ParaView version 
#
set(IGATOOLS_PARAVIEW_MAJOR_VERSION 0)
set(IGATOOLS_PARAVIEW_MINOR_VERSION 1)
set(IGATOOLS_PARAVIEW_PATCH_VERSION 0)
set(IGATOOLS_PARAVIEW_VERSION
  ${IGATOOLS_PARAVIEW_MAJOR_VERSION}.${IGATOOLS_PARAVIEW_MINOR_VERSION}.${IGATOOLS_PARAVIEW_PATCH_VERSION})
#
#-------------------------------------------------------------------------------

set(igatools_paraview_name vtkIgatoolsReader.${IGATOOLS_PARAVIEW_VERSION})
set(igatools_paraview_lib_name vtkIgatoolsReader.${IGATOOLS_PARAVIEW_VERSION})

project (${igatools_paraview_name} CXX)

message("Configuring igatools-ParaView plugin")

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Include all local macros
file(GLOB_RECURSE macro_list ${PROJECT_SOURCE_DIR}/cmake/*.cmake)
foreach(item ${macro_list})
  include (${item})
endforeach(item)
#
#-------------------------------------------------------------------------------

#--------------------------------------------------------------
# Finding the needed libraries.
find_igatools ()

# find_trilinos ()

find_paraview ()
#-------------------------------------------------------------------------------


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Compiler and compile flags
init_cxx_flags()
#
#-------------------------------------------------------------------------------


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Create the configuration, install and documentation targets

create_configuration_target ()

create_install_target()

create_doc_target()
#-------------------------------------------------------------------------------