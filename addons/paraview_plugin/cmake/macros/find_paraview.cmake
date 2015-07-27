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
# Find Paraview library (Required)
#-------------------------------------------------------------------------------
macro(find_paraview)

  find_package (ParaView REQUIRED QUIET)
  include (${PARAVIEW_USE_FILE})

  if (ParaView_SOURCE_DIR)
    include_directories (
      ${PARAVIEW_INCLUDE_DIRS}
      ${PARAVIEW_GUI_INCLUDE_DIRS}
      ${PARAVIEW_KWSYS_INCLUDE_DIRS}
      ${VTK_INCLUDE_DIRS})
  else ()
    find_package (ParaView REQUIRED)
    include (${PARAVIEW_USE_FILE})
  endif ()

  if (PARAVIEW_BUILD_QT_GUI)
    if (PARAVIEW_QT_VERSION VERSION_GREATER "4")
      set (Qt5_FIND_COMPONENTS Widgets)
      include (ParaViewQt5)
    else ()
      include (${QT_USE_FILE})
    endif ()
  endif ()

  message(STATUS "Found ParaView:  version ${ParaView_VERSION}.")
endmacro(find_paraview)
