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
#
# The previous copyright of this file is shown below.
# It has been modified by the igatools authors to fit the igatools framework.
#-+--------------------------------------------------------------------

#-+--------------------------------------------------------------------
# Copyright (c) 2009, Ben Morgan, <Ben.Morgan@warwick.ac.uk>
#
# Redistribution and use is allowed according to the terms of the BSD license.
# For details see the accompanying COPYING-CMAKE-SCRIPTS file.
#-+--------------------------------------------------------------------

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Find xercesc library (Required)
#-------------------------------------------------------------------------------
macro (find_xercesc)

# - Find Xerces-C
# This module tries to find the Xerces-C library and headers.
# Once done this will define
#
#   XERCESC_FOUND - system has Xerces-C headers and libraries
#   XERCESC_INCLUDE_DIRS - the include directories needed for Xerces-C
#   XERCESC_LIBRARIES - the libraries needed to use Xerces-C
#
# Variables used by this module, which can change the default behaviour and
# need to be set before calling find_package:
#
#   XERCESC_ROOT_DIR            Root directory to Xerces-C installation. Will
#                               be used ahead of CMake default path.
#
# The following advanced variables may be used if the module has difficulty
# locating Xerces-C or you need fine control over what is used.
#
#   XERCESC_INCLUDE_DIR
#
#   XERCESC_LIBRARY
#

# Look for the header - preferentially searching below XERCESC_ROOT_DIR
find_path(
    XERCESC_INCLUDE_DIR 
    NAMES xercesc/util/XercesVersion.hpp
    PATHS ${XERCESC_ROOT_DIR}
    PATH_SUFFIXES include
    NO_DEFAULT_PATH
)

# If we didn't find it there, fall back to some standard search paths
find_path(
    XERCESC_INCLUDE_DIR
    NAMES xercesc/util/XercesVersion.hpp
)

# Look for the library, preferentially searching below XERCESC_ROOT_DIR
find_library(
    XERCESC_LIBRARY
    NAMES xerces-c xerces-c_3
    PATHS ${XERCESC_ROOT_DIR}
    PATH_SUFFIXES lib64 lib32 lib
    NO_DEFAULT_PATH
)

find_library(
    XERCESC_LIBRARY
    NAMES xerces-c xerces-c_3
)


include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
    XercesC
    DEFAULT_MSG
    XERCESC_LIBRARY
    XERCESC_INCLUDE_DIR
)

if (XERCESC_FOUND)
    set(XERCESC_LIBRARIES ${XERCESC_LIBRARY})
    set(XERCESC_INCLUDE_DIRS ${XERCESC_INCLUDE_DIR})
else (XERCESC_FOUND)
    set(XERCESC_LIBRARIES)
    set(XERCESC_INCLUDE_DIRS)
endif (XERCESC_FOUND)


mark_as_advanced(
    XERCESC_LIBRARY 
    XERCESC_INCLUDE_DIR
)

endmacro (find_xercesc)