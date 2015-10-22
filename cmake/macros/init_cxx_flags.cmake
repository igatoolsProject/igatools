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
# Compiler and compile flags
macro(init_cxx_flags)
  # Initialize CXXFLAGS.
  set(CMAKE_CXX_FLAGS                "-Wall")
  set(CMAKE_CXX_FLAGS_DEBUG          "-Wall -O0 -g -ftemplate-backtrace-limit=0")
  set(CMAKE_CXX_FLAGS_MINSIZEREL     "-Os -DNDEBUG")
  set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3 -mtune=native -DNDEBUG -g -p -pg")

  # Compiler-specific C++14 activation.
  if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
    execute_process(
      COMMAND ${CMAKE_CXX_COMPILER} -dumpversion OUTPUT_VARIABLE GCC_VERSION)
    set(CMAKE_CXX_FLAGS_RELEASE        "-O3 -mtune=native -DNDEBUG -floop-parallelize-all -ftree-vectorize -fpeel-loops")
    if (NOT (GCC_VERSION VERSION_GREATER 5.1.0 OR GCC_VERSION VERSION_EQUAL 5.1.0))
      message(FATAL_ERROR "To compile igatools is required g++ version 5.1.0 or greater.")
    endif()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -fuse-ld=gold")
  elseif ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14 -fuse-ld=gold -stdlib=libstdc++ -ftemplate-depth=1024")
    set(CMAKE_CXX_FLAGS_RELEASE        "-O3 -mtune=native -DNDEBUG -ftree-vectorize")
  elseif("${CMAKE_CXX_COMPILER_ID}" MATCHES "Intel")
#    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++14")
    message(FATAL_ERROR "Your C++ INTEL compiler is not supported by igatools")
  else()
    message(FATAL_ERROR "Your C++ compiler is not supported by igatools")
  endif()
endmacro(init_cxx_flags)
