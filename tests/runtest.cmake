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

# This script is use to perform a test for the library.
# 
# It runs the following steps:
# 1 - compiles and links the test program
# 2 - execute the program in its working directory
# 3 - compares the output of the program with the expected output
#
# 
# It's called from the tests CMakeLists.txt as follows:
# cmake -DTEST_PROG=${CMAKE_CURRENT_BINARY_DIR}/${dir}/${name}.exe
#       -DSOURCEDIR=${CMAKE_SOURCE_DIR}/${dir}/${name}
#       -DWORK_DIR=${CMAKE_CURRENT_BINARY_DIR}/${dir}/${name}
#       -P ${CMAKE_SOURCE_DIR}/runtest.cmake
#
# TODO: make should be CMAKE_MAKE_COMMAND but is empty!!!


# expected_level_of_success 
# 1: not compile
# 2: compile but run time fail
# 3: compile and run
# in all cases the output is checked

set(expected_level 3)


message ("******************************************************") 
set(level 0)
# Structure check, expected directories and files exist
if(NOT EXISTS ${WORK_DIR})
    message(FATAL_ERROR "** Work directory does not exist.")
endif()
if(NOT EXISTS ${SOURCE_DIR}/expected.txt)    
    message(FATAL_ERROR "** expected.txt does not exist.")
endif()

message ("** Compiling and linking test program")
execute_process(
  COMMAND         make ${EXE_TARGET} 
  RESULT_VARIABLE compile_error)

set(runtime_error "Not run")
if(NOT compile_error)
  set(level 2)
  message ("** Executing test program.")
  execute_process(
    COMMAND           ${TEST_PROG} 
    OUTPUT_FILE       stdout.txt
    ERROR_FILE        stdout.txt
    WORKING_DIRECTORY ${WORK_DIR} 
    RESULT_VARIABLE   runtime_error)
endif(NOT compile_error)


if(NOT runtime_error)
  set(level 3)
endif(NOT runtime_error)

set(level_ok FALSE)
if (level EQUAL expected_level)
  set(level_ok TRUE)
endif()

set(output_error "Not compared")
if (level_ok)
  message ("** Comparing test output with expected one.")

  set(output_file  '')
  if(EXISTS ${WORK_DIR}/output.txt)
    set(output_file ${WORK_DIR}/output.txt)
  elseif(EXISTS ${WORK_DIR}/stdout.txt)
    set(output_file ${WORK_DIR}/stdout.txt)
  endif()

  execute_process(
    COMMAND ${CMAKE_COMMAND} -E compare_files ${output_file} ${SOURCE_DIR}/expected.txt
    RESULT_VARIABLE output_error)
endif()

message("******************************************************") 
message("            Test results")
message("            ------------")
message("Compiling and linking: ${compile_error}")
message("Runnig program:        ${runtime_error}")
message("Expected output:       ${output_error}\n")


if(output_error)
  message ("Test status: FAILED.")
  message ("Aditional Information:")
  if(level_ok)
    message("  Output different from expected.")
    message("  Run something like the following for further details:\n")
    message("  tkdiff ${output_file} ${SOURCE_DIR}/expected.txt\n")
  endif()
  if(compile_error)
    message("  Compile or linking Error. See the test output above for details.\n")
  elseif(runtime_error)
    message("  Run Time Error.")
    message("  You can try to run the executable with:\n   cd ${WORK_DIR} & ${TEST_PROG}\n")
  endif()
  message(FATAL_ERROR " ")
endif()

message ("Test status: PASSED.")
message ("******************************************************") 
