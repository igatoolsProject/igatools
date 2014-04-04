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

# This is a cmake script that runs a test.
# 
# It runs the following steps:
# 1 - compiles and links the test program
# 2 - execute the program in it working directory
# 3 - compares the output of the program with the expected output
#
# 
# It's called from CMakeLists.txt as follows:
# cmake -DTEST_PROG=${CMAKE_CURRENT_BINARY_DIR}/${dir}/${name}.exe
#       -DSOURCEDIR=${CMAKE_SOURCE_DIR}/${dir}/${name}
#       -DWORK_DIR=${CMAKE_CURRENT_BINARY_DIR}/${dir}/${name}
#       -P ${CMAKE_SOURCE_DIR}/runtest.cmake
#
# TODO: make should be CMAKE_MAKE_COMMAND but is empty!!!

message ("******************************************************") 


message ("** Compiling and linking test program")
execute_process(COMMAND make ${EXE_TARGET}	RESULT_VARIABLE compile_error)
if(compile_error)
  message(FATAL_ERROR "** Test FAILED with a Compile Error.")
endif(compile_error)

 
message ("** Executing test program.")
if(NOT EXISTS ${WORK_DIR})
    message(FATAL_ERROR "** Work directory does not exist.")
endif()

execute_process(COMMAND ${TEST_PROG} 
                OUTPUT_FILE stdout.txt
                ERROR_FILE stdout.txt
                WORKING_DIRECTORY ${WORK_DIR} 
    	        RESULT_VARIABLE runtime_error)
if(runtime_error)
    message(FATAL_ERROR "** Test FAILED with a Run Time Error.\n")
endif()


message ("** Comparing test output with expected one ...")
if(NOT EXISTS ${SOURCE_DIR}/expected.txt)    
    message(FATAL_ERROR "** expected.txt does not exist.")
endif()

set(output_file  '')
if(EXISTS ${WORK_DIR}/output.txt)
	set(output_file ${WORK_DIR}/output.txt)
elseif	(EXISTS ${WORK_DIR}/stdout.txt)
	set(output_file ${WORK_DIR}/stdout.txt)
endif()

execute_process(COMMAND ${CMAKE_COMMAND} -E compare_files
    ${output_file} ${SOURCE_DIR}/expected.txt
	RESULT_VARIABLE output_error)
if(output_error)
    message("** Test FAILED with output different from expected.")
    message("** Run something like the following for further details:\n")
    message("tkdiff ${output_file} ${SOURCE_DIR}/expected.txt\n")
    message(FATAL_ERROR " ")
endif()
message ("** Ouputs are similars. Test PASSED.")


message ("******************************************************") 
