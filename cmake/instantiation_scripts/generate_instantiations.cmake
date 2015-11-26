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

# cmake macro to generate the template instantiations
# We look for all *.inst.py files in the source directories
# and we execute them to generate the *.inst in the build/include/...

macro(generate_instantiations)
  
  set(inst_script 
  ${PROJECT_SOURCE_DIR}/cmake/instantiation_scripts/init_instantiation_data.py)
  
  foreach(dir ${source_dirs})
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include/${lib_name}/${dir})
    file(GLOB files ${PROJECT_SOURCE_DIR}/source/${dir}/*.inst.py)
    foreach(py_file ${files})
      get_filename_component(name ${py_file} NAME_WE)
      set (inst_file ${CMAKE_CURRENT_BINARY_DIR}/include/${lib_name}/${dir}/${name}.inst)
      add_custom_command(
        OUTPUT ${inst_file}
        COMMAND PYTHONPATH=${PROJECT_SOURCE_DIR}/cmake/instantiation_scripts
        ${PYTHON_EXECUTABLE} -B ${py_file} out_file=${name}.inst 
        config_file=${CMAKE_CURRENT_BINARY_DIR}/instantiation_table.txt
        max_der_order=${max_der_order}
        nurbs=${NURBS}
        DEPENDS ${py_file}  
        ${inst_script} 
        ${CMAKE_CURRENT_BINARY_DIR}/instantiation_table.txt
        ${PROJECT_SOURCE_DIR}/cmake/instantiation_scripts/generate_inst_table.py
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include/${lib_name}/${dir}
        COMMENT "Generating file ${dir}/${name}.inst")
        set_property(SOURCE ${PROJECT_SOURCE_DIR}/source/${dir}/${name}.cpp 
        PROPERTY OBJECT_DEPENDS ${inst_file})
        execute_process(COMMAND ${SED_EXECUTABLE} -e  \\,.inst,d
        INPUT_FILE  ${PROJECT_SOURCE_DIR}/source/${dir}/${name}.cpp
        OUTPUT_FILE ${CMAKE_CURRENT_BINARY_DIR}/include/igatools/${dir}/${name}-template.h)
    endforeach(py_file ${files}) 	
    endforeach(dir ${source_dirs})  
    
endmacro(generate_instantiations)



# The programmer must provide a file_name.serial.py in the source directory of the library
macro(generate_serial)
  
  set(inst_script 
  ${PROJECT_SOURCE_DIR}/cmake/instantiation_scripts/init_instantiation_data.py)
  
  foreach(dir ${source_dirs})
    file(GLOB files ${PROJECT_SOURCE_DIR}/source/${dir}/*.serial.py)
    foreach(py_file ${files})
      get_filename_component(name ${py_file} NAME_WE)
      set (serial_file ${CMAKE_CURRENT_BINARY_DIR}/include/${lib_name}/${dir}/${name}.serial)
      add_custom_command(
        OUTPUT ${serial_file}
        COMMAND PYTHONPATH=${PROJECT_SOURCE_DIR}/cmake/instantiation_scripts
        ${PYTHON_EXECUTABLE} -B ${py_file} out_file=${name}.serial
        config_file=${CMAKE_CURRENT_BINARY_DIR}/instantiation_table.txt
        max_der_order=${max_der_order}
        nurbs=${NURBS}
        DEPENDS ${py_file}  
        ${inst_script} 
        ${CMAKE_CURRENT_BINARY_DIR}/instantiation_table.txt
        ${PROJECT_SOURCE_DIR}/cmake/instantiation_scripts/generate_inst_table.py
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include/${lib_name}/${dir}
        COMMENT "Generating file ${dir}/${name}.serial")
        set_property(SOURCE ${PROJECT_SOURCE_DIR}/source/${dir}/${name}.cpp 
        APPEND PROPERTY OBJECT_DEPENDS ${serial_file})
    endforeach(py_file ${files}) 	
    endforeach(dir ${source_dirs})  
    
endmacro(generate_serial)