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
      set (target_name inst-${dir}-${name}) 
      add_custom_command(OUTPUT ${inst_file}
	COMMAND PYTHONPATH=${PROJECT_SOURCE_DIR}/cmake/instantiation_scripts
        ${PYTHON_EXECUTABLE} -B ${py_file} out_file=${name}.inst 
        config_file=${CMAKE_CURRENT_BINARY_DIR}/instantiation_table.txt
        max_der_order=${max_der_order}
	DEPENDS ${py_file}  
        ${inst_script} 
        ${CMAKE_CURRENT_BINARY_DIR}/instantiation_table.txt
        ${PROJECT_SOURCE_DIR}/cmake/instantiation_scripts/generate_inst_table.py
	WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/include/${lib_name}/${dir}
	COMMENT "Generating file ${dir}/${name}.inst")
      set_property(SOURCE ${PROJECT_SOURCE_DIR}/source/${dir}/${name}.cpp
	PROPERTY OBJECT_DEPENDS ${inst_file})
    endforeach() 
  endforeach()
endmacro(generate_instantiations)

