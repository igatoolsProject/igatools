#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Creates a file with a table of the spaces that the library
# is suppossed to be compiled for.
#
macro(create_space_table)
  message(STATUS "Generating physical space table.")
  execute_process(COMMAND 
    ${PYTHON_EXECUTABLE} 
    -B ${PROJECT_SOURCE_DIR}/cmake/instantiation_scripts/generate_inst_table.py 
    dim_ref=${dim}
    range_ref=${range}
    rank_ref=${rank}
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    RESULT_VARIABLE res) 
endmacro(create_space_table)