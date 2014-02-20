#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Find sed (Required)
#-------------------------------------------------------------------------------
macro(find_sed)
  find_program(SED_EXECUTABLE sed)
  message(STATUS "Found sed:  ${SED_EXECUTABLE}.") 
endmacro(find_sed)
