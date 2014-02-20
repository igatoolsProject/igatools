#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Find Phython library (Required)
#-------------------------------------------------------------------------------
macro(find_python)
  include(FindPythonInterp)
  find_package(PythonInterp REQUIRED)
endmacro(find_python)
