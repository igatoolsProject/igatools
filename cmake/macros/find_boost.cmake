#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Find Boost library (Required)
#-------------------------------------------------------------------------------
macro(find_boost)
  find_package(Boost 1.48.0 REQUIRED)
  include_directories(${Boost_INCLUDE_DIRS})
endmacro(find_boost)
