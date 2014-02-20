#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Find Trilinos library (Required)
#-------------------------------------------------------------------------------
macro(find_trilinos)
  set(Trilinos_PREFIX $ENV{TRILINOS_PREFIX} CACHE LOCATION 
    "Location where Trilinos library is installed")
  find_package(Trilinos 11 REQUIRED PATHS ${Trilinos_PREFIX})
  message(STATUS "Found Trilinos:  version ${Trilinos_VERSION}.")
  if (NOT (Trilinos_VERSION VERSION_GREATER 11.6))
    message(FATAL_ERROR "Trilinos 11.6.1 or greater is required.")
  endif()
  # Check that individual required Trilinos packages are available
  set(tri_required_packages Tpetra Belos)
  foreach(package ${tri_required_packages})
    list(FIND Trilinos_PACKAGE_LIST  ${package} package_index)
    if (package_index EQUAL -1)
      message(FATAL_ERROR "Trilinos ${package} package not found.")
    endif()
  endforeach(package)
endmacro(find_trilinos)
