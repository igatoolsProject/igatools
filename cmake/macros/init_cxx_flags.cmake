#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Compiler and compile flags
macro(init_cxx_flags)
  # Initialize CXXFLAGS.
  set(CMAKE_CXX_FLAGS                "-Wall -O0 -g")
  set(CMAKE_CXX_FLAGS_DEBUG          "-Wall -O0 -g")
  set(CMAKE_CXX_FLAGS_MINSIZEREL     "-Os -DNDEBUG")
  set(CMAKE_CXX_FLAGS_RELEASE        "-O3 -mtune=native -DNDEBUG")
  set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O3 -mtune=native -DNDEBUG -g")

  # Compiler-specific C++11 activation.
  if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "GNU")
    execute_process(
      COMMAND ${CMAKE_CXX_COMPILER} -dumpversion OUTPUT_VARIABLE GCC_VERSION)
    if (NOT (GCC_VERSION VERSION_GREATER 4.8))
      message(FATAL_ERROR "gcc version 4.8 or greater.")
    endif()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
  elseif ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    message(FATAL_ERROR "Your C++ CLANG compiler is not supported by igatools")
  elseif("${CMAKE_CXX_COMPILER_ID}" MATCHES "Intel")
    message(FATAL_ERROR "Your C++ INTEL compiler is not supported by igatools")
  else()
    message(FATAL_ERROR "Your C++ compiler is not supported by igatools")
  endif()
  set(LINKER_FLAGS "-rdynamic")
endmacro(init_cxx_flags)