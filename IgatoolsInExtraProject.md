Once you have compiled and installed igatools you can integrate it into your project


## Simple programs ##
For creating simple projects (like the tutorial examples) using igatools,
you can use the following configuration instruction.
It assumes all project files are in a single directory and the code in a single .cpp file.

  1. Create a directory for your project
```
mkdir my_project_dir
cd my_project_dir
```
  1. Create a text file CMakeLists.txt with the following content
```
cmake_minimum_required(VERSION 2.8.10)
set (name your_program_name) 
project(${name})

set(CMAKE_PREFIX_PATH $ENV{IGATOOLS_LIB} ${CMAKE_PREFIX_PATH} )
find_package(igatools REQUIRED)

set(CMAKE_CXX_FLAGS ${IGATOOLS_CXX_FLAGS})

include_directories(${IGATOOLS_INCLUDE_DIRS})
link_directories(${IGATOOLS_LIBRARY_DIR})

add_executable(${name} ${name}.cpp )

target_link_libraries (${name} ${IGATOOLS_LIBRARIES})

add_custom_target(run-${name} ./${name} DEPENDS ${name}
     COMMENT "Running program ${name}")
```
> > replacing the "`your_program_name`" string by the name you want to give to your program.
  1. Create a source file `your_program_name.cpp`
  1. Configue
```
cmake .  # don't forget the "dot" at the end
```

At this point you can edit compile and run your program
  1. To compile and link
```
make
```
  1. To execute
```
./your_program_name
```


---


## Advance projects ##
If you are developing a more advance project we suggest that you create a
CMakeLists.txt adapted for your project.
igatools provides a cmake module to detect itself, that upon success it makes
the following variables available:

**IGATOOLS\_VERSION**

> is the version of the library expressed with the format x.y.z where x is the major version number, y is the minor version number and z is the patch version number.

**IGATOOLS\_CXX\_COMPILER**
> is the C++11 compiler used for building the library.

**IGATOOLS\_CXX\_FLAGS**
> are the C++11 compiler flags used for building the library.

**IGATOOLS\_LINKER\_FLAGS**
> are the linker flags used for building the library.

**IGATOOLS\_INCLUDE\_DIRS**
> are the headers directories needed by the library.

**IGATOOLS\_LIBRARY\_DIR**
> is the list the directories containing the libraries needed for the linking phase.

**IGATOOLS\_LIBRARIES**
> is the list of libraries needed for the use of the current library in an external application.