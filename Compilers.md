igatools is **heavily** based on C++11, so you need a compiler that implements the C++11 standard as much as possible.

Usually we develop igatools on Linux and Mac OS X and so far we have the following information:

### Linux environments ###

The following compilers are tested and fully working with igatools:
  * GNU gcc-4.8.1
  * GNU gcc-4.8.2
  * GNU gcc-4.9.1

The following compilers do not (currently) work with igatools (mainly due to some missing C++11 features in the compiler):
  * Intel C++ compilers with version <= 14
  * clang compilers with version <= 3.3

### Mac OS X environments ###

The following compilers are tested and fully working with igatools:
  * GNU gcc-4.8.2 (using [MacPorts](http://www.macports.org))
  * GNU gcc-4.9.1

### Windows environments ###
Not tested