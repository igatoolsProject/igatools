## Igatools coding style ##

Conventional wisdom says that no particular coding style is better than any other, but it is important to have one. It is easier to read and maintain a code base which uses a consistent style.

We strive to keep our programming style and the kind of interfaces we provide as consistent as possible.
To this end, we have adopted a set of coding conventions that we attempt to follow wherever possible.
They have two parts: style issues, and something we call "defensive programming", the latter being an attempt to let our code help us find bugs.
The purpose of the coding style is to keep the library as uniform as possible.
Uniformity reduces the number of bugs we produce because we can, for example, always assume that input arguments come before output arguments of a function call.
They also simplify reading code because some things become clear already by looking  at the style a piece of code is written, without having to look up the exact definition of something.

We use astyle to enforce some of the style guidelines.

We also follow most of the  [Google C++ Style Guide](http://google-styleguide.googlecode.com/svn/trunk/cppguide.xml) with a few exceptions.
Not because we like all the conventions but because there is python script (cpplint.py) that detects violation automatically.
For sake of self consistency we rewrite below what we consider the most important ones.

## Style guides ##
### Function Parameter Ordering ###
  * When defining a function, parameter order is: inputs, then outputs.
  * Avoid a parameter that is both input and output
  * If it does not worse efficiency return by value (move constructor assumed)


### Namespaces ###
  * You may use a using-declaration anywhere in a .cpp file, and in functions, methods or classes in .h files.
  * No using in .h otherwise.

### Source documentation and comments ###
  * Doxygen comments should only be in .h files (not in the .cpp)
  * Avoid comments in between the lines of code, put instead a readable explanation at the beginning. Try to let your code be self documenting.

  * Do not be unnecessarily verbose or state the completely obvious. Notice below that it is not necessary to say "returns false otherwise" because this is implied.
```
  // Returns true if the table cannot hold any more entries.
  bool IsTableFull();
```
  * Do not use `@briefs` or `@param` to document functions, the brief of the function should be implied by its name and parameters names.
  * We will use `@brief` for classes but only after the library itself and its documentation are reasonably stable


### Naming ###

**types** (typedefs and using):
  * prefer the C++11 keyword `using` over the keyword `typedef`;
  * using internal type alias is encourage for long types;
  * for public aliases use CamelCase (ex `Derivatives`);
  * for private alias lower case followed by _t (ex `return_t`);_

**classes**
  * Classes names are camel case (ex `CartesianGrid`)
  * Put one class per file except to classes very similar in name and meaning

**functions**
  * function and member function in lower case
  * private and protected member variables end with an underscore (ex `private_member_`)
  * variables are lowercase

**files**
  * file names are in lowercase with words separated with the underscore (_)
  * file extensions are: .h ; -inline.h ;  .cpp and .inst.py_

### Spaces and new lines ###
We try not to use unnecessary new lines, one case where we require new lines is to separate the implementation of member  functions (in the .cpp file) in this case leave 3 new lines beetwen them.
In the header file only 1 (one) new line should be left to separate
members.

### print\_info ###
All (most) igatools classes have a
` void print_info(LogStream & out) const `
This function is for debugging and tests purposes so some style must be followed
  1. only show the private variable by calling the corresponding print\_info
  1. Do not call class functions to generate output
  1. use out.begin\_item(), out.end\_item() to separate showing variables
  1. Do not by any means show information that may differ in different computers such as memory addresss
  1. Do not call print\_info for pointers (only stack variables)
  1. Do not show the name of the object type


## Style tools ##

### [astyle](http://astyle.sourceforge.net/) ###

Artistic Style is a source code indenter, formatter, and beautifier.
For the igatools code you need to use a version >= 2.03.
The style file can be found in the util directory of the sources.

In order to use astyle with the coding style chosen for the igatools libray, a git's pre-commit hook should be used from any client that is pushing some commits.

Unfortunately, the `git clone` command cannot setup the client-side hooks, so it must be done by hand by any igatools developer on their machine.

For this purpose we provide the proper hook script containing the command necessary to run astyle and apply the style configuration we adopt in igatools.

In order to use this pre-commit hook you should perform the following tasks:
  1. Make sure you have a working version of the astyle program
  1. Manually link the file `$IGATOOLS_SOURCE/utils/pre-commit.astyle` to `$IGATOOLS_SOURCE/.git/hooks/pre-commit`. For the lazy, here the code to perform the previous action with Linux or MacOSX
```
ln -s $IGATOOLS_SOURCE/utils/pre-commit.astyle $IGATOOLS_SOURCE/.git/hooks/pre-commit
```
  1. use git as usual.

**Note:** the hook works in a way such that astyle runs using the files staged for the commit. If one or more files are changed by astyle, then you should run `git commit` a second time without any change in the files between the first and second commit. You can always check which file is modified by astyle (and then if you need to run `git commit` a second time) with the command `git status`.