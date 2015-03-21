## Installing igatools ##
Currently, to install igatools you need to compile it from its sources. Below there are some explanatory steps to guide you through the process.

  * In order to get, configure and compile the sources you need to have the following programs available in your system:
    * mandatory:
      * [git](http://www.git-scm.com) (type `git --version`)
      * [cmake](http://www.cmake.org) (type `cmake --version`)
      * [g++](http://gcc.gnu.org) version >= 4.9.1 (you can type `g++ --version`)
    * optional (if you want to build the documentation):
      * [doxygen](http://www.doxygen.org) (type `doxygen --version`)
      * [graphviz](http://www.graphviz.org) (type `dot -V`)
> If you are a Linux user we assume you know how to install the  appropriate packages for your distribution.
> As we have experienced more difficulties making it work on Mac OSX, we provide some additional instructions [here](MacInstallationNotes.md)
  * Additionally, igatools relies on some libraries that also need to be available in your system. The main ones are listed below, in particular we provide additional instruction for installing Trilinos.
    * mandatory:
      * [Boost](http://www.boost.org)
      * [Trilinos](http://trilinos.sandia.gov) (some installation instructions [here](InstallingTrilinos.md))
    * optional:
      * [Intel Threading Building Blocks (TBB)](https://www.threadingbuildingblocks.org)

For clarity we will assume that you will be installing igatools as a common user (as opposed to root) somewhere in your home directory. More precisely, we will use the following directory structure for the sources, build and install of the library and its related dependencies.
| ~/local/usr/include  | location of the headers |
|:---------------------|:------------------------|
| ~/local/usr/lib      | location of the library |
| ~/local/usr/doc      | location of the documentation |
| ~/local/src          | where the library sources are |
| ~/local/workspace   | where we build the library |

If you agree with this directory structure just copy and paste the following into a terminal
```
cd ~
mkdir local
mkdir local/usr
mkdir local/src
mkdir local/workspace
```

In all the steps below we assume that the following environment variables are defined with your system paths. They are usually defined at the start-up of the shell
|IGATOOLS\_SOURCE      | igatools sources |
|:---------------------|:-----------------|
|IGATOOLS\_WS          | igatools workspace |
|IGATOOLS\_PREFIX      | igatools installation |
|TRILINOS\_PREFIX | trilinos installation |
|TRILINOS\_LIB    | trilinos library |
|IGATOOLS\_LIB         | igatools libray |
|ECLIPSE\_EXE     | eclipse |

We provide cut and paste for the lazy (they assume you are using the Bash shell):
```
echo "
IGATOOLS_SOURCE=~/local/src/igatools
IGATOOLS_PREFIX=~/local/usr
IGATOOLS_WS=~/local/workspace
TRILINOS_PREFIX=~/local/usr/trilinos
TRILINOS_LIB=\$TRILINOS_PREFIX/lib
IGATOOLS_LIB=\$IGATOOLS_PREFIX/lib
ECLIPSE_EXE=~/local/usr/eclipse/eclipse" >> ~/.bash_login
echo "export IGATOOLS_SOURCE IGATOOLS_WS IGATOOLS_PREFIX TRILINOS_PREFIX TRILINOS_LIB IGATOOLS_LIB ECLIPSE_EXE" >> ~/.bash_login
source ~/.bash_login
```
If you use Mac OSX also do
```
echo "  
export CC=/opt/local/bin/gcc-mp-4.8
export CXX=/opt/local/bin/g++-mp-4.8" >> ~/.bash_login
source ~/.bash_login
```
In a Linux system you would most likely use .bashrc instead of .bash\_login.

  * Get the igatools sources from the git repository
```
cd ~/local/src/
git clone https://code.google.com/p/igatools/
cd igatools
git checkout develop
```
  * Configure the library
```
cd $IGATOOLS_WS
mkdir igatools_lib
cd igatools_lib
cmake $IGATOOLS_SOURCE
```
  * At the end of a successful cmake configuration, is a good idea to export the igatools path to the library path. This depends on your system (linux or mac) and shell. For Bash shell you need to add to your .bash\_login or .bashrc
    * for Linux
```
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$TRILINOS_LIB:$IGATOOLS_LIB # Bash on Linux
```
    * for MacOSX
```
export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:$TRILINOS_LIB:$IGATOOLS_LIB # Bash on MacOSX
```
  * Finally, let build and install the library. Remember to replace `N` below by the number of threads you want to use
```
make -jN install #install the library
make doc         #install the documentation
```

At this point, if everything is gone well you should have the igatools library in `$IGATOOLS_LIB` and compiled in **Debug** mode that is safe but slow. If you want to attain maximum performance, see below.

To get started with the library go [here](GettingStarted.md).

## About error checks, performance and debugging (a. k. a. "compilation modes") ##
The library makes use **_defensive programming_** techniques
to detect an easily find bugs at runtime, through the
exception handling mechanism.
We adopt two levels of checks, depending on the _compilation mode_ you have chosen to use.

  * The first level of checks, is active only when igatools is compiled in **"Debug"** mode and it is expensive in terms of running time and size of the compiled library. We perform this kind of checks wherever there may be a chance of error,  e.g. out-of-bound index for accessing to vectors elements,  uninitialized objects, invalid object states, mismatching dimensions, etc.
  * The second level of checks is always active (both in **"Debug"** and in **"Release"** mode) and it is used for checking anomalies that may be introduced by the input data.

The typical workflow for a user of igatools would be to first write his code and test it with the library  in **"Debug"** mode on a small-size problem. When this is working as expected and (virtually) bug-free, link the code with igatools compiled in **"Release"** mode on a real-size problem.


By _default_, igatools will be compiled in **"Debug"** mode, if you want to attain maximum performance (switching off a lot of error checks), you must compile the library in **"Release"** mode. In order to do so, you should add to the cmake command, the option `-DCMAKE_BUILD_TYPE=Release`.



## Remark ##
The igatools repository is in continuous development.
In order to be up to date with the current development you should
periodically update your sources.
```
cd $IGATOOLS_SOURCE
git pull origin develop
cd $IGATOOLS_WS/igatools_lib
cmake $IGATOOLS_SOURCE
make -jN install
make doc
```


For advanced users we include [some custom configurations](CustomConfigurations.md)