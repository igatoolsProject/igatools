
```
sudo port install gcc49 gcc_select 
sudo port select gcc mp-gcc49
```
  1. Install git
```
sudo port install git
```
  1. Install boost
```
sudo port install boost
```
  1. Install doxygen
```
sudo port install doxygen
```
  1. Install cmake
```
sudo port install cmake
```
  1. Install TBB
```
sudo port install tbb
```
  1. Make sure that the CC and CXX environmnet variables  are defined as suggested in [wiki:InstallingIgatools#dirs].


## Compilation problems using g++ 4.8 with MacOSX 10.9 (Mavericks) ##
If you encounter compilation problems using g++ 4.8 on MacOSX 10.9 (Mavericks) maybe this link can help you http://stackoverflow.com/questions/19649421/something-odd-happened-to-c-11-in-mavericks