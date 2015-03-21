## Installing Trilinos ##

Igatools uses [Trilinos](http://trilinos.org) for its linear algebra (matrices, vectors and solvers) capabilities.

Trilinos is a big library (and may not be easy to install), as a matter of facts it does not provide binary packages. You need to compile the library yourself. Below we provide an example that we have successfully used with our library.

In the instructions below we assume the directory layout suggested in the [igatools installation guide](InstallingIgatools.md).
  * Download the [Trilinos](http://trilinos.org/?page_id=85) sources  (version required >= 11.6.1) and save it into the directory `~/local/src`
  * You may need to install the **blas** and **lapack** libraries if not already available in you system.
  * Unpack the tarball and create a build directory (replace `???` with the Trilinos version you are installing)
```
cd ~/local/src
tar -xvzf trilinos-???-Source.tar.gz
cd trilinos-???
mkdir build
cd build
```
  * Assuming you have [defined](InstallingIgatools.md) the path for the Trilinos installation directory, configure Trilinos (by running cmake).
```
cmake \
  -DCMAKE_INSTALL_PREFIX:PATH=$TRILINOS_PREFIX \
  -DCMAKE_BUILD_TYPE:STRING=RELEASE \
  -DBUILD_SHARED_LIBS:BOOL=ON \
  -DTrilinos_ENABLE_Fortran:BOOL=OFF \
  -DTeuchos_ENABLE_COMPLEX:BOOL=ON \
  -DEpetra_ENABLE_THREADS:BOOL=TRUE \
  -DTrilinos_ENABLE_Amesos:BOOL=ON \
  -DAmesos_ENABLE_LAPACK:BOOL=ON \
  -DTrilinos_ENABLE_Amesos2:BOOL=ON \
  -DAmesos2_ENABLE_LAPACK:BOOL=ON \
  -DTrilinos_ENABLE_Belos:BOOL=ON \
  -DTrilinos_ENABLE_Ifpack:BOOL=ON \
  -DTrilinos_ENABLE_Ifpack2:BOOL=ON \
  -DTrilinos_ENABLE_Teuchos:BOOL=ON \
  -DTrilinos_ENABLE_Tpetra:BOOL=ON \
  -DThyra_ENABLE_EXPLICIT_INSTANTIATION:BOOL=ON \
  -DTPL_ENABLE_SuperLU:BOOL=ON \
  ..
```
  * After configuring the building system (i.e. running cmake), compile and install Trilinos
```
# replace N with the number of cores you want to use
make -jN install
```