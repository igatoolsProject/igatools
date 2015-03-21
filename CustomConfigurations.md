### Configuration examples (running cmake) ###
Below you will find some configuration examples used
in different common scenarios.
They should provide a reasonable idea of what is available, you may need to mix and match to fit your own scenario.



Configuration examples:

  * If you want to use Eclipse to develop with the library
```
cmake $IGATOOLS_SOURCE -G"Eclipse CDT4 - Unix Makefiles" -DCMAKE_ECLIPSE_EXECUTABLE=$ECLIPSE_EXE
```
  * Specify the install location
```
cmake $IGATOOLS_SOURCE -DCMAKE_INSTALL_PREFIX=$IGATOOLS_PREFIX
```
  * Specify Trilinos location
```
cmake $IGATOOLS_SOURCE -DCMAKE_INSTALL_PREFIX=$IGATOOLS_PREFIX -DTrilinos_PREFIX=$TRILINOS_PREFIX
```
  * Using quadruple precision **(EXPERIMENTAL)**
```
cmake $IGATOOLS_SOURCE -DCMAKE_INSTALL_PREFIX=$IGATOOLS_PREFIX -DQuadPrecision=ON
```