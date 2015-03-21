# igatools documentation & manual #

If you want to use igatools, we suggest (before embarking in the study of the igatools' classes documentation, as explained below) to read [this report](https://www.dropbox.com/s/jhraxy2krvs83sz/3PV14-1-0.pdf) for a general design overview and some application examples.


### Manually generated documentation ###

The igatools documentation and manual can be generated using [Doxygen](http://www.doxygen.org) from the directory you have chosen to use for building the library (if you have follwed the [installation instructions](InstallingIgatools.md), it is the directory contained assigned to the environment variable `$IGATOOLS_WS`).

First of all you must move into the igatools building directory:
```
cd $IGATOOLS_WS
```
then, if you have not done yet, you should invoke the cmake command before the compilation and installation phase (as explained in the [installation instructions](InstallingIgatools.md)).

Finally you can invoke the documentation generation, by using the command
```
make doc
```

The igatools manual will be created as html and the starting page can be visualized by any browser at the location
```
$IGATOOLS_PREFIX/doc/html/index.html
```


### Online documentation ###

If you have an active internet connection, it is also possible to get the documentation [here](http://www.igatools.altervista.org).