SNCTRL is an optimal control interface for the nonlinear optimization package SNOPT.

SNOPT is required to use SNCTRL, but is not included here.


Configure and install SNCTRL:
```
./configure  --with-snopt=/path/to/SNOPT/library
make
make install
```

8 examples (coded for each of the 3 interfaces) are included in the package.
To build the included examples:
```
make examples
```

Assuming no modifcations to the location of files and libraries, an example can
be run with the following commands
```
cd $SNCTRL/examples
./rocketA
```

