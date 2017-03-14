# vinecopulib
[![Build Status](https://travis-ci.org/tvatter/vinecopulib.svg?branch=master)](https://travis-ci.org/tvatter/vinecopulib)
[![Windows Build status](http://ci.appveyor.com/api/projects/status/github/tvatter/vinecopulib?svg=true)](https://ci.appveyor.com/project/tvatter/vinecopulib)
[![Coverage Status](https://img.shields.io/codecov/c/github/tvatter/vinecopulib/master.svg)](https://codecov.io/github/tvatter/vinecopulib?branch=master)

A C++ library for vine copulas.

## Requirements

cmake, doxygen, GSL, Eigen3, VineCopula

Note that VineCopula is used only for unit-testing the C++ implementation (i.e., comparing with the results from the R package).
A findR.cmake looks for them in the default locations for linux and osx, but problems might occur
with versions installed from R/RStudio. Therefore, prior to building the library,
it is recommended to use:

`sudo Rscript -e 'install.packages("VineCopula", lib="/usr/lib/R/library", repos="http://cran.rstudio.com/")'`

## Compilation

* Create a build folder:

`mkdir build`

* Move to the created folder:

`cd build`

* Create the `MakeFile` via cmake:

    * `cmake .. -DCMAKE_BUILD_TYPE=Debug` for the debug version
    * `cmake .. -DCMAKE_BUILD_TYPE=Release` for the release version (faster)

* Compile the code to generate the executable:

`make` or `make -j n` where `n` is the number of cores to use for the compilation

* A folder bin will be created at the same level as the build directory.

* To build the documentation run
`make doc`
