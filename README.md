# vinecopulib
[![Build Status](https://travis-ci.org/tvatter/vinecoplib.svg?branch=master)](https://travis-ci.org/tvatter/vinecopulib)

A C++ library for vine copulas.

## Requirements

cmake, GSL, Eigen3, Rcpp, RInside, RcppEigen, VineCopula

Note that Rcpp, RInside RcppEigen and VineCopula are used only for unit-testing the C++ implementation (i.e., comparing with the results from the VineCopula R package).
A findR.cmake looks for them in the default locations for linux and osx, but problems might occur
with versions installed from R/RStudio. Therefore, prior to building the library,
it is recommended to use:

`sudo Rscript -e 'install.packages(c("Rcpp","RInside","RcppEigen","VineCopula"), repos="http://cran.rstudio.com/")'`

## Compilation

* Create a build folder:

`mkdir build`

* Move to the created folder:

`cd build`

* Create the `MakeFile` via cmake:

`cmake ..`

* Compile the code to generate the executable:

`make` or `make -j n` where `n` is the number of cores to use for the compilation

* A folder bin will be created at the same level as the build directory.
