# vinecopulib
[![Build Status](https://github.com/vinecopulib/vinecopulib/workflows/Build%20Status/badge.svg)](https://github.com/vinecopulib/vinecopulib/actions)
[![Coverage Status](https://img.shields.io/codecov/c/github/vinecopulib/vinecopulib/master.svg)](https://codecov.io/github/vinecopulib/vinecopulib?branch=master)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

#### What are vine copulas?

Vine copulas are a flexible class of dependence models consisting of bivariate
building blocks (see e.g.,
[Aas et al., 2009](https://mediatum.ub.tum.de/doc/1083600/1083600.pdf)).
You can find a comprehensive list of publications and other materials on
[vine-copula.org](http://www.statistics.ma.tum.de/en/research/vine-copula-models/).

#### What is vinecopulib?

vinecopulib is a header-only C++ library for vine copula models based on
[Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page). It provides
high-performance implementations of the core features of the popular
[VineCopula R library](https://github.com/tnagler/VineCopula), in particular
inference algorithms for both vine copula and bivariate copula models.
Advantages over VineCopula are  
* a stand-alone C++ library with interfaces to both R and Python,
* a sleaker and more modern API,
* shorter runtimes and lower memory consumption, especially in high dimensions,
* nonparametric and multi-parameter families.

#### Status

Version [0.3.2](https://github.com/vinecopulib/vinecopulib/releases) was
released on July 3, 2019. While we did our best to
design a user-friendly API, the library is still under active development and
changes are to be expected. We are also working on interfaces for
[R](https://github.com/vinecopulib/rvinecopulib) and
[Python](https://github.com/vinecopulib/pyvinecopulib).

#### Contact

If you have any questions regarding the library, feel free to
[open an issue](https://github.com/vinecopulib/vinecopulib/issues/new) or
send a mail to <info@vinecopulib.org>.

#### Documentation

For documentation of the library's functionality and
instructions how to use it, check out our
[website](https://vinecopulib.github.io/vinecopulib/) or the `docs/` folder
in this repository.
