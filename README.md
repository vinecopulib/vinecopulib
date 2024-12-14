# vinecopulib

[![Build Status](https://github.com/vinecopulib/vinecopulib/workflows/Build%20Status/badge.svg)](https://github.com/vinecopulib/vinecopulib/actions)
[![Coverage Status](https://img.shields.io/codecov/c/github/vinecopulib/vinecopulib/main.svg)](https://codecov.io/github/vinecopulib/vinecopulib?branch=main)
[![Codacy Badge](https://app.codacy.com/project/badge/Grade/2c8d45ebcb954082b409b4a2bd31af2b)](https://www.codacy.com/gh/vinecopulib/vinecopulib/dashboard?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=vinecopulib/vinecopulib&amp;utm_campaign=Badge_Grade)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Documentation](https://img.shields.io/website/http/vinecopulib.github.io/vinecopulib.svg)](https://vinecopulib.github.io/vinecopulib/)
[![DOI](https://zenodo.org/badge/76354683.svg)](https://zenodo.org/badge/latestdoi/76354683)

#### What is vinecopulib?

[vinecopulib](https://vinecopulib.github.io/vinecopulib/) is a header-only C++ library for vine copula models based on
[Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page). It provides
high-performance implementations of the core features of the popular
[VineCopula R library](https://github.com/tnagler/VineCopula), in particular
inference algorithms for both vine copula and bivariate copula models.
Advantages over VineCopula are  

- a stand-alone C++ library with interfaces to both R and Python,
- a sleaker and more modern API,
- shorter runtimes and lower memory consumption, especially in high dimensions,
- nonparametric and multiparameter families.

The library also has interfaces for
[R](https://github.com/vinecopulib/rvinecopulib) and
[Python](https://github.com/vinecopulib/pyvinecopulib).

#### Contact

If you have any questions regarding the library, feel free to
[open an issue](https://github.com/vinecopulib/vinecopulib/issues/new) or
send a mail to [info@vinecopulib.org](mailto:info@vinecopulib.org).

#### Documentation

For documentation of the library's functionality and
instructions how to use it, check out our
[website](https://vinecopulib.github.io/vinecopulib/) or the `docs/` folder
in this repository.
