# vinecopulib
[![Build Status](https://travis-ci.org/tvatter/vinecopulib.svg?branch=master)](https://travis-ci.org/tvatter/vinecopulib)
[![Windows Build status](http://ci.appveyor.com/api/projects/status/github/tvatter/vinecopulib?svg=true)](https://ci.appveyor.com/project/tvatter/vinecopulib)
[![Coverage Status](https://img.shields.io/codecov/c/github/tvatter/vinecopulib/master.svg)](https://codecov.io/github/tvatter/vinecopulib?branch=master)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT) 
 
#### What are vine copulas?

Vine copulas are a flexible class of dependence models consisting of bivariate
building blocks (see e.g., Aas et al., 2009). You can find a comprehensive list
of publications and other materials on vine-copula.org.

#### What is vinecopulib?

vinecopulib is a C++ library for vine copula models based on 
[Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page). It provides 
high-perfomance implementations of the core features of the popular 
[VineCopula R library](https://github.com/tnagler/VineCopula), in particular 
inference algorithms for both vine copula and bivariate copula models. 
Advantages over VineCopula are  
* interfaces to both both R and python (coming soon)
* a sleaker and more modern API,
* shorter runtimes, especially in high dimensions,
* nonparametric and multi-parameter families.

#### Status
 
The library is under active development and the first release (0.0.0.0) was
on March 24, 2017. The API is still rather unstable. We are also 
working on interfaces for R and python.


# Documentation

Below, we give a brief overview of the most important functionality. The full 
set of classes and methods can be found in the 
[Doxygen documenation](https://tvatter.github.io/vinecopulib/) along with 
minimal documentation.

- [Getting started](#getting-started)
	- [Requirements](#requirements)
	- [How to build the library](#how-to-build-the-library)
- [Bivariate copula models](#bivariate-copula-models)
	- [Set up a custom bivariate copula model](#set-up-a-custom-bivariate-copula-model)
	- [Implemented bivariate copula families](#implemented-bivariate-copula-families)
	- [Fit and select a bivariate copula](#fit-and-select-a-bivariate-copula)
	- [Work with a bivariate copula model](#work-with-a-bivariate-copula-model)
- [Vine copula models](#vine-copula-models)
	- [Set up a custom vine copula model](#set-up-a-custom-vine-copula-model)
	- [How to read the R-vine matrix](#how-to-read-the-r-vine-matrix)
	- [Fit and select a vine copula model](#fit-and-select-a-vine-copula-model)
	- [Work with a vine copula model](#work-with-a-vine-copula-model)


------------------------------------------------

## Getting started

### Requirements

To build the library, you'll need:

   * [a C++11-compatible    compiler](https://en.wikipedia.org/wiki/List_of_compilers#C.2B.2B_compilers),
   * [CMake](https://cmake.org/),
   * [Boost 1.63](http://www.boost.org/),
   * [Eigen 3.3](http://eigen.tuxfamily.org/index.php?title=Main_Page),
   * [NLopt](http://ab-initio.mit.edu/wiki/index.php/NLopt).
   
Optionally, you'll need:
   * [Doxygen](http://www.stack.nl/~dimitri/doxygen/) (to build the documentations),
   * [R](https://www.r-project.org/about.html) and [VineCopula](https://github.com/tnagler/VineCopula) (to build the unit tests).
     
Note that a `findR.cmake` looks for R and VineCopula in the default locations 
for linux and osx, but problems might occur with versions installed from 
R/RStudio. Therefore, prior to building the library, it is recommended to use:

`sudo Rscript -e 'install.packages(c("VineCopula"), lib="/usr/lib/R/library", 
repos="http://cran.rstudio.com/")'`


### How to build the library

The one liner (from the root folder):

`mkdir build && cd build && cmake .. && make && make doc && 
sudo make install && bin/test_all`

| Step | Shell command  |
|-----------------------|------------------------------------|
| Create a build folder  | `mkdir build` |
| Move to the created folder  | `cd build` |
|  Create the `MakeFile` via cmake  |  `cmake .. ` (for the `Release` version) or `cmake .. -DCMAKE_BUILD_TYPE=Debug` (for the `Debug` mode)  |
|  Create the `MakeFile` via cmake while avoiding compiling the unit tests  | `cmake .. -DBUILD_TESTING=OFF`|
| Compile the library | `make` or `make -j n` where `n` is the number of cores |
| Build the documentation (optional)  | `make doc` |
| Install the library on linux/OSX (optional)  | `sudo make install` |
| Run unit tests (optional)  |  `bin/[test_executable]` |

------------------------------------------------

### How to include the library in other projects

Using `make install`, vinecopulib is installed in the usual location of the 
system, namely

* `<prefix>/lib/` (for the shared library),
* `<prefix>/include/` (for the headers),
* `<prefix>/lib/cmake/vinecopulib` (to allow cmake to find the library 
with `find_package`),

where `<prefix>` is e.g. `/usr/` or `/usr/local`. Note that
`make install` only copies `vinecopulib.hpp` in `<prefix>/include/` and 
puts the other headers in a subfolder `<prefix>/include/vinecopulib`, but using 
`#include <vinecopulib.hpp>` is enough to load both bivariate and vine functions.

The easiest way to include vinecopulib in another project (and to avoid writing makefiles) 
is to use CMake. For instance, an example projet where the source code to be linked could contain
* a `CMakeLists.txt` file for the project's setup,
* a subfolder `src` for the source code, containing
    * the source code,
    * another `CMakeLists.txt` file for the project libraries and executables.
 
The top-level `CMakeLists.txt` could be:
```cmake
cmake_minimum_required(VERSION 3.2)

set(CMAKE_CXX_STANDARD 11)

project(Example)

# Setting default folders
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

# C++ compile flags
set(CMAKE_CXX_FLAGS "-std=gnu++11 -Wextra -Wall -Wno-delete-non-virtual-dtor -Werror=return-type -O2 -DNDEBUG")

# Find vinecopulib package and dependencies
find_package(vinecopulib                  REQUIRED)
find_package(Boost 1.63                   REQUIRED)
include(cmake/findEigen3.cmake            REQUIRED)
include(cmake/findNlopt.cmake             REQUIRED)

# Set required variables for includes and libraries
set(external_includes ${VINECOPULIB_INCLUDE_DIR} ${EIGEN3_INCLUDE_DIR} ${NLOPT_INCLUDE_DIR} ${Boost_INCLUDE_DIRS})
set(external_libs ${VINECOPULIB_LIBRARIES} ${NLOPT_LIBRARIES} ${Boost_LIBRARIES})

# Include subdirectory with project sources
add_subdirectory(src)
```
Assuming a single `main.cpp` source file (with `#include <vinecopulib.hpp>` at 
the top), the `CMakeLists.txt` file in `/src/` 
could then be:
```cmake
# Include header files
include_directories(${external_includes})

# Add main executable
add_executable(main main.cpp)

# Link to vinecopulib and dependencies
target_link_libraries(vinecopulib_main ${external_libs})
```

## Bivariate copula models

Bivariate copula models are implemented as the `Bicop` class, and `BicopFamily` 
is a closely related enum class describing the type or "family" of copula.

To use bivariate copula models in your code, include the header 
`vinecopulib/bicop/class.hpp` (or simply `vinecopulib.hpp`) at the top of 
your source.

### Implemented bivariate copula families

| type          | name                  | BicopFamily |
|---------------|-----------------------|-------|
| -             | Independence          | indep     |
| Elliptical    | Gaussian              | gaussian     |
| "             | Student t             | student     |
| Archimedean   | Clayton               | clayton     |
| "             | Gumbel                | gumbel     |
| "             | Frank                 | frank     |
| "             | Joe                   | joe     |  
| "             | BB1                   | bb1     |  
| "             | BB6                   | bb6     |  
| "             | BB7                   | bb7     |  
| "             | BB8                   | bb8    |  
| Nonparametric | Transformation kernel | tll0  | 

Note that several convenience vectors of families are included in the 
sub-namespace `bicop_families`:
* `all` contains all the families
* `parametric` contains the parametric families (all except `tll0`)
* `nonparametric` contains the nonparametric families (`indep` and `tll0`)
* `one_par` contains the parametric families with a single parameter 
(`gaussian`, `clayton`, `gumbel`, `frank`, and `joe`)
* `two_par` contains the parametric families with two parameters
(`student`, `bb1`, `bb6`, `bb7`, and `bb8`)
* `elliptical` contains the elliptical families
* `archimedean` contains the archimedean families
* `BB` contains the BB families
* `rotationless` contains the families that cover both positive and negative
association (`indep`,`gaussian`, `student`, `frank`, `tll0`)
* `lt` families that display lower-tail dependence (`clayton`, `bb1`, `bb7`)
* `ut` families that display upper-tail dependence (`gumbel`, `joe`, `bb1`, 
`bb6`, `bb7`, and `bb8`)
* `itau` families for which estimation by Kendall's tau inversion is available
(`indep`,`gaussian`, `student`,`clayton`, `gumbel`, `frank`, `joe`)
* `flip_by_rotation` families which can be flipped by rotation (currently 
equivalent to `archimedean`)

**Example**
``` cpp
// print all available families
std::cout << "Available families : ";
for (auto family : vinecopulib::bicop_families::all) {
		std::cout << get_family_name(family) << " ";
}
```

### Set up a custom bivariate copula model

There are essentially two ways of setting-up bivariate copulas:
* with known parameters,
* from data (i.e., with estimated parameters).

The constructor with known parameters takes 3 arguments:
* The copula family (default to `indep`)
* The rotation (default to `0`)
* The parameters (default to parameters corresponding to an independence copula)
**Example**
``` cpp
// 90 degree rotated Clayton with default parameter (corresponds to independence)
Bicop clayton(BicopFamily::clayton, 90);

// Gauss copula with parameter 0.5
Bicop gauss(BicopFamily::gaussian, 0,  VecXd::Constant(1, 0.5));
```
The constructor from data takes the same arguments as the select method described 
below.

### Fit and select a bivariate copula

You can either fit the parameters of a given `Bicop` object with `fit()` or
select the best fitting model from a set of families with `select()`. 

**Example**
``` cpp
// create a Gauss copula with parameter 0.5 and simulate 1e3 observations
Bicop model(BicopFamily::gaussian, 0,  VecXd::Constant(1, 0.5));             
auto data = model.simulate(1e3);

// instantiate a gaussian copula with default parameters and fit to data
Bicop fitted(BicopFamily::gaussian);
fitted.fit(data);
std::cout << 
    "estimated parameter: " <<
    fitted.get_parameters() << 
    "\n";
    
// assign another family to the same variable and fit to data
fitted = Bicop(BicopFamily::student);
fitted.fit(data);
std::cout << 
    "estimated parameter: " <<
    fitted.get_parameters() << 
    "\n";

// alternatively, assign to a family and fit automatically
fitted = select(data);
std::cout << 
    "family: " << fitted.get_family_name() <<
    "rotation: " <<  fitted.get_rotation() <<
    "\n";
```
As it's arguably the most important function of the `Bicop` class, it's worth 
understanding the arguments of `select()`:
* `Eigen::Matrix<double, Eigen::Dynamic, 2> data` is the data to fit from and 
it should be in (0,1).
* `std::vector<BicopFamily> family_set` describes the set of family to select 
from. It can take a user specified vector of 
families or any of those mentioned above (default is `bicop_families::all`).
* `std::string method` describes the estimation method. It can take `"mle"` 
(default, for maximum-likelihood estimation) and 
`"itau"` (for Kendall's tau inversion, although only available for families 
included in `bicop_families::itau`). Note that nonparametric families have 
specialized methods for which no specification is required.
* `std::string selection_criterion` describes the criterion to compare the 
families. It can take either `"bic"`(default) or `"aic"`.
* `bool preselect_families` describes a heuristic preselection method (default 
is `true`) based on symmetry properties of the data (e.g., the unrotated 
Clayton won't be preselected if the data displays upper-tail dependence). 

As mentioned above, the arguments of `select()` can be used as arguments to a
 constructor allowing to instantiate a new object directly:
``` cpp
// instantiate an archimedean copula by selecting the "best" family according to
// the BIC and parameters corresponding to the MLE
Bicop best_archimedean(data, bicop_families::archimedean);
std::cout << 
    "family: " << best_archimedean.get_family_name() <<
    "rotation: " <<  best_archimedean.get_rotation() <<
    best_archimedean.get_parameters() << 
    "\n"
    
// instantiate a bivariate copula by selecting the "best" family according to
// the AIC and parameters corresponding to Kendall's tau inversion
Bicop best_itau(data, bicop_families::itau, "itau", "aic");
std::cout << 
    "family: " << best_itau.get_family_name() <<
    "rotation: " <<  best_itau.get_rotation() <<
    best_itau.get_parameters() << 
    "\n"
```


### Work with a bivariate copula model

You can simulate from a model and evaluate the pdf, h-functions, inverse 
h-functions, log-likelihood, AIC, and BIC.

**Example**
``` cpp
// Gauss copula with parameter 0.5
Bicop bicop(BicopFamily::gaussian, 0,  VecXd::Constant(1, 0.5));

// Simulate 100 observations
auto sim_data = bicop.simulate(100);

// Evaluate the pdf
auto pdf  = bicop.pdf(sim_data);

// Evaluate the two h-functions
auto h1   = bicop.hfunc1(sim_data);
auto h2   = bicop.hfunc2(sim_data);

// Evalute the two inverse h-functions
auto hi1  = bicop.hinv1(sim_data);
auto hi2  = bicop.hinv2(sim_data);

// Evaluate the log-likelihood, AIC, and BIC
auto ll   = bicop.loglik(sim_data);
auto aic  = bicop.aic(sim_data);
auto bic  = bicop.bic(sim_data);
```

------------------------------------------------

## Vine copula models

Vine copula models are implemented in the class `Vinecop`. To use this class in
your code, include the header include the header `vinecopulib/vinecop/class.hpp` 
(or simply `vinecopulib.hpp`) at the top of your source. This automatically 
enables all features for bivariate copula models.

### Set up a custom vine copula model

Custom models can be created through the constructor of `Vinecop`. A model is
represented by a `std::vector<std::vector<Bicop>>` containing all 
pair-copulas and an [R-vine matrix](#how-to-read-an-r-vine-matrix).

Similarly to bivariate copulas, there are essentially two ways of setting-up 
vine copulas:
* with known parameters,
* from data (i.e., with estimated parameters).

The constructor with known parameters has two versions:
* one for which the only argument is the dimension, allowing to set-up a D-vine 
with only independence copulas,
* and one for which the two arguments are a matrix of integers (i.e., and 
[R-vine matrix](#how-to-read-an-r-vine-matrix)) and a 
`std::vector<std::vector<Bicop>>` containing all pair-copulas.

**Example**
``` cpp
// specify the dimension of the model
int d = 3;

// instantiate a three dimensional D-vine with independence copulas
Vinecop default_model(d);

// alternatively, instantiate a std::vector<std::vector<Bicop>> object
auto pair_copulas = Vinecop::make_pair_copula_store(3);  

// specify the pair copulas
for (int tree = 0; tree < d - 1; ++tree) {
    for (int edge = 0; edge < d - 1 - tree; ++edge) {
        // 90 degree Clayton with parameter 3.0
        pair_copulas[tree][edge] = Bicop(3, VecXd::Constant(1, 3.0), 90);
    }
}

// specify a structure matrix
Eigen::MatrixXi mat;
mat << 1, 1, 1,
       2, 2, 0,
       3, 0, 0;

// instantiate the custom model
Vinecop custom_model(pair_copulas, mat);
```
The constructors from data take the same arguments as the two select methods 
described [below](#fit-and-select-a-vine-copula-model).

### How to read the R-vine matrix

The R-vine matrix notation in vinecopulib different from the one in VineCopula.
An example matrix is
```
1 1 1 1
2 2 2 0
3 3 0 0
4 0 0 0
```
which encodes the following pair-copulas:

| tree | edge | pair-copulas   |
|------|------|----------------|
| 0    | 0    | `(4, 1)`       |
|      | 1    | `(3, 1)`       |
|      | 2    | `(2, 1)`       |
| 1    | 0    | `(4, 2; 1)`    |
|      | 1    | `(3, 2; 1)`    |
| 2    | 0    | `(4, 3; 2, 1)` |

Denoting by `M[i][j]` the matrix entry in row `i` and column `j` (starting at
0), the pair-copula index for edge `e` in tree `t` of a `d` dimensional vine is
`(M[d - 1 - t][e], M[t][e]; M[t - 1][e], ..., M[0][e])`. Less formally,
1. Start with the counter-diagonal element of column `e` (first conditioned 
   variable).
2. Jump up to the element in row `t` (second conditioned variable).
3. Gather all entries further up in column `e` (conditioning set).

### Fit and select a vine copula model

The method `select_all()` performs parameter estimation and automatic 
model selection when the vine structure is unknown, 
using the sequential procedure proposed by
 [Dissman et al. (2013)](https://mediatum.ub.tum.de/doc/1079277/1079277.pdf). 
 Alternatively, `select_families()`
The only mandatory argument is the data (stored in a `Eigen::MatrixXd`), 
additional arguments allow for customization of the fit 
(**TODO:** link to Doxygen page).

**Example**
``` cpp
MatXd data = simulate_uniform(100, 3);        // dummy data
Vinecop fit_default = Vinecop::select(data);  // fit and select model
```

### Work with a vine copula model

You can simulate from a vine copula model and evaluate its density function.

**Example**
``` cpp
Vinecop model(5);             // 5-dimensional independence vine
auto u = model.simulate(100)  // simulate 100 observations
auto c = model.pdf(u)         // evaluate the density on u
```
