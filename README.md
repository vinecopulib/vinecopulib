# vinecopulib
[![Build Status](https://travis-ci.org/vinecopulib/vinecopulib.svg?branch=master)](https://travis-ci.org/vinecopulib/vinecopulib)
[![Windows Build status](http://ci.appveyor.com/api/projects/status/github/vinecopulib/vinecopulib?branch=master&svg=true)](https://ci.appveyor.com/project/vinecopulib/vinecopulib/branch/master)
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
[Eigen](http://eigen.tuxfamily.org/index.p  hp?title=Main_Page). It provides
high-performance implementations of the core features of the popular
[VineCopula R library](https://github.com/tnagler/VineCopula), in particular
inference algorithms for both vine copula and bivariate copula models.
Advantages over VineCopula are  
* a stand-alone C++ library with interfaces to both R and Python,
* a sleaker and more modern API,
* shorter runtimes and lower memory consumption, especially in high dimensions,
* nonparametric and multi-parameter families.

#### Status

Version [0.2.2](https://github.com/vinecopulib/vinecopulib/releases) was
released on November 9, 2017. While we did our best to
design a user-friendly API, the library is still under active development and
changes are to be expected. We are also working on interfaces for
[R](https://github.com/vinecopulib/rvinecopulib) and
[Python](https://github.com/vinecopulib/pyvinecopulib).

#### Contact

If you have any questions regarding the library, feel free to
[open an issue](https://github.com/vinecopulib/vinecopulib/issues/new) or
send a mail to <info@vinecopulib.org>.


# Documentation

Below, we give a brief overview of the most important functionality. The full
set of classes and methods can be found in the
[API documentation](https://vinecopulib.github.io/vinecopulib/).

- [Getting started](#getting-started)
	- [Requirements](#requirements)
	- [How to build the library](#how-to-build-the-library)
	- [How to include the library in other projects](#how-to-inlude-the-library-in-other-projects)
	- [Namespaces](#namespaces)
- [Bivariate copula models](#bivariate-copula-models)
	- [Implemented bivariate copula families](#implemented-bivariate-copula-families)
	- [Set up a custom bivariate copula model](#set-up-a-custom-bivariate-copula-model)
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

To build the library, you'll need at minimum:

   * [a C++11-compatible compiler (tested with GCC 6.3.0 and Clang 3.5.0 on Linux and AppleClang 8.0.0 on OSX)](https://en.wikipedia.org/wiki/List_of_compilers#C.2B.2B_compilers)
   * [CMake 3.2 (or later)](https://cmake.org/)
   * [Boost 1.56 (or later)](http://www.boost.org/)
   * [Eigen 3.3 (or later)](http://eigen.tuxfamily.org/index.php?title=Main_Page)

Optionally, you'll need:
   * [Doxygen](http://www.stack.nl/~dimitri/doxygen/) (to build the documentations)
   * [R](https://www.r-project.org/about.html) and [VineCopula](https://github.com/tnagler/VineCopula) (to run the unit tests)

Note that:
 
   * The [C++11 thread support library](http://en.cppreference.com/w/cpp/thread), 
   available along with any C++11 compiler on 
   OSX/Windows/most-linux-distributions, is used for multithreading.
   * A `findR.cmake` looks for R and VineCopula in the default locations for 
   linux and osx, but problems might occur with versions installed from
   R/RStudio. Therefore, prior to building the library, it is recommended to 
   use:
   
   `sudo Rscript -e 'install.packages(c("VineCopula"), lib="/usr/lib/R/library",
   repos="http://cran.rstudio.com/")'`

### How to build the library

By default, vinecopulib is header-only. It means that we use the CMake build 
system, but only to build the documentation and unit-tests, and to automate 
installation (i.e., place headers in the usual location). 
If you just want to use vinecopulib, you can use the header files 
(located in the`includes`folder) right away. 

The unix one liner (from the root folder):

`mkdir build && cd build && cmake .. && make && make doc &&
sudo make install && bin/test_all`

Alternatively, we provide an option to precompile compiled the library in 
order to save building time (and memory) via the CMake option 
`VINECOPULIB_SHARED_LIB`. In this case, source files are generated from header 
files and the CMake build system additionally allows to install the 
.dylib/.so/.dll object.

The unix one liner (from the root folder):

`mkdir build && cd build && cmake .. -DVINECOPULIB_SHARED_LIB=ON && make && 
make doc && sudo make install && bin/test_all`

| Step | Shell command  |
|-----------------------|------------------------------------|
| Create a build folder  | `mkdir build` |
| Move to the created folder  | `cd build` |
| Create the `MakeFile` via cmake  |  `cmake .. ` (or `cmake .. -DVINECOPULIB_SHARED_LIB=ON` for the compiled version)  |
| Compile the library | `make` or `make -j n` where `n` is the number of cores |
| Build the documentation (optional)  | `make doc` |
| Install the library on linux/OSX (optional)  | `sudo make install` |
| Run unit tests (optional)  |  `bin/[test_executable]` |

To install the library without unit tests, the `MakeFile` can be created via
 `cmake .. -DBUILD_TESTING=OFF`. Additionally, a `Debug` mode is available via 
 `cmake .. -DCMAKE_BUILD_TYPE=Debug`.

On Windows, CMake will generate Visual Studio files instead of Makefiles,
the following sequence of commands can be used to perform compilation using the command prompt:
```
md build
cd build
cmake ..
cmake --build . --config Debug
cmake --build . --config Release
cmake --build . --config Release --target install
```
Instead of the `cmake --build` commands, the generated `vinecopulib.sln` file can be open in the Visual Studio GUI. Furthermore, 
as for linux systems, the third line can be replaced by
 `cmake .. -DVINECOPULIB_SHARED_LIB=ON` to generate the source files in order 
 to compile vinecopulib in non-header-only mode.

The following CMake flags (given with example values) will likely come handy:
```
-DBOOST_ROOT=c:\local\boost_1_63_0
-DEIGEN3_INCLUDE_DIR=c:\local\eigen-eigen-da9b4e14c255
-DCMAKE_INSTALL_PREFIX=c:\local\vinecopulib-install
-DCMAKE_GENERATOR_PLATFORM=x64
-DBOOST_DEBUG=1
```

------------------------------------------------

### How to include the library in other projects

Using `make install`, vinecopulib is installed in the usual location of the
system, namely

* `<prefix>/include/` (for the headers),
* `<prefix>/lib/` (for the shared library when `VINECOPULIB_SHARED_LIB=ON` is used),
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

**Example**
```cmake
cmake_minimum_required(VERSION 3.2)

set(CMAKE_CXX_STANDARD 11)

project (Example)

# Setting default folders
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

# C++ compile flags
if (NOT WIN32)
  set(CMAKE_CXX_FLAGS "-std=gnu++11 -Wextra -Wall -Wno-delete-non-virtual-dtor -Werror=return-type -O2 -DNDEBUG")
endif()

# Find vinecopulib package and dependencies
find_package(vinecopulib                  REQUIRED)
find_package(Boost 1.56                   REQUIRED)
include(cmake/findEigen3.cmake            REQUIRED)
find_package(Threads                      REQUIRED)

# Set required variables for includes and libraries
# In the second line
#   * VINECOPULIB_LIBRARIES is needed if vinecopulib has been built as a
#     shared lib (does nothing otherwise).
#   * CMAKE_THREAD_LIBS_INIT is needed for some linux systems
#     (but does nothing on OSX/Windows).
set(external_includes ${VINECOPULIB_INCLUDE_DIR} ${EIGEN3_INCLUDE_DIR} ${Boost_INCLUDE_DIRS})
set(external_libs ${VINECOPULIB_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})

# Include subdirectory with project sources
add_subdirectory(src)
```
Assuming a single `main.cpp` source file (with `#include <vinecopulib.hpp>` at
the top), the `CMakeLists.txt` file in `/src/`
could then be:

**Example**
```cmake
# Include header files
include_directories(${external_includes})

# Add main executable
add_executable(main main.cpp)

# Link to vinecopulib if vinecopulib has been built as a shared lib
# and to pthreads on some linux systems (does nothing otherwise)
target_link_libraries(main ${VINECOPULIB_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
```

### Namespaces
In the examples mentioned above, it is assumed that `using namespace vinecopulib;`
is used. While the namespace `vinecopulib` contains the most important
functionalities described below, there are a two others that are available to the user:
-  `bicop_families`: convenience definitions of sets of bivariate copula families
-  `tools_stats`: utilities for statistical analysis

## Bivariate copula models

Bivariate copula models are implemented as the `Bicop` class, and `BicopFamily`
is a closely related enum class describing the type or "family" of copula.

To use bivariate copula models in your code, include the header
`vinecopulib/bicop/class.hpp` (or simply `vinecopulib.hpp`) at the top of
your source file.

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
| Nonparametric | Transformation kernel | tll  |

Note that several convenience vectors of families are included in the
sub-namespace `bicop_families`:
* `all` contains all the families
* `parametric` contains the parametric families (all except `tll`)
* `nonparametric` contains the nonparametric families (`indep` and `tll`)
* `one_par` contains the parametric families with a single parameter
(`gaussian`, `clayton`, `gumbel`, `frank`, and `joe`)
* `two_par` contains the parametric families with two parameters
(`student`, `bb1`, `bb6`, `bb7`, and `bb8`)
* `elliptical` contains the elliptical families
* `archimedean` contains the archimedean families
* `bb` contains the BB families
* `itau` families for which estimation by Kendall's tau inversion is available
(`indep`,`gaussian`, `student`,`clayton`, `gumbel`, `frank`, `joe`)


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
Bicop gauss(BicopFamily::gaussian, 0,  Eigen::VectorXd::Constant(1, 0.5));
```
The constructor from data takes the same arguments as the select method and is
described in the next section.

### Fit and select a bivariate copula

You can either fit the parameters of a given `Bicop` object with `fit()` or
select the best fitting model from a set of families with `select()`.

**Example**
``` cpp
// create a Gauss copula with parameter 0.5 and simulate 1e3 observations
Bicop model(BicopFamily::gaussian, 0,  Eigen::VectorXd::Constant(1, 0.5));             
auto data = model.simulate(1e3);

// instantiate a gaussian copula with default parameters and fit to data
Bicop fitted(BicopFamily::gaussian);
fitted.fit(data);
std::cout <<
    "estimated parameter: " <<
    fitted.get_parameters() <<
    std::endl;

// assign another family to the same variable and fit to data
fitted = Bicop(BicopFamily::student);
fitted.fit(data);
std::cout <<
    "estimated parameter: " <<
    fitted.get_parameters() <<
    std::endl;

// alternatively, assign to a family and fit automatically
fitted.select(data);
std::cout <<
    "family: " << fitted.get_family_name() <<
    "rotation: " <<  fitted.get_rotation() <<
    std::endl;
```
As it's arguably the most important function of the `Bicop` class, it's worth
understanding the second argument of `select()`, namely an object of the class
`FitControlsBicop`, which contain seven data members:
* `std::vector<BicopFamily> family_set` describes the set of family to select
from. It can take a user specified vector of
families or any of those mentioned above (default is `bicop_families::all`).
* `std::string parametric_method` describes the estimation method. It can take
  `"mle"` (default, for maximum-likelihood estimation) and
`"itau"` (for Kendall's tau inversion, although only available for families
included in `bicop_families::itau`).
* `std::string nonparametric_method` describes the degree of the density
approximation for the transformation kernel estimator. It can take
`constant`, `linear` and `quadratic` (default) for approximations of
degree zero, one and two.
* `double nonparametric_mult` a factor with which the smoothing parameters
are multiplied.
* `std::string selection_criterion` describes the criterion to compare the
families. It can take either `"loglik"`, `"aic"`, or `"bic"`(default).
* `bool preselect_families` describes a heuristic preselection method (default
is `true`) based on symmetry properties of the data (e.g., the unrotated
Clayton won't be preselected if the data displays upper-tail dependence).
* `size_t num_threads` number of threads to run in parallel when fitting
several families.

As mentioned [above](#set-up-a-custom-bivariate-copula-model), the arguments
of `select()` can be used as arguments to a
 constructor allowing to instantiate a new object directly:

 **Example**
``` cpp
// instantiate an archimedean copula by selecting the "best" family according to
// the BIC and parameters corresponding to the MLE
Bicop best_archimedean(data, FitControlsBicop(bicop_families::archimedean));
std::cout <<
    "family: " << best_archimedean.get_family_name() <<
    "rotation: " <<  best_archimedean.get_rotation() <<
    best_archimedean.get_parameters() <<
    std::endl

// instantiate a bivariate copula by selecting the "best" family according to
// the AIC and parameters corresponding to Kendall's tau inversion
FitControlsBicop controls(bicop_families::itau, "itau");
controls.set_selection_criterion("aic");
Bicop best_itau(data, controls));
std::cout <<
    "family: " << best_itau.get_family_name() <<
    "rotation: " <<  best_itau.get_rotation() <<
    best_itau.get_parameters() <<
    std::endl
```


### Work with a bivariate copula model

You can simulate from a model and evaluate the pdf, h-functions, inverse
h-functions, log-likelihood, AIC, and BIC.

**Example**
``` cpp
// Gauss copula with parameter 0.5
Bicop bicop(BicopFamily::gaussian, 0,  Eigen::VectorXd::Constant(1, 0.5));

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

Bivariate copula models can also be written to and constructed from JSON files
and `boost::property_tree::ptree` objects:

```
// Gauss copula with parameter 0.5
Bicop bicop(BicopFamily::gaussian, 0,  Eigen::VectorXd::Constant(1, 0.5));

// Save as a ptree object
boost::property_tree::ptree bicop_node = bicop.to_ptree();

// Write into a JSON file
boost::property_tree::write_json("myfile.JSON", bicop_node);

// Equivalently
bicop.to_json("myfile.JSON");

// Then a new Bicop can be constructed from the ptree object
Bicop bicop2(bicop_node);

// Or from the JSON file
Bicop bicop3("myfile.JSON");
```

------------------------------------------------

## Vine copula models

Vine copula models are implemented in the class `Vinecop`. To use this class in
your code, include the header include the header `vinecopulib/vinecop/class.hpp`
(or simply `vinecopulib.hpp`) at the top of your source file. This automatically
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
auto pair_copulas = Vinecop::make_pair_copula_store(d);  

// specify the pair copulas
auto par = Eigen::VectorXd::Constant(1, 3.0);
for (auto& tree : pair_copulas) {
    for (auto& pc : tree) {
        pc = Bicop(BicopFamily::clayton, 270, par);
    }
}

// specify a structure matrix
Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic> mat(3, 3);
mat << 1, 1, 1,
       2, 2, 0,
       3, 0, 0;

// instantiate the custom model
Vinecop custom_model(pair_copulas, mat);
```
The constructors from data take the same arguments as the two select methods
described [below](#fit-and-select-a-vine-copula-model).

### How to read the R-vine matrix

The R-vine matrix notation in vinecopulib is different from the one in VineCopula.
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

Denoting by `M[i, j]` the matrix entry in row `i` and column `j`, the pair-copula index for edge `e` in tree `t` of a `d` dimensional vine is
`(M[d - 1 - t, e], M[t, e]; M[t - 1, e], ..., M[0, e])`. Less formally,
1. Start with the counter-diagonal element of column `e` (first conditioned
   variable).
2. Jump up to the element in row `t` (second conditioned variable).
3. Gather all entries further up in column `e` (conditioning set).

A valid R-vine matrix must satisfy several conditions which are checked
when `RVineMatrix()` is called:
1. The lower right triangle must only contain zeros.
2. The upper left triangle can only contain numbers between 1 and d.
3. The antidiagonal must contain the numbers 1, ..., d.
4. The antidiagonal entry of a column must not be contained in any
   column further to the right.
5. The entries of a column must be contained in all columns to the left.
6. The proximity condition must hold: For all t = 1, ..., d - 2 and
   e = 0, ..., d - t - 1 there must exist an index j > d, such that
   `(M[t, e], {M[0, e], ..., M[t-1, e]})` equals either
   `(M[d-j-1, j], {M[0, j], ..., M[t-1, j]})` or
   `(M[t-1, j], {M[d-j-1, j], M[0, j], ..., M[t-2, j]})`.

Condition 6 already implies conditions 2-5, but is more difficult to
check by hand.

### Fit and select a vine copula model

The method `select_all()` performs parameter estimation and automatic
model selection when the vine structure is unknown (i.e., it modifies the
structure of the object), using the sequential procedure proposed by
 [Dissman et al. (2013)](https://mediatum.ub.tum.de/doc/1079277/1079277.pdf).
 Alternatively, `select_families()` performs parameter estimation and automatic
 model selection for a known structure (i.e., using the structure of the
 object). In both cases, the only mandatory argument is the data
 (stored in a `Eigen::MatrixXd`), while controls argument allow for
 customization of the fit.

**Example**
``` cpp
// specify the dimension of the model
int d = 5;

// simulate dummy data
Eigen::MatrixXd data = tools_stats::simulate_uniform(100, d);

// instantiate a D-vine and select the families
Vinecop model(d);
model.select_families(data);

// alternatively, select the structure along with the families
model.select_all(data);
```

Note that the second argument to `select_all()` and `select_families()` is
similar to the one of `select()` for `Bicop` objects. Objects of the class
`FitControlsVinecop` inherit from `FitControlsBicop` and extend them with 
additional data members to control the structure selection:
* `int truncation_level` describes the tree after which `family_set` is set to
`{BicopFamily::indep}`. In other words, all pair copulas in trees lower than
`truncation_level` (default to none) are "selected" as independence copulas.
* `std::string tree_criterion` describes the criterion used to construct the
minimum spanning tree (see
[Dissman et al. (2013)](https://mediatum.ub.tum.de/doc/1079277/1079277.pdf)).
It can take `"tau"` (default) for Kendall's tau, `"rho"` for Spearman's rho,
or `"hoeffd"` for Hoeffding's D (suited for non-monotonic relationships).
* `double threshold` describes a value (default is 0) of `tree_criterion` under
which the corresponding pair-copula is set to independence.
* `bool select_truncation_level` can be set to true to select the truncation
level automatically (default is `false`).
* `bool select_threshold` can be set to true to select the threshold parameter
automatically (default is `false`).
* `size_t num_threads` number of threads to run in parallel when fitting pair 
copulas within one tree.

As mentioned [above](#set-up-a-custom-vine-copula-model), the arguments
of `select_all()` and `select_families()` can be used as arguments to a
 constructor allowing to instantiate a new object directly:

 **Example**
``` cpp
// specify the dimension of the model
int d = 4;

// simulate dummy data
Eigen::MatrixXd data = simulate_uniform(100, d);

// instantiate a vine from data using the default arguments
Vinecop best_vine(data);

// alternatively, instantiate a structure matrix...
Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic> M;
M << 1, 1, 1, 1,
     2, 2, 2, 0,
     3, 3, 0, 0
     4, 0, 0, 0;

// ... and instantiate a vine copula from data using the custom structure,
// Kendall's tau inversion for parameters
// estimation and a truncation after the second tree
FitControlsVinecop controls(bicop_families::itau, "itau");
controls.set_truncation_level(2);
controls.set_num_threads(4);  // parallelize with 4 threads
Vinecop custom_vine(data, M, controls);
```

### Work with a vine copula model

You can simulate from a vine copula model, evaluate its density, distribution,
log-likelihood, AIC and BIC.

**Example**
``` cpp
// 5-dimensional independence vine
Vinecop model(5);           

// simulate 100 observations
auto data = model.simulate(100)  

// evaluate the density
auto pdf = model.pdf(data)

// evaluate the distribution
auto cdf = model.cdf(data)

// evaluate the log-likelihood
auto ll = model.loglik(data)

// evaluate the AIC
auto aic = model.aic(data)

// evaluate the BIC
auto bic = model.bic(data)
```

Vine copula models can also be written to and constructed from JSON files
and `boost::property_tree::ptree` objects:

```
// 5-dimensional vine copula
Vinecop vinecop(5);

// Save as a ptree object
boost::property_tree::ptree vinecop_node = vinecop.to_ptree();

// Write into a JSON file
boost::property_tree::write_json("myfile.JSON", vinecop_node);

// Equivalently
vinecop.to_json("myfile.JSON");

// Then a new Bicop can be constructed from the ptree object
Vinecop vinecop2(vinecop_node);

// Or from the JSON file
Vinecop vinecop2("myfile.JSON");
```
