# vinecopulib 0.2.2 (November 9, 2017)

NEW FEATURES

   * allow `"loglik"` as selection criterion (#267).

   
BUG FIXES
 
   * make interpolation grid symmetric around (0.5, 0.5) again (#268).


# vinecopulib 0.2.1 (November 7, 2017)

NEW FEATURES

   * faster vine copula estimation and selection by parallelizing further 
     sub-routines (#259).

   * enhanced cross-platform compatibility and addition of a `STRICT_COMPILER`
     option for gcc (#261).
        
   * increased precision of maximum-likelihood estimators (#264).
   
BUG FIXES

   * made arguments of pairwise dependence measures consistent (#258).
   
   * fixed `itau` estimation method for Frank copulas (only allowed for positive
     parameters) (#263).


# vinecopulib 0.2.0 (October 30, 2017)
 
LIBRARY TYPE

   * library is now header only by default (#246), with an option to compile it
     as a shared library (#249).

DEPENDENCIES

   * removed dependency on `NLopt` (#239).

NEW FEATURES
 
   * NA handling (#237, #238).
   
   * parallelized selection/estimation of (pair-) copulas (#240).
   
   * efficient storage and fitting of truncated vines (#248).
   
   * Brent line search for (profile-) maximum-likelihood estimation of 
     one-parameter families (#255).
     
   * more restrictive parameter bounds for Archimedean families, ensuring 
     their numerical stability (#256).

BUG FIXES

   * error thrown whenever `Vinecop` or `Bicop` constructors are called with
    datasets containing a single row (#251).
     
   * ensure `const` correctness of `Vinecop` and `Bicop` member functions 
     (#225).
   
   * made order of inverse Rosenblatt consistent for d = 2 and d > 2 (#232).
   
   * fixed bug in interpolation near upper right corner (#233).
   
   * interpolation grid is now symmetric around (0.5, 0.5) (#234).
   
   * stabilized quadratic tll estimator near zero (#235).
   
   * stabilized Archimedean pdfs (#256).

   
# vinecopulib 0.1.0 (August 23, 2017)

NEW FEATURES
 
   * read/write `Bicop` and `Vinecop` objects (#205) using 
     `boost::property_tree::ptree` with `to_ptree()`, `to_json()`, and 
     constructors taking `const char *filename` or a 
     `boost::property_tree::ptree` for both classes.
     
   * sparse selection of vine copulas (#206) using new data members in 
     `FitControlsVinecop`:
        * `bool select_truncation_level` whether the truncation is selected 
           automatically.
        * `bool select_threshold` whether the threshold parameter is selected 
           automatically.
        * `double threshold` sets a fixed threshold parameter.
        
   * local likelihood estimators (#216) have been implemented by refactoring the 
     `tll0` family into a more general `tll` family, where approximations of 
     degrees zero, one and two can be fitted by setting the new 
     `std::string nonparametric_method` data member of `FitControlsBicop` 
     respectively as `constant`, `linear` and `quadratic` (default).
     
   * Kendall's tau (#211) and normalization (#215) for kernel estimators.
     
   * support for clang compiler on linux (#201, #202, #203).
   
   * option to omit R-vine matrix check in `Vinecop` constructors (#198).
     
BUG FIXES

   * replacing throw `std::string` with throw `std::runtime_error` in 
    `tools_opimization.cpp` (#204).
   
   * ensure valid starting parameters in `Bicop::fit()` (#209, #210).
   
   * fix appveyor and travis problems (#208, #212, #213).

# vinecopulib 0.0.3 (June 7, 2017)

NEW FEATURES

   * new functions `Bicop::cdf()` and `Vinecop::cdf()` for evaluating the
     cumulative distribution function of bivariate and vine copulas (#177,
     #189).

   * the constructor of the `RVineMatrix` class now checks whether it is
     provided with a valid R-vine matrix (#192).

   * extended documentation to build the library under Windows (#188).

   * extended continuous integration tests for Windows (#150, #169).


BUG FIXES

   * vinecopulib.dll is installed to `lib/` instead of `bin/` (#149).

   * more pleasing and portable formatting of error messages (#147, #156, #159,
     #165).

   * fixed bugs in `Bicop::select()` caused by `0`s and `1`s or unsufficient
     data (#173, #180).

   * fixed compatibility issue with CMake 3.8 (#167).

   * fixed uninitialized memory issues on Windows (#169).



# vinecopulib 0.0.2 (March 31, 2017)

MAJOR CHANGES

  * all `tools_xxx` namespaces are no sub-namespaces of `vinecopulib` (#130).

  * header files are encapsulate in an addtional `vinecopulib/` folder, i.e.,
      `include/vinecopulib/subdir/file.hpp` (#126).

  * removed abitility to extract the git revision (#124).

  * new header `misc/tools_interface.hpp` where interface-specific behavior
      can be defined (for example, a custom version of `std::cout`) (#136).

BUG FIXES

   * fix `mat.array() = 0` error on some compilers (#131).

   * add missing `<exception>` header (caused errors on not fully C++11
      compliant compilers, #139).

# vinecopulib 0.0.1 (March 29, 2017)

Initial release.
