# vinecopulib 0.0.3 (June 5, 2017)

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

Initial rlease.
