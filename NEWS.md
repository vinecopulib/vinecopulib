# vinecopulib 0.0.2 (March 31, 2017)

MAJOR CHANGES 

    * all `tools_xxx` namespaces are no sub-namespaces of `vinecopulib` (#130).

    * header files are encapsulate in an addtional `vinecopulib/` folder, i.e.,
      `include/vinecopulib/subdir/file.hpp` (#126).

    * removed abitility to extract the git revision (#124).

    * new header `misc/tools_interface.hpp` where interface-specific behavior
      can be defined (for exaple, a custom version of `std::cout`) (#136).

BUG FIXES

    * fix `mat.array() = 0` error on some compilers (#131).

    * add missing `<exception>` header (caused errors on not fully C++11 
      compliant compilers, #139).
          
# vinecopulib 0.0.1 (March 29, 2017)

Initial rlease.
