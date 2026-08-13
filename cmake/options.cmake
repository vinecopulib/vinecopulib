# Developer-facing switches apply to a top-level build only, so that consumers
# using add_subdirectory or FetchContent keep their own settings.
if(NOT DEFINED VINECOPULIB_IS_TOP_LEVEL)
  if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    set(VINECOPULIB_IS_TOP_LEVEL ON)
  else()
    set(VINECOPULIB_IS_TOP_LEVEL OFF)
  endif()
endif()

if(VINECOPULIB_IS_TOP_LEVEL)
  set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
endif()

# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build (Debug, Release, RelWithDebInfo, MinSizeRel)" FORCE)
endif()

option(VINECOPULIB_PRECOMPILED   "Pre-compiled version"              OFF)
option(VINECOPULIB_SANITIZERS    "Address/UB sanitizers (Debug)"     OFF)
option(BUILD_TESTING             "Build tests."                      ON)
option(VINECOPULIB_CODE_COVERAGE "Code coverage."                    OFF)
option(VINECOPULIB_STRICT_COMPILER "Stricter compiler warnings"      OFF)
option(VINECOPULIB_BUILD_DOC     "Build documentation"               OFF)
option(VINECOPULIB_BUILD_BENCHMARKS "Build benchmark suite"          OFF)
option(VINECOPULIB_NATIVE_ARCH   "Optimize for the build machine's CPU (not redistributable)" OFF)

# The unprefixed names, honoured for one release.
foreach(_old_new "OPT_ASAN;VINECOPULIB_SANITIZERS"
                 "CODE_COVERAGE;VINECOPULIB_CODE_COVERAGE"
                 "STRICT_COMPILER;VINECOPULIB_STRICT_COMPILER"
                 "BUILD_DOC;VINECOPULIB_BUILD_DOC")
  list(GET _old_new 0 _old)
  list(GET _old_new 1 _new)
  if(DEFINED ${_old})
    message(DEPRECATION "${_old} is deprecated; use ${_new}.")
    set(${_new} ${${_old}})
  endif()
endforeach()
