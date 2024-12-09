set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_MACOSX_RPATH 1)

# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build (Debug, Release, RelWithDebInfo, MinSizeRel)" FORCE)
endif()

option(VINECOPULIB_SHARED_LIB    "Pre-compiled version"              OFF)
option(BUILD_SHARED_LIBS         "shared/static lib"                 ON)
option(OPT_ASAN                  "Use adress sanitizer (debug)"      ON)
option(BUILD_TESTING             "Build tests."                      ON)
option(CODE_COVERAGE             "Code coverage."                    OFF)
option(STRICT_COMPILER           "Stricter compiler warnings"        OFF)
option(BUILD_DOC                 "Build documentation"               OFF)
