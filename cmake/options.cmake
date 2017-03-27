set(CMAKE_EXPORT_COMPILE_COMMANDS "ON")
set(BUILD_SHARED_LIBS ON)
set(CMAKE_MACOSX_RPATH 1)

option(WARNINGS_AS_ERRORS   "Compiler warnings as errors"       "OFF")
option(OPT_ASAN             "Use adress sanitizer (debug)"      "ON")
option(BUILD_TESTING        "Build tests."                      "ON")