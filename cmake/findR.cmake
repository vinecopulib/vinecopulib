#
# - Locates an installed R distribution and configures the parity-test header.
#
# Defines the following:
#  RSCRIPT_EXECUTABLE  - Path to the Rscript command, if found
#

if(NOT BUILD_TESTING)
    return()
endif()

set(TEMP_CMAKE_FIND_APPBUNDLE ${CMAKE_FIND_APPBUNDLE})
set(CMAKE_FIND_APPBUNDLE "NEVER")
find_program(RSCRIPT_EXECUTABLE Rscript DOC "Rscript executable.")
set(CMAKE_FIND_APPBUNDLE ${TEMP_CMAKE_FIND_APPBUNDLE})

# option() stores a BOOL, so the default has to be ON/OFF rather than a path.
if(RSCRIPT_EXECUTABLE)
    set(_r_parity_default ON)
else()
    set(_r_parity_default OFF)
endif()
option(VINECOPULIB_R_PARITY_TESTS
       "Cross-check results against the VineCopula R package"
       ${_r_parity_default})

if(VINECOPULIB_R_PARITY_TESTS AND NOT RSCRIPT_EXECUTABLE)
    message(FATAL_ERROR
            "VINECOPULIB_R_PARITY_TESTS is ON but Rscript is not on PATH. The "
            "parity tests need R with VineCopula >= 2.6.2 (from GitHub, not "
            "CRAN); configure with -DVINECOPULIB_R_PARITY_TESTS=OFF to skip "
            "them.")
endif()

if(VINECOPULIB_R_PARITY_TESTS)
    set(VINECOPULIB_HAS_RSCRIPT 1)
else()
    set(VINECOPULIB_HAS_RSCRIPT 0)
endif()

# Always generated: the test sources include it unconditionally and branch on
# VINECOPULIB_HAS_RSCRIPT.
set(r_script_dir "${PROJECT_BINARY_DIR}/test")
configure_file(${CMAKE_CURRENT_LIST_DIR}/templates/rscript.hpp.in
               ${CMAKE_BINARY_DIR}/generated/test/rscript.hpp @ONLY)
