#
# - This module locates an installed R distribution.
#
# Defines the following:
#  RSCRIPT_EXECUTABLE       - Path to the Rscript command
#

if(BUILD_TESTING)

    set(TEMP_CMAKE_FIND_APPBUNDLE ${CMAKE_FIND_APPBUNDLE})
    set(CMAKE_FIND_APPBUNDLE "NEVER")
    find_program(RSCRIPT_EXECUTABLE Rscript DOC "Rscript executable.")
    set(CMAKE_FIND_APPBUNDLE ${TEMP_CMAKE_FIND_APPBUNDLE})

    if(RSCRIPT_EXECUTABLE)
        configure_file(${CMAKE_CURRENT_LIST_DIR}/templates/rscript.hpp.in ${CMAKE_BINARY_DIR}/generated/test/rscript.hpp)
    else()
        # Fail here rather than at build time with a missing rscript.hpp.
        message(FATAL_ERROR
                "BUILD_TESTING requires Rscript on PATH (the R parity tests shell "
                "out to it). Configure with -DBUILD_TESTING=OFF to build without "
                "the test suite.")
    endif()

endif(BUILD_TESTING)
