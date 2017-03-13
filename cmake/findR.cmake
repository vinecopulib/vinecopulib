#
# - This module locates an installed R distribution.
#
# Defines the following:
#  RSCRIPT_EXECUTABLE       - Path to the Rscript command
#

set(TEMP_CMAKE_FIND_APPBUNDLE ${CMAKE_FIND_APPBUNDLE})
set(CMAKE_FIND_APPBUNDLE "NEVER")
find_program(RSCRIPT_EXECUTABLE Rscript DOC "Rscript executable.")
set(CMAKE_FIND_APPBUNDLE ${TEMP_CMAKE_FIND_APPBUNDLE})

if(RSCRIPT_EXECUTABLE)
    configure_file(${CMAKE_CURRENT_LIST_DIR}/rscript.hpp.txt ${CMAKE_BINARY_DIR}/generated/rscript.hpp)
elseif(RSCRIPT_EXECUTABLE)
    message(SEND_ERROR "FindR.cmake requires the following variables to be set: RSCRIPT_EXECUTABLE")
endif(RSCRIPT_EXECUTABLE)

