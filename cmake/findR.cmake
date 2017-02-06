#
# - This module locates an installed R distribution.
#
# Defines the following:
#  R_COMMAND                - Path to R command
#  R_HOME                   - Path to 'R home', as reported by R
#  R_INCLUDE_DIR            - Path to R include directory
#  R_LIBRARY_BASE           - Path to R library
#  R_LIBRARY_BLAS           - Path to Rblas / blas library
#  R_LIBRARY_LAPACK         - Path to Rlapack / lapack library
#  R_LIBRARY_READLINE       - Path to readline library
#  R_LIBRARIES              - Array of: R_LIBRARY_BASE, R_LIBRARY_BLAS, R_LIBRARY_LAPACK, R_LIBRARY_BASE [, R_LIBRARY_READLINE]
#  Rcpp_INCLUDE_DIR         - Path to Rcpp include directory
#  Rcpp_LIBRAIRIES          - Path to the Rcpp library
#  RInside_INCLUDE_DIR      - Path to RInside include directory
#  RcppEigen_INCLUDE_DIR    - Path to RcppEigen include directory
#
#  VTK_R_HOME          - (deprecated, use R_HOME instead) Path to 'R home', as reported by R
#

set(TEMP_CMAKE_FIND_APPBUNDLE ${CMAKE_FIND_APPBUNDLE})
set(CMAKE_FIND_APPBUNDLE "NEVER")
find_program(R_COMMAND R DOC "R executable.")
set(CMAKE_FIND_APPBUNDLE ${TEMP_CMAKE_FIND_APPBUNDLE})

if(NOT R_LIB_ARCH)
    set(R_LIB_ARCH "ThisPathNotExists")
endif()

if(R_COMMAND)
    execute_process(WORKING_DIRECTORY .
            COMMAND ${R_COMMAND} RHOME
            OUTPUT_VARIABLE R_ROOT_DIR
            OUTPUT_STRIP_TRAILING_WHITESPACE)
    # deprecated
    if(VTK_R_HOME)
        set(R_HOME ${VTK_R_HOME} CACHE PATH "R home directory obtained from R RHOME")
    else(VTK_R_HOME)
        set(R_HOME ${R_ROOT_DIR} CACHE PATH "R home directory obtained from R RHOME")
        set(VTK_R_HOME ${R_HOME})
    endif(VTK_R_HOME)
    # /deprecated
    # the following command does nothing currently, but will be used when deprecated code is removed
    set(R_HOME ${R_ROOT_DIR} CACHE PATH "R home directory obtained from R RHOME")

    find_path(R_INCLUDE_DIR R.h
            HINTS ${R_ROOT_DIR}
            PATHS /usr/local/lib /usr/local/lib64 /usr/share
            PATH_SUFFIXES include R/include
            DOC "Path to file R.h")
    if(NOT R_INCLUDE_DIR)
        message(FATAL_ERROR "R.h file not found.")
    endif()

    file(GLOB_RECURSE Rcpp_INCLUDE_DIR "${R_HOME}/*/Rcpp.h")
    if(NOT Rcpp_INCLUDE_DIR)
        file(GLOB_RECURSE Rcpp_INCLUDE_DIR "/usr/local/lib/R/*/Rcpp.h")
        if(NOT Rcpp_INCLUDE_DIR)
            message(FATAL_ERROR "Rcpp.h file not found.")
        endif()
    endif()
    string(REGEX REPLACE "/Rcpp.h$" "" Rcpp_INCLUDE_DIR "${Rcpp_INCLUDE_DIR}")

    file(GLOB_RECURSE RInside_INCLUDE_DIR "${R_HOME}/*/RInside.h")
    if(NOT RInside_INCLUDE_DIR)
        file(GLOB_RECURSE RInside_INCLUDE_DIR "/usr/local/lib/R/*/RInside.h")
        if(NOT RInside_INCLUDE_DIR)
            message(FATAL_ERROR "RInside.h file not found.")
        endif()
    endif()
    string(REGEX REPLACE "/RInside.h$" "" RInside_INCLUDE_DIR "${RInside_INCLUDE_DIR}")

    file(GLOB_RECURSE RcppEigen_INCLUDE_DIR "${R_HOME}/*/RcppEigen.h")
    if(NOT RcppEigen_INCLUDE_DIR)
        file(GLOB_RECURSE RcppEigen_INCLUDE_DIR "/usr/local/lib/R/*/RcppEigen.h")
        if(NOT RInside_INCLUDE_DIR)
            message(FATAL_ERROR "RcppEigen.h file not found.")
        endif()
    endif()
    string(REGEX REPLACE "/RcppEigen.h$" "" RcppEigen_INCLUDE_DIR "${RcppEigen_INCLUDE_DIR}")
    file(GLOB_RECURSE RInside_LIBRAIRIES "${RInside_INCLUDE_DIR}/../*.a")


    find_library(R_LIBRARY_BASE R
            HINTS ${R_ROOT_DIR}/lib ${R_ROOT_DIR}/bin/${R_LIB_ARCH}
            DOC "R library (example libR.a, libR.dylib, etc.).")
    if(NOT R_LIBRARY_BASE)
        set(R_LIBRARY_BASE "")
        message(
                "R library not found. Locations tried:\n"
                "${R_ROOT_DIR}/lib ${R_ROOT_DIR}/bin/${R_LIB_ARCH}"
        )
        if(APPLE)
            message(FATAL_ERROR "")
        endif()
    endif()

    find_library(R_LIBRARY_BLAS NAMES Rblas blas
            HINTS ${R_ROOT_DIR}/lib ${R_ROOT_DIR}/bin/${R_LIB_ARCH}
            DOC "Rblas library (example libRblas.a, libRblas.dylib, etc.).")

    find_library(R_LIBRARY_LAPACK NAMES Rlapack lapack
            HINTS ${R_ROOT_DIR}/lib ${R_ROOT_DIR}/bin/${R_LIB_ARCH}
            DOC "Rlapack library (example libRlapack.a, libRlapack.dylib, etc.).")
#<TEMP>
if (WIN32)
  set(R_LIBRARY_BLAS "")
  set(R_LIBRARY_LAPACK "")
endif()
#</TEMP>

    find_library(R_LIBRARY_READLINE readline
            DOC "(Optional) system readline library. Only required if the R libraries were built with readline support.")

else(R_COMMAND)
    message(SEND_ERROR "FindR.cmake requires the following variables to be set: R_COMMAND")
endif(R_COMMAND)

# Note: R_LIBRARY_BASE is added to R_LIBRARIES twice; this may be due to circular linking dependencies; needs further investigation
set(R_LIBRARIES ${R_LIBRARY_BASE} ${R_LIBRARY_BLAS} ${R_LIBRARY_LAPACK} ${R_LIBRARY_BASE})
if(R_LIBRARY_READLINE)
    set(R_LIBRARIES ${R_LIBRARIES} ${R_LIBRARY_READLINE})
endif(R_LIBRARY_READLINE)
