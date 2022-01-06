# Find the main dependencies
find_package(Eigen3                       REQUIRED)
include(cmake/findR.cmake                 REQUIRED)
find_package(Boost 1.56                   REQUIRED)
find_package(Threads                      REQUIRED)

# Find wdm and install if not found
find_package(wdm QUIET)
if(NOT wdm_FOUND)
  message(STATUS "BLABLA")
  include(FetchContent)
  FetchContent_Declare(
    wdm
    URL https://github.com/tnagler/wdm/archive/refs/heads/master.zip
  )
  FetchContent_MakeAvailable(wdm)
  set(wdm_INCLUDE_DIRS "${wdm_SOURCE_DIR}/include")
endif()

# Set all the external dependencies
set(external_includes ${EIGEN3_INCLUDE_DIR} ${Boost_INCLUDE_DIRS} ${wdm_INCLUDE_DIRS})

# Find doxygen and configure if found
find_package(Doxygen QUIET)
if(DOXYGEN_FOUND)
    configure_file(
        ${CMAKE_CURRENT_SOURCE_DIR}/docs/Doxyfile.in
        ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY
    )
    add_custom_target(doc
        ${DOXYGEN_EXECUTABLE}
        ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        COMMENT "Generating API documentation with Doxygen" VERBATIM
    )
endif (DOXYGEN_FOUND)
