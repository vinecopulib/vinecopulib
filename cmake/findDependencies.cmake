# Find the main dependencies
message(STATUS "CMAKE_INSTALL_PREFIX: ${CMAKE_INSTALL_PREFIX}")
message(STATUS "BOOST_ROOT: ${Boost_ROOT}")
message(STATUS "Eigen3_DIR: ${Eigen3_DIR}")
message(STATUS 'BLABLABLA')
find_package(Eigen3                       REQUIRED)
include(cmake/findR.cmake                 REQUIRED)
find_package(Boost 1.56                   REQUIRED)
find_package(Threads                      REQUIRED)
find_package(wdm                          REQUIRED)

include(FetchContent)
FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG        6b74da4757a549563d7c37c8fae3e704662a043b # release-1.10.0
)

# Download googlestest
if(BUILD_TESTING)
  FetchContent_MakeAvailable(googletest)
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
